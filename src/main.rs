// Standard library imports for file I/O and path handling
use std::{path::PathBuf, sync::Arc};

// Command-line argument parsing library
use clap::Parser;

// Parallel processing
use rayon::prelude::*;

// Async runtime and utilities
use futures::future;

// Import our library modules
use varclock::{
    BedRegion, BgzOutput, analyze_read_allele_content_detailed, load_bed_regions,
    query_bam_records_for_region, query_vcf_variants_for_region, read_spans_variant_position,
};

/// Command-line arguments structure for the bioinformatics scanner tool
///
/// This tool analyzes genomic data by:
/// 1. Reading genomic regions from a BED file
/// 2. Finding genetic variants in those regions from a VCF file
/// 3. Checking which sequencing reads (from BAM file) contain those variants
/// 4. Outputting a summary report
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Input BED file containing genomic regions of interest
    /// BED format: chromosome, start position, end position (0-indexed, half-open intervals)
    #[arg(short, long)]
    bed: PathBuf,

    /// Input VCF file containing genetic variants
    /// VCF format: standardized format for representing genetic variations
    #[arg(short, long)]
    vcf: PathBuf,

    /// Input BAM file containing aligned sequencing reads
    /// BAM format: binary, compressed version of SAM (Sequence Alignment/Map)
    /// Note: Uses -a flag to avoid conflict with -b (bed file)
    #[arg(short = 'a', long)]
    bam: PathBuf,

    /// Output file path for results (will be in BGZ-compressed TSV format)
    /// Contains: read_id, timestamp, allele_match, variant_chrom, variant_pos, variant_ref, variant_alt, variant_description, variant_type, region, region_name, read_start, read_end, mapping_quality
    /// allele_match values: REFERENCE:sequence (read matches reference), VARIANT:sequence (read matches variant allele)
    /// File will be compressed using BGZF format and optionally indexed with tabix
    #[arg(short, long)]
    output: PathBuf,

    /// Create tabix index for compressed output (default: true)
    /// Enables fast random access to compressed genomic data
    #[arg(long, default_value_t = true)]
    create_index: bool,

    /// Number of parallel threads to use for processing (default: number of CPU cores)
    #[arg(short, long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Maximum number of concurrent regions to process simultaneously (default: 16)
    /// Higher values use more memory but may improve throughput for large datasets
    #[arg(long, default_value_t = 16)]
    max_concurrent_regions: usize,

    /// Chunk size for parallel variant processing within each region (default: 50)
    /// Larger chunks reduce overhead but may increase memory usage
    #[arg(long, default_value_t = 50)]
    variant_chunk_size: usize,

    /// Enable debug output with detailed variant analysis information
    #[arg(long, default_value_t = false)]
    debug: bool,

    /// Breakend span tolerance in base pairs (default: 1000)
    /// Breakpoints within this distance from SA alignment spans are considered matches
    #[arg(long, default_value_t = 1000)]
    breakend_span_tolerance: usize,
}

/// Prints detailed information about the current thread pool configuration
/// This helps confirm that parallel processing is working as expected
fn print_thread_pool_info() {
    let num_threads = rayon::current_num_threads();
    let thread_index = rayon::current_thread_index();
    println!("Thread pool configured:");
    println!("  Total worker threads: {num_threads}");
    println!("  Current thread index: {thread_index:?}");

    // Print system information for context
    let cpu_count = num_cpus::get();
    println!("  System CPU cores: {cpu_count}");
    println!(
        "  Thread pool efficiency: {:.1}% of available cores",
        (num_threads as f64 / cpu_count as f64) * 100.0
    );
}

/// Parameters for processing a BED region
struct ProcessRegionParams {
    region_idx: usize,
    bed_region: BedRegion,
    vcf_path: Arc<PathBuf>,
    bam_path: Arc<PathBuf>,
    output_writer: Arc<std::sync::Mutex<BgzOutput>>,
    variant_chunk_size: usize,
    debug: bool,
    breakend_span_tolerance: usize,
}

/// Process a single BED region with asynchronous parallel processing
///
/// This function orchestrates the complete read extraction and variant analysis workflow
/// for a single genomic region. It represents the core computational unit of the pipeline.
async fn process_bed_region_async(
    params: ProcessRegionParams,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    let ProcessRegionParams {
        region_idx,
        bed_region,
        vcf_path,
        bam_path,
        output_writer,
        variant_chunk_size,
        debug,
        breakend_span_tolerance,
    } = params;
    println!(
        "Processing BED region #{}: {} (async with indexed access)",
        region_idx, bed_region.region_string
    );

    // PHASE 1: VARIANT DISCOVERY
    // First, get all variants in this region
    let variants = {
        let vcf_path = Arc::clone(&vcf_path);
        let bed_region = bed_region.clone();
        tokio::task::spawn_blocking(move || {
            query_vcf_variants_for_region(
                &vcf_path,
                &bed_region.chrom,
                bed_region.start_pos,
                bed_region.end_pos,
            )
        })
    }
    .await??;

    println!(
        "  Found {} variants in region (indexed VCF)",
        variants.len()
    );

    // Show detailed variant information including depth data
    if debug {
        for (i, variant) in variants.iter().enumerate().take(5) {
            let depth_info = match (variant.ref_depth, variant.alt_depth, variant.total_depth) {
                (Some(ref_d), Some(alt_d), _) => format!(" [REF={ref_d}, ALT={alt_d}]"),
                (_, _, Some(total_d)) => format!(" [Total={total_d}]"),
                _ => " [No depth info]".to_string(),
            };

            println!(
                "    DEBUG: Variant {}: {}:{} {}>{} ({}){} - desc: '{}'",
                i + 1,
                variant.chrom,
                variant.pos,
                variant.ref_allele,
                variant.alt_allele,
                variant.variant_type,
                depth_info,
                variant.description
            );
        }
        if variants.len() > 5 {
            println!(
                "    DEBUG: ... and {} more variants (showing positions: {})",
                variants.len() - 5,
                variants
                    .iter()
                    .skip(5)
                    .take(10)
                    .map(|v| v.pos.to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }

    // PHASE 2: EARLY TERMINATION CHECK
    // Skip analysis if no variants present (optimization for sparse regions)
    if variants.is_empty() {
        println!(
            "  No variants found in region {}, skipping...",
            bed_region.region_string
        );
        return Ok(0);
    }

    // PHASE 3: VARIANT-SPECIFIC READ EXTRACTION
    let mut lines_written = 0;
    let mut all_output_lines = Vec::<String>::new();

    // Process variants in parallel chunks for better CPU utilization
    let processing_start = std::time::Instant::now();
    let total_variants = variants.len();

    for (chunk_idx, variant_chunk) in variants.chunks(variant_chunk_size).enumerate() {
        let chunk_start = std::time::Instant::now();
        println!(
            "  Processing variant chunk {} ({} variants)...",
            chunk_idx + 1,
            variant_chunk.len()
        );

        // PARALLEL VARIANT PROCESSING:
        let chunk_results: Vec<Vec<String>> = {
            let bam_path = Arc::clone(&bam_path);
            let bed_region = bed_region.clone();

            // Add thread monitoring and timing
            let active_threads = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
            let max_threads_seen = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));

            let results = variant_chunk
                .par_iter()  // Convert to parallel iterator for concurrent processing
                .map(|variant| {
                    // THREAD MONITORING: Track active thread usage
                    let current_thread_count = active_threads.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                    let current_thread_id = rayon::current_thread_index().unwrap_or(999);

                    // Log thread activity for monitoring (can be disabled for performance if needed)
                    if variant.pos % 1000 == 0 {  // Log every 1000th variant to avoid spam
                        println!("      [Thread {}] Processing variant at {}", current_thread_id, variant.pos);
                    }

                    // Update maximum threads seen
                    loop {
                        let current_max = max_threads_seen.load(std::sync::atomic::Ordering::Relaxed);
                        if current_thread_count <= current_max ||
                            max_threads_seen.compare_exchange_weak(
                                current_max,
                                current_thread_count,
                                std::sync::atomic::Ordering::Relaxed,
                                std::sync::atomic::Ordering::Relaxed
                            ).is_ok() {
                            break;
                        }
                    }

                    // VARIANT-SPECIFIC READ EXTRACTION:
                    // Define a focused region around the variant position
                    let variant_start = variant.pos.saturating_sub(1);
                    let variant_end = variant.pos + 1;

                    // Show detailed variant and query information
                    if debug {
                        println!("      DEBUG: Processing variant {}:{} {}>{} ({}) - '{}'",
                                bed_region.chrom, variant.pos, variant.ref_allele, variant.alt_allele,
                                variant.variant_type, variant.description);
                        println!("      DEBUG: Querying BAM region {}:{}-{} for this variant",
                                bed_region.chrom, variant_start, variant_end);
                    }

                    // Extract BAM records that might span this variant position
                    let variant_records = match query_bam_records_for_region(
                        &bam_path,
                        &bed_region.chrom,
                        variant_start,
                        variant_end,
                        debug
                    ) {
                        Ok(records) => {
                            if debug {
                                println!("      DEBUG: Found {} reads for variant {}:{} {}>{}",
                                        records.len(), bed_region.chrom, variant.pos, variant.ref_allele, variant.alt_allele);

                                // Show details of first few reads found including sequences
                                if !records.is_empty() {
                                    println!("      DEBUG: Sample reads found:");
                                    for (i, (record, read_id, start_pos, end_pos, timestamp, mapq)) in records.iter().enumerate().take(3) {
                                        // Extract the sequence from the BAM record
                                        let sequence = {
                                            let seq = record.sequence();
                                            // Convert sequence to string, showing first 50 bases
                                            let mut seq_str = String::new();
                                            let seq_len = seq.len();
                                            let display_len = seq_len.min(50);

                                            for i in 0..display_len {
                                                if let Some(base) = seq.get(i) {
                                                    seq_str.push(match base {
                                                        // Try ASCII encoding first (most common)
                                                        b'A' => 'A',
                                                        b'C' => 'C',
                                                        b'G' => 'G',
                                                        b'T' => 'T',
                                                        b'N' => 'N',
                                                        // Also try the standard BAM 4-bit encoding as fallback
                                                        0 => '=',  // Reference skip
                                                        1 => 'A',  // A = 1
                                                        2 => 'C',  // C = 2
                                                        4 => 'G',  // G = 4
                                                        8 => 'T',  // T = 8
                                                        15 => 'N', // N = 15
                                                        3 => 'G',  // Sometimes G = 3
                                                        5 => 'T',  // Sometimes T = 5
                                                        _ => {
                                                            // For any other value, use as ASCII if printable
                                                            if (32..=126).contains(&base) {
                                                                base as char
                                                            } else {
                                                                '?'
                                                            }
                                                        },
                                                    });
                                                } else {
                                                    seq_str.push('?');
                                                }
                                            }

                                            if seq_len > 50 {
                                                format!("{seq_str}... (len={seq_len})")
                                            } else {
                                                format!("{seq_str} (len={seq_len})")
                                            }
                                        };

                                        println!("        Read {}: {} at {}:{}-{} (mapq={}, ts={})",
                                                i + 1, read_id, bed_region.chrom, start_pos, end_pos, mapq,
                                                if timestamp.len() > 20 { &timestamp[..20] } else { timestamp });
                                        println!("          Sequence: {sequence}");
                                    }
                                    if records.len() > 3 {
                                        println!("        ... and {} more reads", records.len() - 3);
                                    }
                                }
                            }
                            records
                        },
                        Err(e) => {
                            eprintln!("Warning: Failed to query reads for variant at {}:{}: {:?}",
                                    bed_region.chrom, variant.pos, e);
                            active_threads.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                            return Vec::new();
                        }
                    };

                    // SEQUENCE-BASED VARIANT ANALYSIS:
                    let mut reads_processed = 0;
                    let mut reads_spanning = 0;
                    let mut reads_with_definitive_match = 0;

                    let result: Vec<String> = variant_records
                        .iter()
                        .filter_map(|(record, read_id, start_pos, end_pos, timestamp, mapping_quality)| {
                            reads_processed += 1;

                            // DEBUG: Show read details for first few reads
                            if reads_processed <= 3 {
                                // Extract a portion of the read sequence for debug
                                let seq_preview = {
                                    let seq = record.sequence();
                                    let mut seq_str = String::new();
                                    let display_len = seq.len().min(30); // Show first 30 bases

                                    for i in 0..display_len {
                                        if let Some(base) = seq.get(i) {
                                            seq_str.push(match base {
                                                // Try ASCII encoding first (most common)
                                                b'A' => 'A',
                                                b'C' => 'C',
                                                b'G' => 'G',
                                                b'T' => 'T',
                                                b'N' => 'N',
                                                // Also try the standard BAM 4-bit encoding as fallback
                                                0 => '=',  // Reference skip
                                                1 => 'A',  // A = 1
                                                2 => 'C',  // C = 2
                                                4 => 'G',  // G = 4
                                                8 => 'T',  // T = 8
                                                15 => 'N', // N = 15
                                                3 => 'G',  // Sometimes G = 3
                                                5 => 'T',  // Sometimes T = 5
                                                _ => {
                                                    // For any other value, use as ASCII if printable
                                                    if (32..=126).contains(&base) {
                                                        base as char
                                                    } else {
                                                        '?'
                                                    }
                                                },
                                            });
                                        } else {
                                            seq_str.push('?');
                                        }
                                    }

                                    if seq.len() > 30 {
                                        format!("{seq_str}...")
                                    } else {
                                        seq_str
                                    }
                                };

                                if debug {
                                    println!("        DEBUG: Read {}: {} spans {}:{}-{}, analyzing for variant {}:{} {}>{}",
                                            reads_processed, read_id, bed_region.chrom, start_pos, end_pos,
                                            bed_region.chrom, variant.pos, variant.ref_allele, variant.alt_allele);
                                    println!("        DEBUG: Read {reads_processed} sequence: {seq_preview}");
                                }
                            }

                            // STEP 1: Check if read spans the variant position (quick filter)
                            if !read_spans_variant_position(*start_pos, *end_pos, variant.pos) {
                                if debug && reads_processed <= 3 {
                                    println!("        DEBUG: Read {reads_processed} does NOT span variant position");
                                }
                                return None; // Skip reads that don't span the variant position
                            }

                            reads_spanning += 1;
                            if debug && reads_processed <= 3 {
                                println!("        DEBUG: Read {reads_processed} SPANS variant position");
                            }

                            // STEP 2: Analyze read sequence to determine variant content
                            let allele_result = analyze_read_allele_content_detailed(
                                record,
                                variant.pos,
                                &variant.ref_allele,
                                &variant.alt_allele,
                                debug,
                                breakend_span_tolerance
                            );

                            if debug && reads_processed <= 3 {
                                let result_description = match &allele_result {
                                    varclock::AlleleMatch::Reference(seq) => format!("REFERENCE match (found: {seq})"),
                                    varclock::AlleleMatch::Variant(seq) => format!("VARIANT match (found: {seq})"),
                                    varclock::AlleleMatch::Other(seq) => format!("OTHER sequence (found: {seq})"),
                                    varclock::AlleleMatch::NoSpan => "NO SPAN (read doesn't span variant)".to_string(),
                                    varclock::AlleleMatch::Deletion => "DELETION (variant in deletion)".to_string(),
                                    varclock::AlleleMatch::Indeterminate => "INDETERMINATE (couldn't analyze)".to_string(),
                                };
                                println!("        DEBUG: Read {reads_processed} allele analysis: {result_description}");
                                println!("        DEBUG: Expected - REF: '{}', ALT: '{}'",
                                        variant.ref_allele, variant.alt_allele);
                            }

                            // Include reads with definitive matches and OTVAR (skip only indeterminate/no span/deletion)
                            match &allele_result {
                                varclock::AlleleMatch::Reference(_) |
                                varclock::AlleleMatch::Variant(_) |
                                varclock::AlleleMatch::Other(_) => {
                                    // Include these in output
                                }
                                varclock::AlleleMatch::NoSpan |
                                varclock::AlleleMatch::Deletion |
                                varclock::AlleleMatch::Indeterminate => {
                                    if debug && reads_processed <= 3 {
                                        println!("        DEBUG: Read {reads_processed} excluded (no span/deletion/indeterminate)");
                                    }
                                    return None;
                                }
                            }

                            reads_with_definitive_match += 1;
                            if debug && reads_processed <= 3 {
                                let match_type = match &allele_result {
                                    varclock::AlleleMatch::Reference(_) => "REFERENCE",
                                    varclock::AlleleMatch::Variant(_) => "VARIANT",
                                    varclock::AlleleMatch::Other(_) => "OTHER/OTVAR",
                                    _ => "UNKNOWN"
                                };
                                println!("        DEBUG: Read {reads_processed} included as {match_type}");
                            }

                            let allele_match = allele_result.to_output_string();

                            // OUTPUT FORMAT: Tab-separated values with read and variant metadata
                            Some(format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                read_id,                    // Read identifier
                                timestamp,                  // Sequencing timestamp
                                allele_match,               // VARIANT:seq|REFERENCE:seq - which allele the read contains
                                variant.chrom,              // Variant chromosome
                                variant.pos,                // Variant position (1-based)
                                variant.ref_allele,         // Reference allele
                                variant.alt_allele,         // Alternative allele
                                variant.description,        // Variant description (e.g., "A>G")
                                variant.variant_type,       // Variant type (SNV/INS/DEL/etc.)
                                bed_region.region_string,   // Genomic region (chr:start-end)
                                bed_region.region_name,     // Region name/identifier
                                start_pos,                  // Read start position
                                end_pos,                    // Read end position
                                mapping_quality             // Mapping quality score
                            ))
                        })
                        .collect();  // Collect results for this variant

                    // Report filtering statistics
                    if debug {
                        println!("      DEBUG: Variant {}:{} - {} reads total, {} spanning, {} with definitive matches, {} output lines",
                                bed_region.chrom, variant.pos, reads_processed, reads_spanning, reads_with_definitive_match, result.len());
                    }

                    // SANITY CHECK: Compare our results with VCF depth information
                    let ref_count = result.iter().filter(|line| line.contains("REFERENCE:")).count();
                    let alt_count = result.iter().filter(|line| line.contains("VARIANT:")).count();
                    let other_count = result.iter().filter(|line| line.contains("OTHER:")).count();

                    if let (Some(vcf_ref), Some(vcf_alt)) = (variant.ref_depth, variant.alt_depth) {
                        println!("      SANITY CHECK: VCF reports REF={vcf_ref} ALT={vcf_alt}, we found REF={ref_count} ALT={alt_count} OTHER={other_count}");

                        let ref_diff = (ref_count as i32 - vcf_ref as i32).abs();
                        let alt_diff = (alt_count as i32 - vcf_alt as i32).abs();

                        if ref_diff > 0 || alt_diff > 0 {
                            println!("      WARNING: Read count mismatch! REF diff: {ref_diff}, ALT diff: {alt_diff}");
                        } else {
                            println!("      PASS: Read counts match VCF exactly!");
                        }
                    } else if let Some(vcf_total) = variant.total_depth {
                        println!("      SANITY CHECK: VCF reports total depth={}, we found {} supporting reads",
                                vcf_total, result.len());

                        let total_diff = (result.len() as i32 - vcf_total as i32).abs();
                        if total_diff > vcf_total as i32 / 2 { // Allow 50% difference for total depth
                            println!("      WARNING: Large difference in total depth! Diff: {total_diff}");
                        }
                    } else {
                        println!("      SANITY CHECK: No VCF depth information available for comparison");
                    }

                    // Show actual allele matches for this variant
                    if debug {
                        if !result.is_empty() {
                            println!("      DEBUG: Allele matches for variant {}:{}:", bed_region.chrom, variant.pos);
                            for (i, line) in result.iter().enumerate().take(3) {
                                // Extract just the read_id and allele_match columns for cleaner output
                                let parts: Vec<&str> = line.split('\t').collect();
                                if parts.len() >= 3 {
                                    println!("        Match {}: Read {} -> {}", i + 1, parts[0], parts[2]);
                                }
                            }
                            if result.len() > 3 {
                                println!("        ... and {} more matches", result.len() - 3);
                            }
                        } else {
                            println!("      DEBUG: No allele matches found for variant {}:{}", bed_region.chrom, variant.pos);
                        }
                    }

                    // THREAD MONITORING: Clean up thread count when done
                    active_threads.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                    result
                })
                .collect(); // Collect results for all variants in this chunk

            // THREAD USAGE REPORT: Show how much parallelism was achieved
            let max_concurrent = max_threads_seen.load(std::sync::atomic::Ordering::Relaxed);
            println!(
                "    Thread usage: max {max_concurrent} concurrent threads (thread IDs active)"
            );

            results
        };

        // Collect all lines from this chunk
        for lines in chunk_results {
            all_output_lines.extend(lines);
        }

        // CHUNK PERFORMANCE REPORT: Show timing and thread usage
        let chunk_duration = chunk_start.elapsed();
        println!(
            "    Chunk {} completed in {:.2}s ({} variants processed)",
            chunk_idx + 1,
            chunk_duration.as_secs_f64(),
            variant_chunk.len()
        );
    }

    // OVERALL PROCESSING PERFORMANCE REPORT
    let total_duration = processing_start.elapsed();
    println!(
        "  All {} variants processed in {:.2}s (avg {:.3}s per variant)",
        total_variants,
        total_duration.as_secs_f64(),
        total_duration.as_secs_f64() / total_variants as f64
    );

    // PHASE 5: ATOMIC OUTPUT WRITING
    {
        let mut writer = output_writer.lock().unwrap();
        for line in &all_output_lines {
            writer.write_line(line)?;
            lines_written += 1;
        }
    }

    println!("  Processed region #{region_idx} with {lines_written} output lines (async cached)");
    Ok(lines_written)
}

/// Main function that orchestrates the bioinformatics analysis pipeline
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments using clap derive macro
    let args = Args::parse();

    // Configure Rayon thread pool for CPU-intensive work
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    println!(
        "Starting VarClock analysis with {} threads (async mode)...",
        args.threads
    );

    // Print detailed thread pool information for monitoring
    print_thread_pool_info();
    println!("Input files:");
    println!("  BED file: {:?}", args.bed);
    println!("  VCF file: {:?}", args.vcf);
    println!("  BAM file: {:?}", args.bam);
    println!("  Output file: {:?} (BGZ-compressed)", args.output);
    println!("  Create index: {}", args.create_index);

    // Ensure output file has .gz extension for BGZ format
    let output_path = if args.output.extension().is_none_or(|ext| ext != "gz") {
        let mut new_path = args.output.clone();
        match new_path.extension() {
            Some(ext) => {
                let new_ext = format!("{}.gz", ext.to_string_lossy());
                new_path.set_extension(new_ext);
            }
            None => {
                new_path.set_extension("tsv.gz");
            }
        }
        println!("  Note: Added .gz extension to output file: {new_path:?}");
        new_path
    } else {
        args.output.clone()
    };

    // Create the BGZ output file and write TSV header row
    let mut bgz_output = BgzOutput::new(output_path.clone(), args.create_index)?;

    if let Err(e) = bgz_output.write_line("read_id\ttimestamp\tallele_match\tvariant_chrom\tvariant_pos\tvariant_ref\tvariant_alt\tvariant_description\tvariant_type\tregion\tregion_name\tread_start\tread_end\tmapping_quality") {
        println!("ERROR: Failed to write header to output file: {e:?}");
        return Err(format!("Failed to write header: {e:?}").into());
    }
    println!("Created BGZ-compressed output file and wrote header");

    // ===== BED FILE PROCESSING =====
    // Load all BED regions first
    let bed_regions = load_bed_regions(&args.bed, args.debug)?;

    // Show loaded BED regions
    if args.debug {
        println!("DEBUG: Loaded BED regions:");
        for (i, region) in bed_regions.iter().enumerate().take(10) {
            println!(
                "  Region {}: {} ({}:{}-{}, name: '{}')",
                i + 1,
                region.region_string,
                region.chrom,
                region.start_pos,
                region.end_pos,
                region.region_name
            );
        }
        if bed_regions.len() > 10 {
            println!("  ... and {} more regions", bed_regions.len() - 10);
        }
    }

    println!(
        "\n=== Processing {} BED regions with async parallelization ===",
        bed_regions.len()
    );

    // Create Arc-wrapped shared resources for async processing
    let vcf_path = Arc::new(args.vcf);
    let bam_path = Arc::new(args.bam);
    let output_mutex = Arc::new(std::sync::Mutex::new(bgz_output));

    // Global thread activity monitoring (for overall system insights)
    let global_active_tasks = Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let max_concurrent_tasks = Arc::new(std::sync::atomic::AtomicUsize::new(0));

    // Process BED regions with controlled concurrency
    let max_concurrent_regions = args.max_concurrent_regions;
    let variant_chunk_size = args.variant_chunk_size;

    let mut total_lines_written = 0;

    // Process regions in chunks to manage memory and concurrency
    let mut regions_processed = 0;

    for chunk in bed_regions.chunks(max_concurrent_regions) {
        let tasks: Vec<_> = chunk
            .iter()
            .enumerate()
            .map(|(chunk_idx, bed_region)| {
                let vcf_path = Arc::clone(&vcf_path);
                let bam_path = Arc::clone(&bam_path);
                let output_mutex = Arc::clone(&output_mutex);
                let bed_region = bed_region.clone();
                let debug = args.debug;
                let breakend_span_tolerance = args.breakend_span_tolerance;

                // Global task monitoring
                let global_tasks = Arc::clone(&global_active_tasks);
                let max_tasks = Arc::clone(&max_concurrent_tasks);

                // Calculate actual region index across all chunks
                let region_idx = regions_processed + chunk_idx + 1;

                tokio::spawn(async move {
                    // Track active task count
                    let current_tasks =
                        global_tasks.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;

                    // Update maximum concurrent tasks seen
                    loop {
                        let current_max = max_tasks.load(std::sync::atomic::Ordering::Relaxed);
                        if current_tasks <= current_max
                            || max_tasks
                                .compare_exchange_weak(
                                    current_max,
                                    current_tasks,
                                    std::sync::atomic::Ordering::Relaxed,
                                    std::sync::atomic::Ordering::Relaxed,
                                )
                                .is_ok()
                        {
                            break;
                        }
                    }

                    let params = ProcessRegionParams {
                        region_idx,
                        bed_region,
                        vcf_path,
                        bam_path,
                        output_writer: output_mutex,
                        variant_chunk_size,
                        debug,
                        breakend_span_tolerance,
                    };
                    let result = process_bed_region_async(params).await;

                    // Decrement active task count when done
                    global_tasks.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                    result
                })
            })
            .collect();

        // Wait for all tasks in this chunk to complete
        let results = future::join_all(tasks).await;

        // Process results and count lines written
        for result in results {
            match result {
                Ok(Ok(lines_written)) => {
                    total_lines_written += lines_written;
                }
                Ok(Err(e)) => {
                    eprintln!("Error processing region: {e:?}");
                }
                Err(e) => {
                    eprintln!("Task panicked: {e:?}");
                }
            }
        }

        regions_processed += chunk.len();

        // GLOBAL TASK MONITORING REPORT
        let max_tasks_seen = max_concurrent_tasks.load(std::sync::atomic::Ordering::Relaxed);
        let remaining_tasks = global_active_tasks.load(std::sync::atomic::Ordering::Relaxed);

        println!(
            "Completed chunk with {} regions, {} total regions processed, {} total lines written",
            chunk.len(),
            regions_processed,
            total_lines_written
        );
        println!(
            "  Task concurrency: max {max_tasks_seen} concurrent async tasks, {remaining_tasks} currently active"
        );
    }

    println!("Total lines written to output: {total_lines_written}");

    // Finalize the BGZ output and create index
    println!("Finalizing BGZ-compressed output...");
    let bgz_output = Arc::try_unwrap(output_mutex)
        .map_err(|_| "Failed to unwrap output mutex")?
        .into_inner()
        .map_err(|_| "Failed to acquire output mutex")?;
    bgz_output.finalize()?;
    println!("BGZ output finalized and indexed successfully");

    println!("\n=== Analysis Complete ===");
    println!("Processed {} BED regions", bed_regions.len());

    // FINAL THREAD AND PERFORMANCE SUMMARY
    let final_max_tasks = max_concurrent_tasks.load(std::sync::atomic::Ordering::Relaxed);
    let final_max_threads = rayon::current_num_threads();
    println!("=== Thread Usage Summary ===");
    println!("  Rayon thread pool: {final_max_threads} worker threads configured");
    println!("  Maximum concurrent async tasks: {final_max_tasks} (region-level parallelism)");
    println!(
        "  Hybrid parallelism: {final_max_tasks} async tasks × {final_max_threads} Rayon threads per task"
    );
    println!(
        "  Theoretical max concurrency: {} threads",
        final_max_tasks * final_max_threads
    );
    if bed_regions.is_empty() {
        println!("WARNING: No BED regions found! Check if your BED file is empty or malformed.");
        println!("Your BED file should contain lines like:");
        println!("chr1\t10000\t20000\tregion_name");
        println!("chr1\t30000\t40000\tregion_name");
        println!("...");
        return Err("No BED regions to process. Analysis cannot continue.".into());
    }

    // ===== SUCCESSFUL COMPLETION =====
    Ok(())
}
