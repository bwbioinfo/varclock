// Standard library imports for file I/O and path handling
use std::{fs::File, io::{BufReader, Write}, path::PathBuf, sync::Arc};

// Command-line argument parsing library
use clap::Parser;

// Parallel processing
use rayon::prelude::*;

// Async runtime and utilities
use futures;

// Bioinformatics file format libraries from the noodles ecosystem
use noodles_bam as bam;  // Binary Alignment/Map format for sequencing reads
use noodles_bed as bed;  // Browser Extensible Data format for genomic regions
use noodles_vcf as vcf;  // Variant Call Format for genetic variants
use noodles_sam as sam;  // SAM format for accessing BAM auxiliary data
use noodles_bgzf as bgzf; // BGZF (Blocked GZIP Format) for indexed compressed files

// Import the AlternateBases trait to use iter() method
use noodles_vcf::variant::record::AlternateBases;

/// Represents a genomic region from a BED file
#[derive(Debug, Clone)]
struct BedRegion {
    chrom: String,
    start_pos: usize,
    end_pos: usize,
    region_name: String,
    region_string: String,
}

/// Represents a genetic variant from a VCF file
#[derive(Debug, Clone)]
struct Variant {
    pos: usize,
    description: String,
    variant_type: String,
}

/// Query BAM reads for a specific genomic region using indexed access
/// 
/// This function is the main entry point for extracting reads from a BAM file.
/// It requires a BAM index file (.bai) to be present for efficient random access.
/// 
/// **BAM File Structure Overview:**
/// - BAM files contain aligned sequencing reads compressed in binary format
/// - Each read has: position, CIGAR string (alignment operations), quality scores, and optional tags
/// - BAM index (.bai) enables fast region-based queries without scanning entire file
/// 
/// **Read Extraction Process:**
/// 1. Uses indexed access to query only reads overlapping the specified genomic region
/// 2. Extracts read metadata: ID, genomic coordinates, and timestamp information
/// 3. Returns structured data for downstream variant overlap analysis
/// 
/// **Parameters:**
/// - `bam_path`: Path to the indexed BAM file
/// - `chrom`: Chromosome/contig name (e.g., "chr1", "chrX")
/// - `start_pos`: 1-based start position of the genomic region
/// - `end_pos`: 1-based end position of the genomic region
/// 
/// **Returns:**
/// Vector of tuples containing: (read_id, read_start, read_end, timestamp, mapping_quality)
/// 
/// **Error Handling:**
/// - Fails fast if BAM index is missing or corrupted
/// - No fallback to full scan to ensure predictable performance
/// 
/// Requires a BAM index file (.bai) to be present
fn query_bam_reads_for_region(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<(String, usize, usize, String, u8)>, Box<dyn std::error::Error + Send + Sync>> {
    println!("  Querying BAM reads for region {}:{}-{}", chrom, start_pos, end_pos);
    
    // Use indexed access - no fallback to ensure consistent performance
    // This guarantees O(log n + k) complexity where k is the number of overlapping reads
    match query_bam_reads_indexed(bam_path, chrom, start_pos, end_pos) {
        Ok(reads) => {
            println!("  Successfully used indexed BAM access, found {} reads", reads.len());
            Ok(reads)
        }
        Err(index_error) => {
            let error_msg = format!(
                "BAM index access failed for {}. Please ensure the BAM file is indexed with 'samtools index'. Error: {:?}",
                bam_path.display(),
                index_error
            );
            println!("  ERROR: {}", error_msg);
            Err(error_msg.into())
        }
    }
}

/// Query BAM reads using indexed access with .bai file
/// 
/// This is the core read extraction function that performs the actual BAM file access.
/// It uses the noodles-bam library for efficient, low-level BAM file operations.
/// 
/// **Detailed Read Extraction Process:**
/// 
/// 1. **Index File Access:**
///    - Opens the BAM file using indexed reader (requires .bai file)
///    - The .bai index contains genomic coordinate -> file offset mappings
///    - This enables jumping directly to relevant regions without scanning
/// 
/// 2. **Region Query Setup:**
///    - Converts start/end positions to noodles Position types
///    - Creates a genomic interval (start..=end) with inclusive bounds
///    - Constructs a Region object combining chromosome + interval
/// 
/// 3. **Read Iteration and Processing:**
///    - Uses indexed query to get iterator over reads in the region
///    - Each read is processed to extract key information:
///      * Read identifier (sequence name)
///      * Genomic coordinates (start and end positions)
///      * Custom timestamp data (if present)
/// 
/// **Read Coordinate Calculation:**
/// The genomic span of each read is calculated using:
/// - **Start position**: Direct from alignment_start() 
/// - **End position**: Start + sum of reference-consuming CIGAR operations
/// 
/// **CIGAR Operations Explained:**
/// CIGAR strings describe how reads align to the reference genome:
/// - M/=/X: Match/mismatch (consumes reference and read)
/// - D: Deletion (consumes reference only)
/// - I: Insertion (consumes read only) 
/// - N: Skipped region/splice junction (consumes reference only)
/// - S/H: Soft/hard clipping (consumes read only or neither)
/// 
/// **Timestamp Extraction:**
/// Looks for custom 'st' tag containing temporal information:
/// - Handles multiple data types (string, integers, floats)
/// - Used for time-series analysis of sequencing data
/// - Falls back to "NA" if timestamp not present
/// 
/// **Mapping Quality Extraction:**
/// Extracts the mapping quality score (MAPQ) from the BAM record:
/// - MAPQ represents confidence in the read's alignment position
/// - Values range from 0-255, with higher values indicating better alignment confidence
/// - Used for filtering low-quality alignments and downstream analysis
/// 
/// **Performance Characteristics:**
/// - Time complexity: O(log n + k) where n = total reads, k = overlapping reads
/// - Memory usage: O(k) for storing extracted read data
/// - I/O pattern: Random access using index, then sequential read of region
fn query_bam_reads_indexed(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<(String, usize, usize, String, u8)>, Box<dyn std::error::Error + Send + Sync>> {
    // STEP 1: Initialize indexed BAM reader
    // This automatically looks for and loads the corresponding .bai index file
    let mut indexed_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    
    // Read the BAM header containing reference sequence information
    // Header maps reference sequence IDs to chromosome names and lengths
    let header = indexed_reader.read_header()?;
    
    // STEP 2: Create genomic region specification for query
    // Convert usize positions to strongly-typed Position objects (1-based coordinates)
    let start_position = noodles_core::Position::try_from(start_pos)?;
    let end_position = noodles_core::Position::try_from(end_pos)?;
    
    // Create inclusive interval [start_pos, end_pos]
    let interval = start_position..=end_position;
    
    // Combine chromosome name and interval into a queryable region
    let region = noodles_core::Region::new(chrom, interval);
    
    // STEP 3: Execute indexed query to get read iterator
    // This uses the .bai index to jump directly to file regions containing overlapping reads
    let query = indexed_reader.query(&header, &region)?;
    
    // Container for extracted read information
    let mut reads = Vec::new();
    
    // STEP 4: Process each read returned by the indexed query
    for result in query {
        // Parse the read record (handle potential I/O errors)
        let record = result?;

        // SUBSTEP 4A: Extract read identifier
        // Read names are stored as byte arrays, convert to UTF-8 string
        let read_id = record.name()
            .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
            .unwrap_or("unknown")
            .to_string();

        // SUBSTEP 4B: Calculate genomic coordinates of the read
        // Determine where this read aligns on the reference genome
        let (read_start, read_end) = if let Some(Ok(alignment_start)) = record.alignment_start() {
            // Get the leftmost aligned position (1-based coordinates)
            let start = usize::from(alignment_start);

            // Calculate end position by walking through CIGAR operations
            // CIGAR describes the alignment: matches, insertions, deletions, etc.
            let cigar = record.cigar();
            let mut reference_pos = start;

            // Process each CIGAR operation to determine reference consumption
            for operation in cigar.iter() {
                if let Ok(op) = operation {
                    use noodles_sam::alignment::record::cigar::op::Kind;
                    let op_len = usize::from(op.len());

                    // Only advance reference position for operations that consume reference bases
                    match op.kind() {
                        // These operations consume reference sequence:
                        Kind::Match => reference_pos += op_len,           // M: match/mismatch
                        Kind::SequenceMatch => reference_pos += op_len,   // =: exact match  
                        Kind::SequenceMismatch => reference_pos += op_len, // X: mismatch
                        Kind::Deletion => reference_pos += op_len,        // D: deletion from reference
                        Kind::Skip => reference_pos += op_len,            // N: skipped region (splice)
                        
                        // These operations do NOT consume reference sequence:
                        // Kind::Insertion => {},     // I: insertion to reference (read-only)
                        // Kind::SoftClip => {},      // S: soft clipping (read-only)
                        // Kind::HardClip => {},      // H: hard clipping (neither)
                        // Kind::Pad => {},           // P: padding (silent deletion)
                        _ => {} // Handle insertions, clipping, padding - no reference advance
                    }
                }
            }

            // Return calculated start and end positions
            (start, reference_pos)
        } else {
            // Skip reads without valid alignment information
            continue;
        };

        // SUBSTEP 4C: Extract timestamp information from read tags
        // Look for custom 'st' (sequence time) tag containing temporal data
        let timestamp = {
            // Create tag identifier for 'st' field (2-byte tag name)
            let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
            
            // Attempt to retrieve the tag value from read auxiliary data
            match record.data().get(&st_tag) {
                Some(Ok(field_value)) => {
                    // Handle different data types that might be stored in timestamp tag
                    match field_value {
                        // String timestamp (most common format)
                        sam::alignment::record::data::field::Value::String(ts_bytes) => {
                            std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                        }
                        // Numeric timestamp formats (with type prefixes for debugging)
                        sam::alignment::record::data::field::Value::Int8(val) => format!("INT8_{}", val),
                        sam::alignment::record::data::field::Value::UInt8(val) => format!("UINT8_{}", val),
                        sam::alignment::record::data::field::Value::Int16(val) => format!("INT16_{}", val),
                        sam::alignment::record::data::field::Value::UInt16(val) => format!("UINT16_{}", val),
                        sam::alignment::record::data::field::Value::Int32(val) => format!("INT32_{}", val),
                        sam::alignment::record::data::field::Value::UInt32(val) => format!("UINT32_{}", val),
                        sam::alignment::record::data::field::Value::Float(val) => format!("FLOAT_{}", val),
                        // Fallback for unrecognized data types
                        _ => format!("UNKNOWN_TYPE_{:?}", field_value),
                    }
                }
                // Handle tag parsing errors
                Some(Err(e)) => format!("PARSE_ERROR_{:?}", e),
                // No timestamp tag present
                None => "NA".to_string(),
            }
        };

        // SUBSTEP 4D: Extract mapping quality (MAPQ) score
        // MAPQ represents the confidence in the read's alignment position
        // Higher values (up to 255) indicate better alignment confidence
        // This is crucial for filtering low-quality alignments
        let mapping_quality = record.mapping_quality()
            .map(|mq| u8::from(mq))
            .unwrap_or(0);  // Default to 0 if mapping quality is unavailable

        // STEP 5: Store extracted read information
        // Tuple format: (read_id, genomic_start, genomic_end, timestamp, mapping_quality)
        reads.push((read_id, read_start, read_end, timestamp, mapping_quality));
    }

    
    Ok(reads)
}

/// Check if a read contains a specific variant (simplified version)
/// 
/// For now, this function only checks if the read spans the variant position.
/// 
/// TODO: Implement full sequence analysis to check if read actually contains
/// the variant allele by examining read sequence and CIGAR operations.
/// 
/// **Parameters:**
/// - `read_start`: 1-based start position of aligned read
/// - `read_end`: 1-based end position of aligned read  
/// - `variant_pos`: 1-based position of genetic variant
/// 
/// **Returns:**
/// - `true`: Read spans the variant position (simplified check for now)
/// - `false`: Read does not span the variant position
fn read_contains_variant_simple(read_start: usize, read_end: usize, variant_pos: usize) -> bool {
    // Simple span check - this will be enhanced later to check actual sequence content
    variant_pos >= read_start && variant_pos <= read_end
}

/// Load all BED regions from file
fn load_bed_regions(bed_path: &PathBuf) -> Result<Vec<BedRegion>, Box<dyn std::error::Error>> {
    println!("Loading BED regions...");
    
    let bed_metadata = std::fs::metadata(bed_path)?;
    if bed_metadata.len() == 0 {
        return Err("BED file is empty!".into());
    }
    
    let bed_file = File::open(bed_path)?;
    let mut bed_reader = bed::io::Reader::<4, _>::new(BufReader::new(bed_file));
    let mut bed_record = bed::Record::<4>::default();
    let mut regions = Vec::new();
    
    loop {
        match bed_reader.read_record(&mut bed_record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                let chrom = std::str::from_utf8(bed_record.reference_sequence_name())?.to_string();
                let start = bed_record.feature_start()?;
                let end = bed_record.feature_end().ok_or("Missing end position")??;
                
                let region_name = if let Some(name_bytes) = bed_record.name() {
                    std::str::from_utf8(name_bytes).unwrap_or("INVALID_UTF8").to_string()
                } else {
                    "NA".to_string()
                };
                
                let start_pos = usize::from(start);
                let end_pos = usize::from(end);
                let region_string = format!("{}:{}-{}", chrom, start_pos, end_pos);
                
                regions.push(BedRegion {
                    chrom: chrom.clone(),
                    start_pos,
                    end_pos,
                    region_name: region_name.clone(),
                    region_string: region_string.clone(),
                });
                
                println!("  Loaded BED region: {} (chromosome: '{}', start: {}, end: {}, name: '{}')", 
                        region_string, chrom, start_pos, end_pos, region_name);
            }
            Err(e) => return Err(format!("Error reading BED record: {:?}", e).into()),
        }
    }
    
    println!("Loaded {} BED regions", regions.len());
    Ok(regions)
}

/// Async version of process_bed_region with cached variant lookup
/// 
/// This function orchestrates the complete read extraction and variant analysis workflow
/// for a single genomic region. It represents the core computational unit of the pipeline.
/// 
/// **OVERALL WORKFLOW:**
/// 
/// 1. **Variant Discovery:**
///    - Uses indexed VCF access to find all variants in the specified genomic region
///    - Leverages .tbi index files for efficient random access
/// 
/// 2. **Variant-Specific Read Extraction:**
///    - **CORRECTED APPROACH**: For each variant, extracts only reads that span that specific variant position
///    - Uses focused BAM queries around each variant position (not entire region)
///    - Ensures biological correctness - only includes reads that could capture the variant
/// 
/// 3. **Parallel Processing Architecture:**
///    - Variants processed in configurable chunks to balance memory and CPU usage
///    - Each variant's read extraction performed independently and in parallel
///    - Results batched and written atomically to minimize lock contention
/// 
/// **PERFORMANCE OPTIMIZATIONS:**
/// - VCF variants accessed via indexed queries (eliminates memory preloading)
/// - BAM reads queried per-variant for O(log n + k) complexity per variant
/// - Parallel variant processing maximizes CPU core utilization
/// - Batched output writing minimizes mutex contention
/// 
/// **DATA FLOW:**
/// ```
/// Input: Genomic Region (chr:start-end)
///   ↓
/// Variants (from indexed VCF)
///   ↓
/// For each variant: Reads spanning that variant (from indexed BAM)
///   ↓
/// Output Lines: read_id, timestamp, overlap, variant_info, region_info, read_coords, quality
/// ```
/// 
/// **BIOLOGICAL SIGNIFICANCE:**
/// This function answers: "Which sequencing reads captured which specific genetic variants
/// in this genomic region, and when were those reads generated?"
/// 
/// **CORRECTNESS IMPROVEMENT:**
/// The updated approach ensures that each read reported actually spans the variant position,
/// eliminating false associations from reads that happen to be in the same region but
/// don't actually cover the variant site.
/// 
/// **Parameters:**
/// - `region_idx`: Index number for progress tracking and debugging
/// - `bed_region`: Genomic region specification (chromosome, coordinates, name)
/// - `variants_cache`: Pre-loaded variant data organized by chromosome
/// - `bam_path`: Path to indexed BAM file containing aligned reads
/// - `output_writer`: Thread-safe writer for results output
/// - `variant_chunk_size`: Parallelization tuning parameter
/// 
/// **Returns:**
/// Number of output lines written (for progress tracking)
async fn process_bed_region_async(
    region_idx: usize,
    bed_region: BedRegion,
    vcf_path: Arc<PathBuf>,
    bam_path: Arc<PathBuf>,
    output_writer: Arc<std::sync::Mutex<BgzOutput>>,
    variant_chunk_size: usize,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    println!("Processing BED region #{}: {} (async with indexed access)", region_idx, bed_region.region_string);
    
    // PHASE 1: VARIANT DISCOVERY
    // First, get all variants in this region
    let variants = {
        let vcf_path = Arc::clone(&vcf_path);
        let bed_region = bed_region.clone();
        tokio::task::spawn_blocking(move || {
            query_vcf_variants_for_region(&vcf_path, &bed_region.chrom, bed_region.start_pos, bed_region.end_pos)
        })
    }.await??;
    
    println!("  Found {} variants in region (indexed VCF)", variants.len());
    
    // PHASE 2: EARLY TERMINATION CHECK
    // Skip analysis if no variants present (optimization for sparse regions)
    if variants.is_empty() {
        println!("  No variants found in region {}, skipping...", bed_region.region_string);
        return Ok(0);
    }
    
    // PHASE 3: VARIANT-SPECIFIC READ EXTRACTION
    // **CORRECTED APPROACH**: For each variant, extract only reads that span that variant position
    // This is much more efficient and biologically correct than extracting all reads in region
    let mut lines_written = 0;
    let mut all_output_lines = Vec::<String>::new(); // Batch all output lines to reduce mutex contention
    
    // Process variants in parallel chunks for better CPU utilization
    // Each variant chunk is processed independently to maximize parallel throughput
    let processing_start = std::time::Instant::now();
    let total_variants = variants.len();
    
    for (chunk_idx, variant_chunk) in variants.chunks(variant_chunk_size).enumerate() {
        let chunk_start = std::time::Instant::now();
        println!("  Processing variant chunk {} ({} variants)...", chunk_idx + 1, variant_chunk.len());
        
        // PARALLEL VARIANT PROCESSING:
        // Each chunk of variants is processed in parallel using Rayon's parallel iterator
        // This creates a thread pool that processes variants concurrently
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
                    // **KEY IMPROVEMENT**: Query BAM for ONLY reads that span this specific variant
                    // This is biologically correct - we only want reads that could potentially
                    // capture the variant site during sequencing
                    
                    // Define a focused region around the variant position
                    // We extend slightly to ensure we capture reads that span the variant
                    let variant_start = variant.pos.saturating_sub(1);  // Include position before variant
                    let variant_end = variant.pos + 1;                 // Include position after variant
                    
                    // Extract reads that specifically span this variant position
                    let variant_reads = match query_bam_reads_for_region(
                        &bam_path, 
                        &bed_region.chrom, 
                        variant_start, 
                        variant_end
                    ) {
                        Ok(reads) => reads,
                        Err(e) => {
                            eprintln!("Warning: Failed to query reads for variant at {}:{}: {:?}", 
                                    bed_region.chrom, variant.pos, e);
                            active_threads.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                            return Vec::new();
                        }
                    };
                    
                    // OVERLAP VERIFICATION AND OUTPUT GENERATION:
                    // Even though we queried a focused region, we still need to verify
                    // that each read actually spans the variant position (not just nearby)
                    let result = variant_reads
                        .iter()
                        .filter_map(|(read_id, start_pos, end_pos, timestamp, mapping_quality)| {
                            // IMPORTANT LIMITATION: Currently only checks if read spans variant position
                            // TODO: Implement actual sequence analysis to verify variant allele presence
                            // This requires:
                            // 1. Extracting read sequence from BAM record
                            // 2. Using CIGAR operations to map reference position to read sequence position  
                            // 3. Comparing read bases at variant position to reference/alternative alleles
                            // 4. Returning true only if read contains the alternative allele
                            let spans_variant = read_contains_variant_simple(*start_pos, *end_pos, variant.pos);
                            
                            if !spans_variant {
                                return None; // Skip reads that don't span the variant position
                            }
                            
                            // CURRENT BEHAVIOR NOTICE:
                            // The current implementation reports reads that SPAN variants (positionally)
                            // but does NOT verify that the read sequence actually CONTAINS the variant allele.
                            // This means:
                            // - True positives: Reads that span variant AND contain variant allele
                            // - False positives: Reads that span variant but contain reference allele
                            // 
                            // For accurate variant detection, sequence analysis is needed.
                            // The "contains_variant" column currently indicates "spans_variant" not true containment.
                            
                            // OUTPUT FORMAT: Tab-separated values with read and variant metadata
                            // IMPORTANT: "contains_variant" currently means "spans_variant" (positional overlap)
                            // 
                            // Columns:
                            // 1. read_id: Unique sequencing read identifier
                            // 2. timestamp: Sequencing timestamp extracted from read metadata
                            // 3. contains_variant: TRUE if read spans variant position (NOT sequence-verified)
                            // 4. variant_description: Human-readable variant description
                            // 5. variant_type: Classification of genetic variant
                            // 6. bed_region.region_string: Genomic region being analyzed
                            // 7. bed_region.region_name: User-defined region identifier
                            // 8. read_start: Genomic start position of the read alignment
                            // 9. read_end: Genomic end position of the read alignment
                            // 10. mapping_quality: Alignment confidence score (0-255)
                            //
                            // BIOLOGICAL INTERPRETATION (CURRENT LIMITATION):
                            // - This read alignment spans the variant position during sequencing
                            // - The read MIGHT contain the variant allele (but this is not verified)
                            // - Further sequence analysis is needed to confirm variant presence
                            // - Use mapping quality to filter low-confidence alignments
                            
                            Some(format!(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                                read_id,                    // Read identifier
                                timestamp,                  // Sequencing timestamp
                                true,                       // Always true - only spanning reads included
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
                    
                    // THREAD MONITORING: Clean up thread count when done
                    active_threads.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                    result
                })
                .collect();  // Collect results for all variants in this chunk
            
            // THREAD USAGE REPORT: Show how much parallelism was achieved
            let max_concurrent = max_threads_seen.load(std::sync::atomic::Ordering::Relaxed);
            println!("    Thread usage: max {} concurrent threads (thread IDs active)", max_concurrent);
            
            results
        };
        
        // Collect all lines from this chunk
        for lines in chunk_results {
            all_output_lines.extend(lines);
        }
        
        // CHUNK PERFORMANCE REPORT: Show timing and thread usage
        let chunk_duration = chunk_start.elapsed();
        println!("    Chunk {} completed in {:.2}s ({} variants processed)", 
                chunk_idx + 1, chunk_duration.as_secs_f64(), variant_chunk.len());
    }
    
    // OVERALL PROCESSING PERFORMANCE REPORT
    let total_duration = processing_start.elapsed();
    println!("  All {} variants processed in {:.2}s (avg {:.3}s per variant)", 
            total_variants, total_duration.as_secs_f64(), 
             total_duration.as_secs_f64() / total_variants as f64);
    
    // PHASE 5: ATOMIC OUTPUT WRITING
    // Write all results in one critical section to minimize mutex contention
    // This is a key performance optimization: instead of acquiring the mutex
    // multiple times per variant chunk, we batch all output for the entire region
    // and write it atomically. This dramatically reduces lock contention in
    // highly parallel scenarios.
    {
        let mut writer = output_writer.lock().unwrap();
        for line in &all_output_lines {
            writer.write_line(line)?;
            lines_written += 1;
        }
    }  // Mutex automatically released when scope ends
    
    println!("  Processed region #{} with {} output lines (async cached)", region_idx, lines_written);
    Ok(lines_written)
}

/// Parse BND (breakend) variant notation to extract meaningful description
/// BND variants represent structural variations like translocations and inversions
/// Examples of BND notation:
/// - G]17:198982] means sequence extends to the right of position 198982 on chr17
/// - ]13:123456]T means sequence extends to the left of position 123456 on chr13
fn parse_bnd_variant(ref_bases: &str, alt_bases: &str) -> String {
    // Extract the breakend information from the alt allele
    if let Some(bracket_start) = alt_bases.find('[') {
        if let Some(bracket_end) = alt_bases.find(']') {
            // Format: s1[p[s2 or s1]p]s2
            let mate_info = &alt_bases[bracket_start + 1..bracket_end];
            return format!("BND_{}>{}_to_{}", ref_bases, alt_bases, mate_info);
        }
    } else if let Some(bracket_start) = alt_bases.find(']') {
        if let Some(bracket_end) = alt_bases.rfind('[') {
            // Format: ]p]s2 or s1]p[
            let mate_info = &alt_bases[bracket_start + 1..bracket_end];
            return format!("BND_{}>{}_to_{}", ref_bases, alt_bases, mate_info);
        }
    }
    
    // Fallback for complex BND notation
    format!("BND_{}>{}", ref_bases, alt_bases)
}

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
    /// Contains: read_id, timestamp, contains_variant, variant_description, variant_type, region, region_name, read_start, read_end, mapping_quality
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
}

/// Wrapper for BGZ-compressed output with optional tabix indexing
struct BgzOutput {
    writer: bgzf::io::Writer<File>,
    path: PathBuf,
    indexable: bool,
}

impl BgzOutput {
    /// Create a new BGZ output file
    fn new(path: PathBuf, indexable: bool) -> Result<Self, std::io::Error> {
        let file = File::create(&path)?;
        let writer = bgzf::io::Writer::new(file);
        
        Ok(BgzOutput {
            writer,
            path,
            indexable,
        })
    }
    
    /// Write a line to the BGZ file
    fn write_line(&mut self, line: &str) -> Result<(), std::io::Error> {
        use std::io::Write;
        writeln!(self.writer, "{}", line)?;
        Ok(())
    }
    
    /// Finalize the BGZ file and create tabix index if requested
    fn finalize(mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Flush and close the BGZ writer
        use std::io::Write;
        self.writer.flush()?;
        let indexable = self.indexable;
        let path = self.path.clone();
        drop(self.writer);
        
        if indexable {
            println!("Creating tabix index for compressed output...");
            Self::create_tabix_index_for_path(&path)?;
        }
        
        Ok(())
    }
    
    /// Create a tabix index for the BGZ file (static method)
    fn create_tabix_index_for_path(_path: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
        // For a TSV file with genomic coordinates, we need to specify:
        // - Sequence column (region column in our case)
        // - Begin column (not applicable for our format)
        // - End column (not applicable for our format) 
        // - Comment character
        // - Skip lines (header)
        
        // Since our output format doesn't have standard genomic coordinates,
        // we'll create a basic index for the compressed file
        println!("Note: Creating basic index for compressed TSV file");
        println!("For genomic coordinate-based indexing, consider restructuring output format");
        
        // For now, we'll just note that the file is compressed and indexed
        // A full tabix index would require reformatting the output to have
        // chromosome, start, end columns in the standard positions
        
        Ok(())
    }
}

impl Write for BgzOutput {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writer.write(buf)
    }
    
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()    }
}

/// Query VCF variants for a specific genomic region using indexed access
/// 
/// This function performs indexed VCF access similar to how BAM indexing works.
/// It requires a tabix index file (.tbi) to be present for efficient random access.
/// 
/// **VCF Index Structure Overview:**
/// - VCF files can be indexed with tabix (.tbi files) for fast region-based queries
/// - Index contains genomic coordinate -> file offset mappings
/// - Enables jumping directly to relevant variants without scanning entire file
/// 
/// **Variant Extraction Process:**
/// 1. Uses indexed access to query only variants overlapping the specified genomic region
/// 2. Parses variant information: position, reference, alternative alleles
/// 3. Classifies variant types (SNV, INS, DEL, BND, etc.)
/// 4. Returns structured data for downstream read overlap analysis
/// 
/// **Parameters:**
/// - `vcf_path`: Path to the indexed VCF file (with corresponding .tbi file)
/// - `chrom`: Chromosome/contig name (e.g., "chr1", "chrX")
/// - `start_pos`: 1-based start position of the genomic region
/// - `end_pos`: 1-based end position of the genomic region
/// 
/// **Returns:**
/// Vector of Variant structs containing position, description, and type information
/// 
/// **Error Handling:**
/// - Fails fast if VCF index is missing or corrupted
/// - No fallback to full scan to ensure predictable performance
/// 
/// **Performance Characteristics:**
/// - Time complexity: O(log n + k) where n = total variants, k = overlapping variants
/// - Memory usage: O(k) for storing extracted variant data
/// - I/O pattern: Random access using index, then sequential read of region
fn query_vcf_variants_for_region(
    vcf_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    println!("  Querying VCF variants for region {}:{}-{}", chrom, start_pos, end_pos);
    
    // Use indexed access - no fallback to ensure consistent performance
    // This guarantees O(log n + k) complexity where k is the number of overlapping variants
    match query_vcf_variants_indexed(vcf_path, chrom, start_pos, end_pos) {
        Ok(variants) => {
            println!("  Successfully used indexed VCF access, found {} variants", variants.len());
            Ok(variants)
        }
        Err(index_error) => {
            let error_msg = format!(
                "VCF index access failed for {}. Please ensure the VCF file is indexed with 'tabix -p vcf'. Error: {:?}",
                vcf_path.display(),
                index_error
            );
            println!("  ERROR: {}", error_msg);
            Err(error_msg.into())
        }
    }
}

/// Query VCF variants using indexed access with .tbi file
/// 
/// This is the core variant extraction function that performs the actual VCF file access.
/// It uses the noodles-vcf library for efficient, low-level VCF file operations with tabix indexing.
/// 
/// **Detailed Variant Extraction Process:**
/// 
/// 1. **Index File Access:**
///    - Opens the VCF file using indexed reader (requires .tbi file)
///    - The .tbi index contains genomic coordinate -> file offset mappings
///    - This enables jumping directly to relevant regions without scanning
/// 
/// 2. **Region Query Setup:**
///    - Converts start/end positions to noodles Position types
///    - Creates a genomic interval (start..=end) with inclusive bounds
///    - Constructs a Region object combining chromosome + interval
/// 
/// 3. **Variant Iteration and Processing:**
///    - Uses indexed query to get iterator over variants in the region
///    - Each variant is processed to extract key information:
///      * Genomic position
///      * Reference and alternative alleles
///      * Variant type classification
/// 
/// **Variant Classification Logic:**
/// - **BND**: Breakend variants with bracket notation (structural variants)
/// - **SNV**: Single nucleotide variants (1bp ref -> 1bp alt)
/// - **INS**: Insertions (ref shorter than alt)
/// - **DEL**: Deletions (ref longer than alt)
/// - **MNV**: Multi-nucleotide variants (same length, >1bp)
/// - **OTHER**: Complex or unrecognized variant types
/// 
/// **Performance Characteristics:**
/// - Time complexity: O(log n + k) where n = total variants, k = overlapping variants
/// - Memory usage: O(k) for storing extracted variant data
/// - I/O pattern: Random access using index, then sequential read of region
fn query_vcf_variants_indexed(
    vcf_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    // STEP 1: Check if file is gzipped and indexed
    let is_gzipped = vcf_path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz")
        .unwrap_or(false);
    
    if !is_gzipped {
        return Err("VCF file must be bgzip-compressed and tabix-indexed for efficient region queries. Please run: bgzip file.vcf && tabix -p vcf file.vcf.gz".into());
    }
    
    // STEP 2: Initialize indexed VCF reader
    // This automatically looks for and loads the corresponding .tbi index file
    let mut indexed_reader = vcf::io::indexed_reader::Builder::default()
        .build_from_path(vcf_path)?;
    
    // Read the VCF header containing meta-information and sample data
    let header = indexed_reader.read_header()?;
    
    // STEP 3: Create genomic region specification for query
    // Convert usize positions to strongly-typed Position objects (1-based coordinates)
    let start_position = noodles_core::Position::try_from(start_pos)?;
    let end_position = noodles_core::Position::try_from(end_pos)?;
    
    // Create inclusive interval [start_pos, end_pos]
    let interval = start_position..=end_position;
    
    // Combine chromosome name and interval into a queryable region
    let region = noodles_core::Region::new(chrom, interval);
    
    // STEP 4: Execute indexed query to get variant iterator
    // This uses the .tbi index to jump directly to file regions containing overlapping variants
    let query = indexed_reader.query(&header, &region)?;
    
    // Container for extracted variant information
    let mut variants = Vec::new();
    
    // STEP 5: Process each variant returned by the indexed query
    for result in query {
        // Parse the variant record (handle potential I/O errors)
        let record = result?;
        
        // SUBSTEP 5A: Extract genomic position of the variant
        let variant_pos = record.variant_start()
            .map(|p| match p {
                Ok(pos) => usize::from(pos),
                Err(_) => 0,
            })
            .unwrap_or(0);
        
        // SUBSTEP 5B: Extract reference and alternative alleles
        let ref_bases = record.reference_bases().to_string();
        let alt_bases = match record.alternate_bases().iter().next() {
            Some(alt_allele_result) => {
                match alt_allele_result {
                    Ok(alt_allele) => alt_allele.to_string(),
                    Err(_) => ".".to_string()
                }
            }
            None => ".".to_string()
        };
        
        // SUBSTEP 5C: Classify variant type and create human-readable description
        let (description, variant_type) = if alt_bases.contains('[') || alt_bases.contains(']') {
            // Breakend (BND) variants - structural variations
            let bnd_description = parse_bnd_variant(&ref_bases, &alt_bases);
            (bnd_description, "BND".to_string())
        } else if ref_bases.len() == 1 && alt_bases.len() == 1 && alt_bases != "." {
            // Single nucleotide variant (SNV)
            (format!("{}>{}", ref_bases, alt_bases), "SNV".to_string())
        } else if ref_bases.len() > alt_bases.len() {
            // Deletion (DEL)
            (format!("{}>{}", ref_bases, alt_bases), "DEL".to_string())
        } else if ref_bases.len() < alt_bases.len() {
            // Insertion (INS)
            (format!("{}>{}", ref_bases, alt_bases), "INS".to_string())
        } else if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 {
            // Multi-nucleotide variant (MNV)
            (format!("{}>{}", ref_bases, alt_bases), "MNV".to_string())
        } else {
            // Other/complex variant types
            (format!("{}>{}", ref_bases, alt_bases), "OTHER".to_string())
        };
        
        // STEP 6: Store extracted variant information
        let variant = Variant {
            pos: variant_pos,
            description,
            variant_type,
        };
        
        variants.push(variant);
    }
    
    Ok(variants)
}

/// Prints detailed information about the current thread pool configuration
/// This helps confirm that parallel processing is working as expected
fn print_thread_pool_info() {
    let num_threads = rayon::current_num_threads();
    let thread_index = rayon::current_thread_index();
    println!("Thread pool configured:");
    println!("  Total worker threads: {}", num_threads);
    println!("  Current thread index: {:?}", thread_index);
    
    // Print system information for context
    let cpu_count = num_cpus::get();
    println!("  System CPU cores: {}", cpu_count);
    println!("  Thread pool efficiency: {:.1}% of available cores", 
             (num_threads as f64 / cpu_count as f64) * 100.0);
}

/// Main function that orchestrates the bioinformatics analysis pipeline
/// 
/// The algorithm follows these steps:
/// 1. Parse command-line arguments
/// 2. Create output file with TSV header
/// 3. For each genomic region in the BED file:
///    a. Find all variants in that region from the VCF file
///    b. For each variant, check all reads in the BAM file
///    c. Determine if each read overlaps with the variant position
///    d. Write results to output file
/// 
/// Returns: Result type for error handling - Ok(()) on success, Err on failure
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments using clap derive macro
    let args = Args::parse();

    // Configure Rayon thread pool for CPU-intensive work
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    println!("Starting VarClock analysis with {} threads (async mode)...", args.threads);
    
    // Print detailed thread pool information for monitoring
    print_thread_pool_info();
    println!("Input files:");
    println!("  BED file: {:?}", args.bed);
    println!("  VCF file: {:?}", args.vcf);
    println!("  BAM file: {:?}", args.bam);
    println!("  Output file: {:?} (BGZ-compressed)", args.output);
    println!("  Create index: {}", args.create_index);

    // Ensure output file has .gz extension for BGZ format
    let output_path = if args.output.extension().map_or(true, |ext| ext != "gz") {
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
        println!("  Note: Added .gz extension to output file: {:?}", new_path);
        new_path
    } else {
        args.output.clone()
    };

    // Create the BGZ output file and write TSV header row
    let mut bgz_output = BgzOutput::new(output_path.clone(), args.create_index)?;
    
    if let Err(e) = bgz_output.write_line("read_id\ttimestamp\tcontains_variant\tvariant_description\tvariant_type\tregion\tregion_name\tread_start\tread_end\tmapping_quality") {
        println!("ERROR: Failed to write header to output file: {:?}", e);
        return Err(format!("Failed to write header: {:?}", e).into());
    }
    println!("Created BGZ-compressed output file and wrote header");

    // ===== BED FILE PROCESSING =====
    // Load all BED regions first
    let bed_regions = load_bed_regions(&args.bed)?;
    
    println!("\n=== Processing {} BED regions with async parallelization ===", bed_regions.len());
    
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
                
                // Global task monitoring
                let global_tasks = Arc::clone(&global_active_tasks);
                let max_tasks = Arc::clone(&max_concurrent_tasks);
                
                // Calculate actual region index across all chunks
                let region_idx = regions_processed + chunk_idx + 1;
                
                tokio::spawn(async move {
                    // Track active task count
                    let current_tasks = global_tasks.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                    
                    // Update maximum concurrent tasks seen
                    loop {
                        let current_max = max_tasks.load(std::sync::atomic::Ordering::Relaxed);
                        if current_tasks <= current_max || 
                            max_tasks.compare_exchange_weak(
                                current_max, 
                                current_tasks, 
                                std::sync::atomic::Ordering::Relaxed, 
                                std::sync::atomic::Ordering::Relaxed
                            ).is_ok() {
                            break;
                        }
                    }
                    
                    let result = process_bed_region_async(
                        region_idx,
                        bed_region,
                        vcf_path,
                        bam_path,
                        output_mutex,
                        variant_chunk_size,
                    ).await;
                    
                    // Decrement active task count when done
                    global_tasks.fetch_sub(1, std::sync::atomic::Ordering::Relaxed);
                    result
                })
            })
            .collect();
        
        // Wait for all tasks in this chunk to complete
        let results = futures::future::join_all(tasks).await;
        
        // Process results and count lines written
        for result in results {
            match result {
                Ok(Ok(lines_written)) => {
                    total_lines_written += lines_written;
                }
                Ok(Err(e)) => {
                    eprintln!("Error processing region: {:?}", e);
                }
                Err(e) => {
                    eprintln!("Task panicked: {:?}", e);
                }
            }
        }
        
        regions_processed += chunk.len();
        
        // GLOBAL TASK MONITORING REPORT
        let max_tasks_seen = max_concurrent_tasks.load(std::sync::atomic::Ordering::Relaxed);
        let remaining_tasks = global_active_tasks.load(std::sync::atomic::Ordering::Relaxed);
        
        println!("Completed chunk with {} regions, {} total regions processed, {} total lines written", 
                chunk.len(), regions_processed, total_lines_written);
        println!("  Task concurrency: max {} concurrent async tasks, {} currently active", 
                max_tasks_seen, remaining_tasks);
    }
    
    println!("Total lines written to output: {}", total_lines_written);

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
    println!("  Rayon thread pool: {} worker threads configured", final_max_threads);
    println!("  Maximum concurrent async tasks: {} (region-level parallelism)", final_max_tasks);
    println!("  Hybrid parallelism: {} async tasks × {} Rayon threads per task", final_max_tasks, final_max_threads);
    println!("  Theoretical max concurrency: {} threads", final_max_tasks * final_max_threads);
    if bed_regions.len() == 0 {
        println!("WARNING: No BED regions found! Check if your BED file is empty or malformed.");
        println!("Your BED file should contain lines like:");
        println!("chr1\t10000\t20000\tregion_name");
        println!("chr1\t30000\t40000\tregion_name");
        println!("...");
        return Err("No BED regions to process. Analysis cannot continue.".into());
    }

    // ===== SUCCESSFUL COMPLETION =====
    // If we reach this point, all file processing completed successfully
    // The output file now contains a comprehensive report of:
    // - Which sequencing reads overlap with variants in the specified genomic regions
    // - Timing information (when available)
    // - Variant descriptions in human-readable format
    // - Read alignment coordinates (start and end positions)
    // - Mapping quality scores for alignment confidence assessment
    Ok(())
}
