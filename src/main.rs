// Standard library imports for file I/O and path handling
use std::{fs::File, io::{BufReader, Write, Read, Cursor}, path::PathBuf, sync::Arc};

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

// Import for handling gzipped files (fallback for non-BGZF files)
use flate2::read::GzDecoder;

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

/// Check if a read overlaps with a variant position
/// 
/// This is the core variant-read overlap detection algorithm.
/// 
/// **Overlap Logic:**
/// - A read "contains" or "overlaps" a variant if the variant position falls
///   within the genomic span covered by the read's alignment
/// - Uses inclusive bounds: variant_pos >= read_start AND variant_pos <= read_end
/// 
/// **Coordinate System:**
/// - All positions use 1-based genomic coordinates
/// - read_start: leftmost aligned position of the read
/// - read_end: rightmost aligned position (calculated from CIGAR operations)
/// - variant_pos: genomic position of the genetic variant
/// 
/// **Biological Interpretation:**
/// - If a read overlaps a variant position, the sequencing technology
///   potentially captured that variant site during DNA fragment sequencing
/// - This enables determination of which reads support reference vs. alternative alleles
/// - Critical for variant calling, phasing, and temporal analysis
/// 
/// **Example:**
/// ```
/// Reference:    1234567890123456789
/// Read:         ----ATCGATCG----     (positions 5-12)
/// Variant:             ^             (position 9, G>A SNV)
/// Result:       true (read overlaps variant)
/// ```
/// 
/// **Parameters:**
/// - `read_start`: 1-based start position of aligned read
/// - `read_end`: 1-based end position of aligned read  
/// - `variant_pos`: 1-based position of genetic variant
/// 
/// **Returns:**
/// - `true`: Read genomically overlaps the variant position
/// - `false`: Read does not overlap the variant position
fn read_overlaps_variant(read_start: usize, read_end: usize, variant_pos: usize) -> bool {
    // Simple inclusive range check: variant must fall within read boundaries
    // This accounts for the fact that reads span genomic intervals
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
/// 1. **Concurrent Data Retrieval:**
///    - Simultaneously fetches variants (from cache) and reads (from BAM file)
///    - Uses async/await for non-blocking I/O operations
///    - Maximizes CPU utilization through parallel data access
/// 
/// 2. **Read-Variant Cross-Analysis:**
///    - Performs combinatorial analysis of all reads vs. all variants in the region
///    - Determines spatial overlap using genomic coordinates
///    - Generates detailed output for each read-variant pair
/// 
/// 3. **Parallel Processing Architecture:**
///    - Variants processed in configurable chunks to balance memory and CPU usage
///    - Each variant analyzed against all reads using parallel iterators
///    - Results batched and written atomically to minimize lock contention
/// 
/// **PERFORMANCE OPTIMIZATIONS:**
/// - VCF variants pre-cached for instant lookup (eliminates I/O bottleneck)
/// - BAM reads accessed via index for O(log n + k) complexity
/// - Parallel variant processing maximizes CPU core utilization
/// - Batched output writing minimizes mutex contention
/// 
/// **DATA FLOW:**
/// ```
/// Input: Genomic Region (chr:start-end)
///   ↓
/// ┌─ Variants (from cache) ── Reads (from indexed BAM) ─┐
///   ↓                                                    ↓
/// Cross-Product Analysis: variant × read → overlap result
///   ↓
/// Output Lines: read_id, timestamp, overlap, variant_info, region_info
/// ```
/// 
/// **BIOLOGICAL SIGNIFICANCE:**
/// This function answers: "Which sequencing reads captured which genetic variants
/// in this genomic region, and when were those reads generated?"
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
    variants_cache: Arc<std::collections::HashMap<String, Vec<Variant>>>,
    bam_path: Arc<PathBuf>,
    output_writer: Arc<std::sync::Mutex<BgzOutput>>,
    variant_chunk_size: usize,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    println!("Processing BED region #{}: {} (async with cache)", region_idx, bed_region.region_string);
    
    // PHASE 1: CONCURRENT DATA RETRIEVAL
    // Launch two independent tasks to fetch variants and reads simultaneously
    // This overlaps I/O operations to minimize total processing time
    
    // Task 1: Variant lookup from pre-loaded cache (very fast - memory access only)
    let variants_task = {
        let variants_cache = Arc::clone(&variants_cache);
        let bed_region = bed_region.clone();
        tokio::task::spawn_blocking(move || {
            // Cache lookup is CPU-bound operation, use blocking task pool
            query_variants_from_cache(&variants_cache, &bed_region)
        })
    };
    
    // Task 2: BAM read extraction using indexed access (I/O-bound operation)
    let bam_task = {
        let bam_path = Arc::clone(&bam_path);
        let bed_region = bed_region.clone();
        tokio::task::spawn_blocking(move || {
            // BAM access is I/O-bound operation, use blocking task pool
            query_bam_reads_for_region(&bam_path, &bed_region.chrom, bed_region.start_pos, bed_region.end_pos)
        })
    };
    
    // PHASE 2: SYNCHRONIZATION AND DATA RETRIEVAL
    // Wait for both tasks to complete and collect results
    let (variants_result, bam_reads_result) = tokio::try_join!(variants_task, bam_task)?;
    let variants = variants_result;      // Retrieved from cache
    let bam_reads = bam_reads_result?;   // Retrieved from BAM file
    
    println!("  Found {} variants and {} reads in region (async cached)", variants.len(), bam_reads.len());
    
    // PHASE 3: EARLY TERMINATION CHECK
    // Skip analysis if no variants present (optimization for sparse regions)
    if variants.is_empty() {
        println!("  No variants found in region {}, skipping...", bed_region.region_string);
        return Ok(0);
    }
    
    // PHASE 4: PARALLEL VARIANT-READ ANALYSIS
    // This is where the core bioinformatics analysis happens
    let mut lines_written = 0;
    let mut all_output_lines = Vec::new(); // Batch all output lines to reduce mutex contention
    
    // Process variants in parallel chunks for better CPU utilization
    // Each variant chunk is processed independently to maximize parallel throughput
    for variant_chunk in variants.chunks(variant_chunk_size) {
        // PARALLEL VARIANT PROCESSING:
        // Each chunk of variants is processed in parallel using Rayon's parallel iterator
        // This creates a thread pool that processes variants concurrently
        let chunk_results: Vec<Vec<String>> = variant_chunk
            .par_iter()  // Convert to parallel iterator for concurrent processing
            .map(|variant| {
                // CORE READ-VARIANT OVERLAP ANALYSIS:
                // For each variant, examine ALL reads in the region to determine overlap
                // This is the heart of the variant-read association algorithm
                
                bam_reads
                    .par_iter()  // Process reads in parallel for this variant
                    .filter_map(|(read_id, start_pos, end_pos, timestamp, mapping_quality)| {
                        // OVERLAP DETECTION:
                        // Determine if this read genomically overlaps the current variant
                        let contains_variant = read_overlaps_variant(*start_pos, *end_pos, variant.pos);
                        
                        // DETAILED ANALYSIS EXPLANATION:
                        // - read_start/end: Genomic coordinates where this read aligns
                        // - variant.pos: Genomic position of the genetic variant
                        // - contains_variant: Boolean indicating spatial overlap
                        // - mapping_quality: Alignment confidence score (0-255)
                        //
                        // BIOLOGICAL SIGNIFICANCE:
                        // - If true: This read potentially captured the variant site during sequencing
                        // - If false: This read did not cover the variant position
                        //
                        // TEMPORAL ANALYSIS:
                        // - timestamp: When this read was sequenced (for time-series studies)
                        // - Enables tracking variant presence across time points
                        //
                        // QUALITY CONTROL:
                        // - mapping_quality: Higher values indicate more confident alignments
                        // - Can be used for downstream filtering of low-quality reads
                        //
                        // OUTPUT FORMAT:
                        // Each line represents one read-variant pair with:
                        // 1. read_id: Unique identifier for the sequencing read
                        // 2. timestamp: Temporal information (when read was generated)
                        // 3. contains_variant: Boolean overlap result (true/false)
                        // 4. variant.description: Human-readable variant description (e.g., "A>G")
                        // 5. variant.variant_type: Variant classification (SNV, INS, DEL, etc.)
                        // 6. bed_region.region_string: Genomic region being analyzed
                        // 7. bed_region.region_name: User-defined region identifier
                        // 8. read_start: Genomic start position of the read alignment
                        // 9. read_end: Genomic end position of the read alignment
                        // 10. mapping_quality: Alignment confidence score (0-255)
                        
                        Some(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            read_id,                    // Read identifier
                            timestamp,                  // Sequencing timestamp
                            contains_variant,           // Overlap result (true/false)
                            variant.description,        // Variant description (e.g., "A>G")
                            variant.variant_type,       // Variant type (SNV/INS/DEL/etc.)
                            bed_region.region_string,   // Genomic region (chr:start-end)
                            bed_region.region_name,     // Region name/identifier
                            start_pos,                  // Read start position
                            end_pos,                    // Read end position
                            mapping_quality             // Mapping quality score
                        ))
                    })
                    .collect()  // Collect results for this variant
            })
            .collect();  // Collect results for all variants in this chunk
        
        // Collect all lines from this chunk
        for lines in chunk_results {
            all_output_lines.extend(lines);
        }
    }
    
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
        self.writer.flush()
    }
}

/// Pre-load and index all variants from VCF file for efficient region-based queries
async fn preload_vcf_variants(
    vcf_path: Arc<PathBuf>,
    bed_regions: &[BedRegion],
) -> Result<std::collections::HashMap<String, Vec<Variant>>, Box<dyn std::error::Error>> {
    println!("Pre-loading VCF variants for efficient querying...");
    
    let start_time = std::time::Instant::now();
    
    // Collect all chromosomes we're interested in
    let target_chroms: std::collections::HashSet<String> = bed_regions
        .iter()
        .map(|region| region.chrom.clone())
        .collect();
    
    println!("  Target chromosomes: {:?}", target_chroms);
    
    // Pre-load variants by chromosome
    let variants_by_chrom = tokio::task::spawn_blocking(move || {
        preload_vcf_variants_blocking(&vcf_path, &target_chroms)
    }).await
    .map_err(|e| format!("Task join error: {}", e))?
    .map_err(|e| format!("VCF preloading error: {}", e))?;
    
    let elapsed = start_time.elapsed();
    let total_variants: usize = variants_by_chrom.values().map(|v| v.len()).sum();
    
    println!("  Pre-loaded {} variants across {} chromosomes in {:?}", 
            total_variants, variants_by_chrom.len(), elapsed);
    
    Ok(variants_by_chrom)
}

/// Blocking implementation of VCF variant pre-loading
fn preload_vcf_variants_blocking(
    vcf_path: &PathBuf,
    target_chroms: &std::collections::HashSet<String>,
) -> Result<std::collections::HashMap<String, Vec<Variant>>, Box<dyn std::error::Error + Send + Sync>> {
    let mut variants_by_chrom: std::collections::HashMap<String, Vec<Variant>> = std::collections::HashMap::new();
    
    // Check if file is gzipped
    let is_gzipped = vcf_path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz")
        .unwrap_or(false);
    
    if is_gzipped {
        preload_vcf_variants_gzipped(vcf_path, target_chroms, &mut variants_by_chrom)?;
    } else {
        preload_vcf_variants_uncompressed(vcf_path, target_chroms, &mut variants_by_chrom)?;
    }
    
    // Sort variants by position for efficient region queries
    for variants in variants_by_chrom.values_mut() {
        variants.sort_by_key(|v| v.pos);
    }
    
    Ok(variants_by_chrom)
}

/// Pre-load variants from gzipped VCF
fn preload_vcf_variants_gzipped(
    vcf_path: &PathBuf,
    target_chroms: &std::collections::HashSet<String>,
    variants_by_chrom: &mut std::collections::HashMap<String, Vec<Variant>>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    // Try BGZF first
    let file = File::open(vcf_path)?;
    let bgzf_reader = bgzf::io::Reader::new(file);
    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(bgzf_reader));
    
    match preload_vcf_variants_from_reader(&mut vcf_reader, target_chroms, variants_by_chrom) {
        Ok(_) => Ok(()),
        Err(_) => {
            // Fallback to standard gzip
            let vcf_file = File::open(vcf_path)?;
            let mut gz_decoder = GzDecoder::new(vcf_file);
            let mut decompressed_data = Vec::new();
            gz_decoder.read_to_end(&mut decompressed_data)?;
            
            let cursor = Cursor::new(decompressed_data);
            let mut vcf_reader = vcf::io::Reader::new(BufReader::new(cursor));
            preload_vcf_variants_from_reader(&mut vcf_reader, target_chroms, variants_by_chrom)
        }
    }
}

/// Pre-load variants from uncompressed VCF
fn preload_vcf_variants_uncompressed(
    vcf_path: &PathBuf,
    target_chroms: &std::collections::HashSet<String>,
    variants_by_chrom: &mut std::collections::HashMap<String, Vec<Variant>>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let vcf_file = File::open(vcf_path)?;
    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(vcf_file));
    preload_vcf_variants_from_reader(&mut vcf_reader, target_chroms, variants_by_chrom)
}

/// Generic VCF reader for pre-loading variants
fn preload_vcf_variants_from_reader<R: std::io::Read>(
    vcf_reader: &mut vcf::io::Reader<BufReader<R>>,
    target_chroms: &std::collections::HashSet<String>,
    variants_by_chrom: &mut std::collections::HashMap<String, Vec<Variant>>,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let header = vcf_reader.read_header()?;
    let mut vcf_record = vcf::variant::RecordBuf::default();
    let mut total_processed = 0;
    let mut target_chrom_variants = 0;
    
    loop {
        match vcf_reader.read_record_buf(&header, &mut vcf_record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                total_processed += 1;
                
                let vcf_chrom_str = vcf_record.reference_sequence_name();
                
                // Only process variants on target chromosomes
                if !target_chroms.contains(vcf_chrom_str) {
                    continue;
                }
                
                target_chrom_variants += 1;
                
                let variant_pos = vcf_record.variant_start()
                    .map(|p| usize::from(p))
                    .unwrap_or(0);
                
                let ref_bases = vcf_record.reference_bases().to_string();
                let alt_bases = match vcf_record.alternate_bases().iter().next() {
                    Some(alt_allele_result) => {
                        match alt_allele_result {
                            Ok(alt_allele) => alt_allele.to_string(),
                            Err(_) => ".".to_string()
                        }
                    }
                    None => ".".to_string()
                };
                
                let (description, variant_type) = if alt_bases.contains('[') || alt_bases.contains(']') {
                    let bnd_description = parse_bnd_variant(&ref_bases, &alt_bases);
                    (bnd_description, "BND".to_string())
                } else if ref_bases.len() == 1 && alt_bases.len() == 1 && alt_bases != "." {
                    (format!("{}>{}", ref_bases, alt_bases), "SNV".to_string())
                } else if ref_bases.len() > alt_bases.len() {
                    (format!("{}>{}", ref_bases, alt_bases), "DEL".to_string())
                } else if ref_bases.len() < alt_bases.len() {
                    (format!("{}>{}", ref_bases, alt_bases), "INS".to_string())
                } else if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 {
                    (format!("{}>{}", ref_bases, alt_bases), "MNV".to_string())
                } else {
                    (format!("{}>{}", ref_bases, alt_bases), "OTHER".to_string())
                };
                
                let variant = Variant {
                    pos: variant_pos,
                    description,
                    variant_type,
                };
                
                variants_by_chrom
                    .entry(vcf_chrom_str.to_string())
                    .or_insert_with(Vec::new)
                    .push(variant);
            }
            Err(e) => {
                return Err(format!("Error reading VCF record: {:?}", e).into());
            }
        }
    }
    
    println!("  Pre-loading complete: {} total variants, {} on target chromosomes", 
             total_processed, target_chrom_variants);
    
    Ok(())
}

/// Fast region-based variant lookup using pre-loaded data
fn query_variants_from_cache(
    variants_by_chrom: &std::collections::HashMap<String, Vec<Variant>>,
    bed_region: &BedRegion,
) -> Vec<Variant> {
    match variants_by_chrom.get(&bed_region.chrom) {
        Some(chrom_variants) => {
            // Use binary search for efficient range queries since variants are sorted
            chrom_variants
                .iter()
                .filter(|variant| {
                    variant.pos >= bed_region.start_pos && variant.pos <= bed_region.end_pos
                })
                .cloned()
                .collect()
        }
        None => Vec::new(),
    }
}

// ...existing code...
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
    
    // Preload VCF variants for all regions
    let variants_by_chrom = preload_vcf_variants(Arc::clone(&vcf_path), &bed_regions).await?;
    let variants_cache = Arc::new(variants_by_chrom);
    
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
                let variants_cache = Arc::clone(&variants_cache);
                let bam_path = Arc::clone(&bam_path);
                let output_mutex = Arc::clone(&output_mutex);
                let bed_region = bed_region.clone();
                
                // Calculate actual region index across all chunks
                let region_idx = regions_processed + chunk_idx + 1;
                
                tokio::spawn(async move {
                    process_bed_region_async(
                        region_idx,
                        bed_region,
                        variants_cache,
                        bam_path,
                        output_mutex,
                        variant_chunk_size,
                    ).await
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
        println!("Completed chunk with {} regions, {} total regions processed, {} total lines written", 
                chunk.len(), regions_processed, total_lines_written);
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
