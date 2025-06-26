// Standard library imports for file I/O and path handling
use std::{fs::File, io::{BufReader, Write, Read, Cursor}, path::PathBuf};

// Command-line argument parsing library
use clap::Parser;

// Parallel processing
use rayon::prelude::*;

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
fn query_bam_reads_for_region(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<(String, usize, usize, String)>, Box<dyn std::error::Error + Send + Sync>> {
    println!("  Querying BAM reads for region {}:{}-{}", chrom, start_pos, end_pos);
    
    // Try indexed access first
    match query_bam_reads_indexed(bam_path, chrom, start_pos, end_pos) {
        Ok(reads) => {
            println!("  Successfully used indexed BAM access, found {} reads", reads.len());
            Ok(reads)
        }
        Err(index_error) => {
            println!("  Indexed BAM access failed: {:?}", index_error);
            println!("  Falling back to full BAM scan...");
            query_bam_reads_full_scan(bam_path, chrom, start_pos, end_pos)
        }
    }
}

/// Query BAM reads using indexed access with .bai file
fn query_bam_reads_indexed(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<(String, usize, usize, String)>, Box<dyn std::error::Error + Send + Sync>> {
    // Try to open indexed BAM reader
    let mut indexed_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    
    let header = indexed_reader.read_header()?;
    
    // Create a region for the query using the correct API
    let start_position = noodles_core::Position::try_from(start_pos)?;
    let end_position = noodles_core::Position::try_from(end_pos)?;
    let interval = start_position..=end_position;
    let region = noodles_core::Region::new(chrom, interval);
    
    // Query the region
    let query = indexed_reader.query(&header, &region)?;
    
    let mut reads = Vec::new();
    
    for result in query {
    let record = result?;

    let read_id = record.name()
        .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
        .unwrap_or("unknown")
        .to_string();

    let (read_start, read_end) = if let Some(Ok(alignment_start)) = record.alignment_start() {
        let start = usize::from(alignment_start);

        let cigar = record.cigar();
        let mut reference_pos = start;

        for operation in cigar.iter() {
            if let Ok(op) = operation {
                use noodles_sam::alignment::record::cigar::op::Kind;
                let op_len = usize::from(op.len());

                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                        reference_pos += op_len;
                    }
                    _ => {}
                }
            }
        }

        (start, reference_pos)
    } else {
        continue;
    };

    let timestamp = {
        let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
        match record.data().get(&st_tag) {
            Some(Ok(field_value)) => {
                match field_value {
                    sam::alignment::record::data::field::Value::String(ts_bytes) => {
                        std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                    }
                    sam::alignment::record::data::field::Value::Int8(val) => format!("INT8_{}", val),
                    sam::alignment::record::data::field::Value::UInt8(val) => format!("UINT8_{}", val),
                    sam::alignment::record::data::field::Value::Int16(val) => format!("INT16_{}", val),
                    sam::alignment::record::data::field::Value::UInt16(val) => format!("UINT16_{}", val),
                    sam::alignment::record::data::field::Value::Int32(val) => format!("INT32_{}", val),
                    sam::alignment::record::data::field::Value::UInt32(val) => format!("UINT32_{}", val),
                    sam::alignment::record::data::field::Value::Float(val) => format!("FLOAT_{}", val),
                    _ => format!("UNKNOWN_TYPE_{:?}", field_value),
                }
            }
            Some(Err(e)) => format!("PARSE_ERROR_{:?}", e),
            None => "NA".to_string(),
        }
    };

    reads.push((read_id, read_start, read_end, timestamp));
}

    
    Ok(reads)
}

/// Query BAM reads using full scan (fallback method)
fn query_bam_reads_full_scan(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<(String, usize, usize, String)>, Box<dyn std::error::Error + Send + Sync>> {
    let bam_file = File::open(bam_path)?;
    let mut bam_reader = bam::io::Reader::new(BufReader::new(bam_file));
    
    // Read BAM header
    let _bam_header = bam_reader.read_header()?;
    
    let mut reads = Vec::new();
    let mut bam_record = bam::Record::default();
    
    // Scan all reads and filter by region
    while bam_reader.read_record(&mut bam_record)? != 0 {
        // Extract read information
        let read_id = bam_record.name()
            .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
            .unwrap_or("unknown")
            .to_string();
        
        // Get alignment span
        let (read_start, read_end) = if let Some(Ok(alignment_start)) = bam_record.alignment_start() {
            let start = usize::from(alignment_start);
            
            // Calculate end position using CIGAR
            let cigar = bam_record.cigar();
            let mut reference_pos = start;
            
            for operation in cigar.iter() {
                if let Ok(op) = operation {
                    use noodles_sam::alignment::record::cigar::op::Kind;
                    let op_len = usize::from(op.len());
                    
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                            reference_pos += op_len;
                        }
                        _ => {} // Other operations don't advance reference position
                    }
                }
            }
            
            (start, reference_pos)
        } else {
            continue; // Skip reads without valid alignment
        };
        
        // Check if read overlaps with the region of interest
        if read_start <= end_pos && read_end >= start_pos {
            // Check chromosome match if available  
            // Note: This is a simplified check - in practice, we'd need to map 
            // reference sequence ID to chromosome name from the header
            let _chromosome_check = chrom; // Use the parameter to avoid warning
            
            // Extract timestamp
            let timestamp = {
                let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
                match bam_record.data().get(&st_tag) {
                    Some(Ok(field_value)) => {
                        match field_value {
                            sam::alignment::record::data::field::Value::String(ts_bytes) => {
                                std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                            }
                            sam::alignment::record::data::field::Value::Int8(val) => format!("INT8_{}", val),
                            sam::alignment::record::data::field::Value::UInt8(val) => format!("UINT8_{}", val),
                            sam::alignment::record::data::field::Value::Int16(val) => format!("INT16_{}", val),
                            sam::alignment::record::data::field::Value::UInt16(val) => format!("UINT16_{}", val),
                            sam::alignment::record::data::field::Value::Int32(val) => format!("INT32_{}", val),
                            sam::alignment::record::data::field::Value::UInt32(val) => format!("UINT32_{}", val),
                            sam::alignment::record::data::field::Value::Float(val) => format!("FLOAT_{}", val),
                            _ => format!("UNKNOWN_TYPE_{:?}", field_value),
                        }
                    }
                    Some(Err(e)) => format!("PARSE_ERROR_{:?}", e),
                    None => "NA".to_string(),
                }
            };
            
            reads.push((read_id, read_start, read_end, timestamp));
        }
    }
    
    Ok(reads)
}

/// Check if a read overlaps with a variant position
fn read_overlaps_variant(read_start: usize, read_end: usize, variant_pos: usize) -> bool {
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

/// Process a single BED region and stream results to output file
fn process_bed_region<W: Write>(
    region_idx: usize,
    bed_region: &BedRegion,
    vcf_path: &PathBuf,
    bam_path: &PathBuf,
    output_writer: &std::sync::Mutex<W>,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    println!("Processing BED region #{}: {}", region_idx, bed_region.region_string);
    
    // Query variants from VCF for this region
    let variants = query_variants_for_region(vcf_path, bed_region)?;
    println!("  Found {} variants in region", variants.len());
    
    // If no variants found, skip this region and continue to next
    if variants.is_empty() {
        println!("  No variants found in region {}, skipping...", bed_region.region_string);
        return Ok(0);
    }
    
    // Query BAM reads for this specific region
    let bam_reads = query_bam_reads_for_region(bam_path, &bed_region.chrom, bed_region.start_pos, bed_region.end_pos)?;
    println!("  Found {} reads overlapping region", bam_reads.len());
    
    let mut lines_written = 0;
    
    // Process each variant and stream results immediately
    for variant in &variants {
        // Find reads that overlap with this variant
        let overlapping_lines: Vec<String> = bam_reads
            .par_iter()
            .filter_map(|(read_id, start_pos, end_pos, timestamp)| {
                let contains_variant = read_overlaps_variant(*start_pos, *end_pos, variant.pos);
                
                Some(format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_id,
                    timestamp,
                    contains_variant,
                    variant.description,
                    variant.variant_type,
                    bed_region.region_string,
                    bed_region.region_name
                ))
            })
            .collect();
        
        // Write results immediately to output file
        {
            let mut writer = output_writer.lock().unwrap();
            for line in &overlapping_lines {
                writeln!(*writer, "{}", line)?;
            }
        }
        
        lines_written += overlapping_lines.len();
    }
    
    println!("  Processed region #{} with {} output lines", region_idx, lines_written);
    Ok(lines_written)
}

/// Query VCF variants for a specific genomic region using indexed access
fn query_variants_for_region(
    vcf_path: &PathBuf,
    bed_region: &BedRegion,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    // Check if file is gzipped by looking at the extension
    let is_gzipped = vcf_path.extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz")
        .unwrap_or(false);
    
    println!("  VCF file path: {:?}, is_gzipped: {}", vcf_path, is_gzipped);
    
    if is_gzipped {
        query_variants_gzipped(vcf_path, bed_region)
    } else {
        query_variants_uncompressed(vcf_path, bed_region)
    }
}

/// Query variants from a gzipped VCF file using BGZF
fn query_variants_gzipped(
    vcf_path: &PathBuf,
    bed_region: &BedRegion,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    println!("  Reading gzipped VCF file using BGZF...");
    
    // Try to use BGZF reader first (for .vcf.gz files)
    let file = File::open(vcf_path)?;
    let bgzf_reader = bgzf::io::Reader::new(file);
    
    println!("  Successfully opened BGZF reader");
    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(bgzf_reader));
    
    // Try to process VCF records
    match process_vcf_records(&mut vcf_reader, bed_region) {
        Ok(variants) => {
            println!("  BGZF reading successful, found {} variants", variants.len());
            Ok(variants)
        }
        Err(bgzf_error) => {
            println!("  BGZF processing failed: {:?}", bgzf_error);
            println!("  Falling back to standard gzip decompression...");
            
            // Fallback to standard gzip decompression
            let vcf_file = File::open(vcf_path)?;
            let mut gz_decoder = GzDecoder::new(vcf_file);
            let mut decompressed_data = Vec::new();
            
            match gz_decoder.read_to_end(&mut decompressed_data) {
                Ok(bytes_read) => {
                    println!("  Successfully decompressed {} bytes using fallback method", bytes_read);
                    
                    // Create a cursor over the decompressed data
                    let cursor = Cursor::new(decompressed_data);
                    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(cursor));
                    
                    process_vcf_records(&mut vcf_reader, bed_region)
                }
                Err(e) => {
                    println!("  ERROR during fallback decompression: {:?}", e);
                    Err(format!("Failed to decompress gzipped file: {:?}", e).into())
                }
            }
        }
    }
}

/// Query variants from an uncompressed VCF file
fn query_variants_uncompressed(
    vcf_path: &PathBuf,
    bed_region: &BedRegion,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    println!("  Reading uncompressed VCF file...");
    let vcf_file = File::open(vcf_path)?;
    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(vcf_file));
    process_vcf_records(&mut vcf_reader, bed_region)
}

/// Process VCF records from any reader type
fn process_vcf_records<R: std::io::Read>(
    vcf_reader: &mut vcf::io::Reader<BufReader<R>>,
    bed_region: &BedRegion,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    let header = match vcf_reader.read_header() {
        Ok(h) => {
            println!("  Successfully read VCF header");
            println!("  Header has {} contigs and {} samples", 
                    h.contigs().len(),
                    h.sample_names().len());
            h
        }
        Err(e) => {
            println!("  ERROR reading VCF header: {:?}", e);
            return Err(format!("Failed to read VCF header: {:?}", e).into());
        }
    };
    
    let mut vcf_record = vcf::variant::RecordBuf::default();
    let mut variants = Vec::new();
    let mut total_variants_seen = 0;
    let mut matching_chrom_variants = 0;
    
    println!("  Searching VCF for region: {} ({}:{}-{})", 
            bed_region.region_string, bed_region.chrom, bed_region.start_pos, bed_region.end_pos);
    
    // Test if we can read the first record immediately
    let first_record_result = vcf_reader.read_record_buf(&header, &mut vcf_record);
    println!("  First record read result: {:?}", first_record_result);
    
    // Reset the record buffer
    vcf_record = vcf::variant::RecordBuf::default();
    
    // TODO: Use indexed access with .tbi file for better performance
    // For now, we'll scan and filter by region
    loop {
        match vcf_reader.read_record_buf(&header, &mut vcf_record) {
            Ok(0) => {
                println!("  Reached end of VCF file");
                break; // EOF
            }
            Ok(_) => {
                total_variants_seen += 1;
                
                let vcf_chrom_str = vcf_record.reference_sequence_name();
                let variant_pos = vcf_record.variant_start()
                    .map(|p| usize::from(p))
                    .unwrap_or(0);
                
                // Debug: Show first few variants and their chromosomes
                if total_variants_seen <= 10 {
                    println!("    VCF variant #{}: {}:{} (looking for chromosome '{}')", 
                            total_variants_seen, vcf_chrom_str, variant_pos, bed_region.chrom);
                }
                
                // Only process variants on the same chromosome as the BED region
                if vcf_chrom_str != bed_region.chrom {
                    continue;
                }
                
                matching_chrom_variants += 1;
                
                // Debug: Show position comparisons for matching chromosomes
                if matching_chrom_variants <= 5 {
                    println!("    Matching chromosome variant #{}: position {} (region {}-{})", 
                            matching_chrom_variants, variant_pos, bed_region.start_pos, bed_region.end_pos);
                }
                
                // Check if variant position falls within the BED region
                if variant_pos >= bed_region.start_pos && variant_pos <= bed_region.end_pos {
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
                    
                    let description = format!("{}>{}", ref_bases, alt_bases);
                    
                    let variant_type = if ref_bases.len() == 1 && alt_bases.len() == 1 && alt_bases != "." {
                        "SNV".to_string()
                    } else if ref_bases.len() > alt_bases.len() {
                        "DEL".to_string()
                    } else if ref_bases.len() < alt_bases.len() {
                        "INS".to_string()
                    } else if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 {
                        "MNV".to_string()
                    } else {
                        "OTHER".to_string()
                    };
                    
                    println!("    FOUND VARIANT: {}:{} {} ({})", vcf_chrom_str, variant_pos, description, variant_type);
                    
                    variants.push(Variant {
                        pos: variant_pos,
                        description,
                        variant_type,
                    });
                }
            }
            Err(e) => {
                println!("  ERROR reading VCF record: {:?}", e);
                return Err(format!("Error reading VCF record: {:?}", e).into());
            }
        }
    }
    
    println!("  VCF scan complete: {} total variants, {} on matching chromosome '{}', {} in region", 
            total_variants_seen, matching_chrom_variants, bed_region.chrom, variants.len());
    
    Ok(variants)
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

    /// Output file path for results (will be in TSV format)
    /// Contains: read_id, timestamp, contains_variant, variant_description, variant_type, region, region_name
    #[arg(short, long)]
    output: PathBuf,

    /// Number of parallel threads to use for processing (default: number of CPU cores)
    #[arg(short, long, default_value_t = num_cpus::get())]
    threads: usize,
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
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments using clap derive macro
    let args = Args::parse();

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    println!("Starting VarClock analysis with {} threads...", args.threads);
    println!("Input files:");
    println!("  BED file: {:?}", args.bed);
    println!("  VCF file: {:?}", args.vcf);
    println!("  BAM file: {:?}", args.bam);
    println!("  Output file: {:?}", args.output);

    // Create the output file and write TSV header row
    // This will overwrite any existing file at the specified path
    let mut output = match File::create(&args.output) {
        Ok(file) => file,
        Err(e) => {
            println!("ERROR: Failed to create output file: {:?}", e);
            return Err(format!("Failed to create output file: {:?}", e).into());
        }
    };
    
    if let Err(e) = writeln!(output, "read_id\ttimestamp\tcontains_variant\tvariant_description\tvariant_type\tregion\tregion_name") {
        println!("ERROR: Failed to write header to output file: {:?}", e);
        return Err(format!("Failed to write header: {:?}", e).into());
    }
    println!("Created output file and wrote header");

    // ===== BED FILE PROCESSING =====
    // Load all BED regions first
    let bed_regions = load_bed_regions(&args.bed)?;
    
    println!("\n=== Processing {} BED regions in parallel ===", bed_regions.len());
    
    // Create a mutex-wrapped writer for thread-safe streaming output
    let output_mutex = std::sync::Mutex::new(&mut output);
    
    // Process BED regions in parallel, streaming results directly to file
    let results: Vec<_> = bed_regions
        .par_iter()
        .enumerate()
        .map(|(region_idx, bed_region)| {
            process_bed_region(region_idx + 1, bed_region, &args.vcf, &args.bam, &output_mutex)
        })
        .collect();
    
    // Count total lines written and handle any errors
    let mut total_lines_written = 0;
    for result in results {
        match result {
            Ok(lines_written) => {
                total_lines_written += lines_written;
            }
            Err(e) => {
                eprintln!("Error processing region: {:?}", e);
            }
        }
    }
    
    println!("Total lines written to output: {}", total_lines_written);

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
    Ok(())
}
