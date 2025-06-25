// Standard library imports for file I/O and path handling
use std::{fs::File, io::{BufReader, Write}, path::PathBuf, sync::Arc};

// Command-line argument parsing library
use clap::Parser;

// Async runtime and parallel processing
use tokio;
use rayon::prelude::*;

// Bioinformatics file format libraries from the noodles ecosystem
use noodles_bam as bam;  // Binary Alignment/Map format for sequencing reads
use noodles_bed as bed;  // Browser Extensible Data format for genomic regions
use noodles_vcf as vcf;  // Variant Call Format for genetic variants
use noodles_sam as sam;  // SAM format for accessing BAM auxiliary data

// Import for handling gzipped files
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
    ref_bases: String,
    alt_bases: String,
    description: String,
    variant_type: String,
}

/// Represents a BAM read with relevant information
#[derive(Debug, Clone)]
struct BamRead {
    read_id: String,
    timestamp: String,
    contains_variant: bool,
}

/// Preload all BAM reads to avoid re-reading the file for each variant
async fn preload_bam_reads(bam_path: &PathBuf) -> Result<Vec<(String, usize, usize, String)>, Box<dyn std::error::Error>> {
    println!("Preloading all BAM reads...");
    let bam_file = File::open(bam_path)?;
    let mut bam_reader = bam::io::Reader::new(BufReader::new(bam_file));
    
    // Read BAM header
    let _bam_header = bam_reader.read_header()?;
    
    let mut reads = Vec::new();
    let mut bam_record = bam::Record::default();
    let mut read_count = 0;
    
    while bam_reader.read_record(&mut bam_record)? != 0 {
        read_count += 1;
        if read_count % 100000 == 0 {
            println!("  Loaded {} BAM reads...", read_count);
        }
        
        // Extract read information
        let read_id = bam_record.name()
            .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
            .unwrap_or("unknown")
            .to_string();
        
        // Get alignment span
        let (start_pos, end_pos) = if let Some(Ok(alignment_start)) = bam_record.alignment_start() {
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
        
        reads.push((read_id, start_pos, end_pos, timestamp));
    }
    
    println!("Preloaded {} BAM reads", reads.len());
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
                    chrom,
                    start_pos,
                    end_pos,
                    region_name,
                    region_string,
                });
            }
            Err(e) => return Err(format!("Error reading BED record: {:?}", e).into()),
        }
    }
    
    println!("Loaded {} BED regions", regions.len());
    Ok(regions)
}

/// Process a single BED region: find variants and check which reads overlap
fn process_bed_region(
    region_idx: usize,
    bed_region: &BedRegion,
    vcf_path: &PathBuf,
    bam_reads: &Arc<Vec<(String, usize, usize, String)>>,
) -> Result<Vec<String>, Box<dyn std::error::Error + Send + Sync>> {
    println!("Processing BED region #{}: {}", region_idx, bed_region.region_string);
    
    // Load variants from VCF for this region
    let variants = load_variants_for_region(vcf_path, bed_region)?;
    println!("  Found {} variants in region", variants.len());
    
    let mut output_lines = Vec::new();
    
    // Process each variant in parallel
    let variant_results: Vec<_> = variants
        .par_iter()
        .flat_map(|variant| {
            // Find reads that overlap with this variant
            let overlapping_reads: Vec<_> = bam_reads
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
            
            overlapping_reads
        })
        .collect();
    
    output_lines.extend(variant_results);
    
    println!("  Processed region #{} with {} output lines", region_idx, output_lines.len());
    Ok(output_lines)
}

/// Load variants from VCF that fall within a specific BED region
fn load_variants_for_region(
    vcf_path: &PathBuf,
    bed_region: &BedRegion,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    let vcf_file = File::open(vcf_path)?;
    let gz_decoder = GzDecoder::new(vcf_file);
    let mut vcf_reader = vcf::io::Reader::new(BufReader::new(gz_decoder));
    
    let header = vcf_reader.read_header()?;
    let mut vcf_record = vcf::variant::RecordBuf::default();
    let mut variants = Vec::new();
    
    loop {
        match vcf_reader.read_record_buf(&header, &mut vcf_record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                let variant_pos = vcf_record.variant_start()
                    .map(|p| usize::from(p))
                    .unwrap_or(0);
                
                // Check if variant is in this region
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
                    
                    variants.push(Variant {
                        pos: variant_pos,
                        ref_bases,
                        alt_bases,
                        description,
                        variant_type,
                    });
                }
            }
            Err(e) => return Err(format!("Error reading VCF record: {:?}", e).into()),
        }
    }
    
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
#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
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

    // ===== PRELOAD BAM READS =====
    // Load all BAM reads once to avoid re-reading for each variant
    let bam_reads = preload_bam_reads(&args.bam).await?;
    let bam_reads = Arc::new(bam_reads); // Share between threads

    // ===== BED FILE PROCESSING =====
    // Load all BED regions first
    let bed_regions = load_bed_regions(&args.bed)?;
    
    println!("\n=== Processing {} BED regions in parallel ===", bed_regions.len());
    
    // Process BED regions in parallel
    let results: Vec<_> = bed_regions
        .par_iter()
        .enumerate()
        .map(|(region_idx, bed_region)| {
            process_bed_region(region_idx + 1, bed_region, &args.vcf, &bam_reads)
        })
        .collect();

    // Write all results to output file
    for result in results {
        match result {
            Ok(region_results) => {
                for line in region_results {
                    writeln!(output, "{}", line)?;
                }
            }
            Err(e) => {
                eprintln!("Error processing region: {:?}", e);
            }
        }
    }

    // ===== BED FILE PROCESSING =====
    // Open and parse the BED file containing genomic regions of interest
    // BED files define genomic intervals (chromosome, start, end positions)
    
    println!("\n=== Processing BED file ===");
    
    // Check if BED file exists and get its size
    let bed_metadata = match std::fs::metadata(&args.bed) {
        Ok(metadata) => metadata,
        Err(e) => {
            println!("ERROR: Cannot access BED file: {:?}", e);
            return Err(format!("Cannot access BED file: {:?}", e).into());
        }
    };
    println!("BED file size: {} bytes", bed_metadata.len());
    if bed_metadata.len() == 0 {
        return Err("BED file is empty! Please provide a valid BED file with genomic regions.".into());
    }
    
    let bed_file = match File::open(&args.bed) {
        Ok(file) => file,
        Err(e) => {
            println!("ERROR: Failed to open BED file: {:?}", e);
            return Err(format!("Failed to open BED file: {:?}", e).into());
        }
    };
    println!("Successfully opened BED file: {:?}", args.bed);
    
    // Use BufReader for efficient reading of potentially large files
    // Specify BED4 format (4 columns: chr, start, end, name) to match your file format
    let mut bed_reader = bed::io::Reader::<4, _>::new(BufReader::new(bed_file));
    
    // Create a default BED record to reuse for memory efficiency
    // BED4 records contain: reference sequence name, start position, end position, name
    let mut bed_record = bed::Record::<4>::default();
    
    let mut bed_region_count = 0;
    
    // Process each BED record (genomic region) in the file
    // read_record() returns the number of bytes read (0 when EOF reached)
    println!("Starting to read BED records...");
    
    loop {
        match bed_reader.read_record(&mut bed_record) {
            Ok(0) => {
                // EOF reached
                println!("Reached end of BED file");
                break;
            }
            Ok(bytes_read) => {
                bed_region_count += 1;
                println!("Processing BED region #{} (read {} bytes)", bed_region_count, bytes_read);
            }
            Err(e) => {
                println!("ERROR reading BED record: {:?}", e);
                return Err(format!("Failed to read BED record #{}: {:?}", bed_region_count + 1, e).into());
            }
        }
        // Extract chromosome/contig name from the BED record
        // Convert from bytes to UTF-8 string with error handling
        let chrom = match std::str::from_utf8(bed_record.reference_sequence_name()) {
            Ok(s) => s.to_string(),
            Err(e) => {
                println!("ERROR: Invalid UTF-8 in chromosome name: {:?}", e);
                return Err(format!("Invalid chromosome name in BED record #{}: {:?}", bed_region_count, e).into());
            }
        };
        
        // Extract start position (0-indexed in BED format)
        // Returns Result<Position, Error> which we unwrap with ?
        let start = match bed_record.feature_start() {
            Ok(pos) => pos,
            Err(e) => {
                println!("ERROR: Invalid start position in BED record #{}: {:?}", bed_region_count, e);
                return Err(format!("Invalid start position in BED record #{}: {:?}", bed_region_count, e).into());
            }
        };
        
        // Extract end position (exclusive in BED format - half-open interval [start, end))
        // feature_end() returns Option<Result<Position, Error>>
        // .ok_or() converts None to an error, ?? unwraps both Option and Result
        let end = match bed_record.feature_end() {
            Some(Ok(pos)) => pos,
            Some(Err(e)) => {
                println!("ERROR: Invalid end position in BED record #{}: {:?}", bed_region_count, e);
                return Err(format!("Invalid end position in BED record #{}: {:?}", bed_region_count, e).into());
            }
            None => {
                println!("ERROR: Missing end position in BED record #{}", bed_region_count);
                return Err(format!("Missing end position in BED record #{}", bed_region_count).into());
            }
        };
            
        // Extract the region name from the BED4 name column (4th column)
        let region_name = if let Some(name_bytes) = bed_record.name() {
            std::str::from_utf8(name_bytes).unwrap_or("INVALID_UTF8")
        } else {
            "NA" // Fallback if name is missing
        };
        
        // Convert genomic positions to usize for easier arithmetic
        // noodles Position types implement Into<usize>
        let start_pos = usize::from(start);
        let end_pos = usize::from(end);
        
        // Create human-readable string representation of the genomic region
        // Format: "chromosome:start-end" (converting to 1-indexed for display)
        let region_string = format!("{}:{}-{}", chrom, start_pos, end_pos);
        
        println!("  Region: {} (chr: {}, start: {}, end: {}, name: {})", 
                 region_string, chrom, start_pos, end_pos, region_name);

        // ===== VCF FILE PROCESSING =====
        // For each BED region, scan the entire VCF file to find overlapping variants
        // Note: This is inefficient for large files - ideally we'd use indexed access
        
        println!("  Opening VCF file: {:?}", args.vcf);
        let vcf_file = File::open(&args.vcf)?;
        
        // Expect gzipped VCF files for production use
        // Check if the VCF file has .gz extension
        if !args.vcf.extension().map(|ext| ext == "gz").unwrap_or(false) {
            println!("WARNING: VCF file does not have .gz extension. This tool expects gzipped VCF files.");
            println!("If your VCF is not gzipped, please compress it first with: bgzip your_file.vcf");
        }
        
        // Create VCF reader for gzipped files
        let gz_decoder = GzDecoder::new(vcf_file);
        let mut vcf_reader = vcf::io::Reader::new(BufReader::new(gz_decoder));
        println!("  Reading VCF header...");
        let header = match vcf_reader.read_header() {
            Ok(h) => {
                println!("  VCF header read successfully");
                h
            }
            Err(e) => {
                println!("ERROR reading VCF header: {:?}", e);
                return Err(format!("Failed to read VCF header: {:?}", e).into());
            }
        };
        
        // Create reusable VCF record buffer for memory efficiency
        // RecordBuf can store variant information: position, ref allele, alt alleles, etc.
        let mut vcf_record = vcf::variant::RecordBuf::default();
        
        let mut vcf_variant_count = 0;
        let mut matching_variants = 0;
        
        // Process each variant record in the VCF file
        // Re-open the gzipped file for record reading
        let vcf_file_for_records = File::open(&args.vcf)?;
        let gz_decoder = GzDecoder::new(vcf_file_for_records);
        let mut vcf_reader = vcf::io::Reader::new(BufReader::new(gz_decoder));
        
        // Skip the header since we already read it
        let _header = vcf_reader.read_header()?;
        
        loop {
            match vcf_reader.read_record_buf(&header, &mut vcf_record) {
                Ok(0) => {
                    // EOF reached
                    if vcf_variant_count <= 5 {
                        println!("    Reached end of VCF file");
                    }
                    break;
                }
                Ok(bytes_read) => {
                    vcf_variant_count += 1;
                    if vcf_variant_count <= 5 || vcf_variant_count % 1000 == 0 {
                        println!("    Processing VCF variant #{} (read {} bytes)", vcf_variant_count, bytes_read);
                    }
                }
                Err(e) => {
                    println!("ERROR reading VCF record: {:?}", e);
                    return Err(format!("Failed to read VCF record #{}: {:?}", vcf_variant_count + 1, e).into());
                }
            }
            
            // Extract the genomic position where this variant occurs
            let variant_pos = vcf_record.variant_start()
                .map(|p| usize::from(p))
                .unwrap_or(0);
            
            // Check if this variant falls within the current BED region of interest
            if variant_pos >= start_pos && variant_pos <= end_pos {
                matching_variants += 1;
                
                // Extract the reference allele
                let ref_bases = vcf_record.reference_bases().to_string();
                
                // Extract alternate alleles
                let alt_bases = match vcf_record.alternate_bases().iter().next() {
                    Some(alt_allele_result) => {
                        match alt_allele_result {
                            Ok(alt_allele) => alt_allele.to_string(),
                            Err(_) => ".".to_string()
                        }
                    }
                    None => ".".to_string()
                };
                
                // Create variant description
                let variant_description = format!("{}>{}", ref_bases, alt_bases);
                
                // Determine variant type
                let variant_type = if ref_bases.len() == 1 && alt_bases.len() == 1 && alt_bases != "." {
                    "SNV"
                } else if ref_bases.len() > alt_bases.len() {
                    "DEL"
                } else if ref_bases.len() < alt_bases.len() {
                    "INS"
                } else if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 {
                    "MNV"
                } else {
                    "OTHER"
                };
                
                println!("    Found matching variant #{}: {} at position {} ({})", 
                         matching_variants, variant_description, variant_pos, region_string);

                // BAM file processing
                println!("    Opening BAM file: {:?}", args.bam);
                let bam_file = File::open(&args.bam)?;
                let mut bam_reader = bam::io::Reader::new(BufReader::new(bam_file));
                
                println!("    Reading BAM header...");
                let _bam_header = bam_reader.read_header()?;
                println!("    BAM header read successfully");
                
                let mut bam_record = bam::Record::default();
                let mut bam_read_count = 0;
                let mut variant_supporting_reads = 0;
                
                // Process each sequencing read in the BAM file
                while bam_reader.read_record(&mut bam_record)? != 0 {
                    bam_read_count += 1;
                    if bam_read_count <= 5 || bam_read_count % 10000 == 0 {
                        println!("      Processing BAM read #{}", bam_read_count);
                    }
                    
                    // Extract read ID
                    let read_id = bam_record.name()
                        .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
                        .unwrap_or("unknown");
                    
                    // ===== VARIANT OVERLAP DETECTION =====
                    // Use CIGAR-based approach for accurate overlap detection
                    let contains_variant = if let Some(Ok(alignment_start)) = bam_record.alignment_start() {
                        let start_pos = usize::from(alignment_start);
                        
                        // Parse CIGAR string to determine actual reference span
                        let cigar = bam_record.cigar();
                        let mut reference_pos = start_pos;
                        let mut covers_variant = false;
                        
                        // Iterate through each CIGAR operation
                        for operation in cigar.iter() {
                            match operation {
                                Ok(op) => {
                                    use noodles_sam::alignment::record::cigar::op::Kind;
                                    let op_len = usize::from(op.len());
                                    
                                    match op.kind() {
                                        // Operations that consume both read and reference
                                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                                            // Check if variant falls within this operation's span
                                            if variant_pos >= reference_pos && variant_pos < reference_pos + op_len {
                                                covers_variant = true;
                                                break;
                                            }
                                            reference_pos += op_len;
                                        }
                                        // Operations that consume reference but not read
                                        Kind::Deletion | Kind::Skip => {
                                            reference_pos += op_len;
                                        }
                                        // Operations that consume read but not reference
                                        Kind::Insertion | Kind::SoftClip => {
                                            continue;
                                        }
                                        // Operations that don't consume either
                                        Kind::HardClip | Kind::Pad => {
                                            continue;
                                        }
                                    }
                                }
                                Err(_) => {
                                    // If we can't parse a CIGAR operation, skip this read
                                    break;
                                }
                            }
                        }
                        
                        covers_variant
                    } else {
                        false
                    };

                    // Extract timestamp from BAM auxiliary data
                    let timestamp = {
                        let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
                        match bam_record.data().get(&st_tag) {
                            Some(Ok(field_value)) => {
                                match field_value {
                                    sam::alignment::record::data::field::Value::String(ts_bytes) => {
                                        std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                                    }
                                    sam::alignment::record::data::field::Value::Int8(val) => {
                                        format!("INT8_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt8(val) => {
                                        format!("UINT8_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Int16(val) => {
                                        format!("INT16_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt16(val) => {
                                        format!("UINT16_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Int32(val) => {
                                        format!("INT32_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt32(val) => {
                                        format!("UINT32_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Float(val) => {
                                        format!("FLOAT_{}", val)
                                    }
                                    _ => {
                                        format!("UNKNOWN_TYPE_{:?}", field_value)
                                    }
                                }
                            }
                            Some(Err(e)) => {
                                format!("PARSE_ERROR_{:?}", e)
                            }
                            None => {
                                "NA".to_string()
                            }
                        }
                    };

                    // Output generation
                    if contains_variant {
                        variant_supporting_reads += 1;
                        if variant_supporting_reads <= 5 {
                            println!("      Read {} supports variant {} (timestamp: {})", 
                                   read_id, variant_description, timestamp);
                        }
                    }
                    
                    writeln!(
                        output,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_id,
                        timestamp,
                        contains_variant,
                        variant_description,
                        variant_type,
                        region_string,
                        region_name
                    )?;
                }
                
                println!("    Processed {} BAM reads for variant {}, {} support the variant", 
                         bam_read_count, variant_description, variant_supporting_reads);
            }
        }
        
        println!("  Processed {} VCF variants, {} matched the BED region", vcf_variant_count, matching_variants);
    }
    
    println!("\n=== Analysis Complete ===");
    println!("Processed {} BED regions", bed_region_count);
    if bed_region_count == 0 {
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
