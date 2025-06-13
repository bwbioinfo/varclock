// Standard library imports for file I/O and path handling
use std::{fs::File, io::{BufReader, Write}, path::PathBuf};

// Command-line argument parsing library
use clap::Parser;

// Bioinformatics file format libraries from the noodles ecosystem
use noodles_bam as bam;  // Binary Alignment/Map format for sequencing reads
use noodles_bed as bed;  // Browser Extensible Data format for genomic regions
use noodles_vcf as vcf;  // Variant Call Format for genetic variants
use noodles_sam as sam;  // SAM format for accessing BAM auxiliary data

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

    // Create the output file and write TSV header row
    // This will overwrite any existing file at the specified path
    let mut output = File::create(&args.output)?;
    writeln!(output, "read_id\ttimestamp\tcontains_variant\tvariant_description\tvariant_type\tregion\tregion_name")?;

    // ===== BED FILE PROCESSING =====
    // Open and parse the BED file containing genomic regions of interest
    // BED files define genomic intervals (chromosome, start, end positions)
    
    let bed_file = File::open(&args.bed)?;
    // Use BufReader for efficient reading of potentially large files
    // Specify BED3 format (3 columns: chr, start, end) - most basic BED format
    let mut bed_reader = bed::io::Reader::<3, _>::new(BufReader::new(bed_file));
    
    // Create a default BED record to reuse for memory efficiency
    // BED3 records contain: reference sequence name, start position, end position
    let mut bed_record = bed::Record::<3>::default();
    
    // Process each BED record (genomic region) in the file
    // read_record() returns the number of bytes read (0 when EOF reached)
    while bed_reader.read_record(&mut bed_record)? != 0 {
        // Extract chromosome/contig name from the BED record
        // Convert from bytes to UTF-8 string with error handling
        let chrom = std::str::from_utf8(bed_record.reference_sequence_name())?.to_string();
        
        // Extract start position (0-indexed in BED format)
        // Returns Result<Position, Error> which we unwrap with ?
        let start = bed_record.feature_start()?;
        
        // Extract end position (exclusive in BED format - half-open interval [start, end))
        // feature_end() returns Option<Result<Position, Error>>
        // .ok_or() converts None to an error, ?? unwraps both Option and Result
        let end = bed_record.feature_end()
            .ok_or("Missing end position")??;
            
        // BED3 format doesn't include names, so we use a placeholder
        let region_name = "NA"; // BED3 doesn't have names
        
        // Convert genomic positions to usize for easier arithmetic
        // noodles Position types implement Into<usize>
        let start_pos = usize::from(start);
        let end_pos = usize::from(end);
        
        // Create human-readable string representation of the genomic region
        // Format: "chromosome:start-end" (converting to 1-indexed for display)
        let region_string = format!("{}:{}-{}", chrom, start_pos, end_pos);

        // ===== VCF FILE PROCESSING =====
        // For each BED region, scan the entire VCF file to find overlapping variants
        // Note: This is inefficient for large files - ideally we'd use indexed access
        
        let vcf_file = File::open(&args.vcf)?;
        // Create VCF reader with buffered I/O for performance
        let mut vcf_reader = vcf::io::Reader::new(BufReader::new(vcf_file));
        
        // Read and parse the VCF header which contains metadata about the file
        // Header includes format version, sample information, and field definitions
        let header = vcf_reader.read_header()?;
        
        // Create reusable VCF record buffer for memory efficiency
        // RecordBuf can store variant information: position, ref allele, alt alleles, etc.
        let mut vcf_record = vcf::variant::RecordBuf::default();
        
        // Process each variant record in the VCF file
        // read_record_buf() returns number of bytes read (0 at EOF)
        while vcf_reader.read_record_buf(&header, &mut vcf_record)? != 0 {
            // Extract the genomic position where this variant occurs
            // variant_start() returns Option<Position> - the 1-indexed position of the variant
            let variant_pos = vcf_record.variant_start()
                .map(|p| usize::from(p))  // Convert Position to usize
                .unwrap_or(0);            // Default to 0 if position is missing
            
            // ===== VARIANT FILTERING =====
            // Check if this variant falls within the current BED region of interest
            // Uses inclusive bounds check: start_pos <= variant_pos <= end_pos
            if variant_pos >= start_pos && variant_pos <= end_pos {
                // Extract the reference allele (what the genome "normally" has at this position)
                let ref_bases = vcf_record.reference_bases().to_string();
                
                // Extract alternate alleles (the variant forms observed in samples)
                // This is a complex operation due to noodles API limitations
                let alt_bases = {
                    // Use debug formatting as a workaround to access alternate alleles
                    // Format will be something like ["A", "T"] for multiple alternates
                    let alt_str = format!("{:?}", vcf_record.alternate_bases());
                    
                    if alt_str.is_empty() || alt_str == "[]" {
                        // No alternate alleles found
                        ".".to_string()
                    } else {
                        // Parse the debug string to extract the first alternate allele
                        // Remove brackets, split by comma, take first element, clean quotes
                        alt_str.trim_start_matches('[').trim_end_matches(']')
                            .split(',').next().unwrap_or(".").trim().trim_matches('"').to_string()
                    }
                };
                
                // Create human-readable variant description in "REF>ALT" format
                // Example: "A>T" for an A-to-T substitution
                let variant_description = format!("{}>{}", ref_bases, alt_bases);

                // ===== BAM FILE PROCESSING =====
                // For each variant found, scan all sequencing reads to see which ones contain it
                // Note: This approach reads the entire BAM file for each variant - very inefficient!
                // A production tool would use BAM indexing for random access by genomic region
                
                let bam_file = File::open(&args.bam)?;
                // Create BAM reader with buffered I/O
                let mut bam_reader = bam::io::Reader::new(BufReader::new(bam_file));
                
                // Read BAM header containing reference sequence information and metadata
                // We don't use the header data but it's required for proper file parsing
                let _bam_header = bam_reader.read_header()?;
                
                // Create reusable BAM record buffer for memory efficiency
                // BAM records contain: read name, alignment position, CIGAR string, sequence, etc.
                let mut bam_record = bam::Record::default();
                
                // Process each sequencing read (alignment record) in the BAM file
                while bam_reader.read_record(&mut bam_record)? != 0 {
                    // ===== READ IDENTIFICATION =====
                    // Extract the read identifier/name from the BAM record
                    // Read names are typically assigned by the sequencing instrument
                    let read_id = bam_record.name()
                        .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))  // Convert bytes to string
                        .unwrap_or("unknown");  // Handle missing names gracefully
                    
                    // ===== VARIANT OVERLAP DETECTION =====
                    // Determine if this sequencing read overlaps with the variant position
                    // This is simplified logic - a full implementation would parse CIGAR strings
                    let contains_variant = if let Some(Ok(alignment_start)) = bam_record.alignment_start() {
                        // Get the genomic position where this read's alignment begins
                        let start_pos = usize::from(alignment_start);
                        
                        // SIMPLIFIED: Use fixed read length instead of parsing CIGAR
                        // Real implementation should calculate from CIGAR string operations
                        // CIGAR describes how the read aligns (matches, insertions, deletions, etc.)
                        let alignment_len = 100; // Simplified fixed length assumption
                        
                        // Calculate the end position of the read's alignment
                        let end_pos = start_pos + alignment_len;
                        
                        // Check if variant position falls within the read's alignment span
                        // This is a simple overlap test: read_start <= variant_pos <= read_end
                        variant_pos >= start_pos && variant_pos <= end_pos
                    } else {
                        // If alignment start is missing or invalid, assume no overlap
                        false
                    };

                    // ===== TIMESTAMP EXTRACTION =====
                    // Extract timing information from the read's auxiliary tags
                    // Nanopore reads often include start time in 'st:Z:<timestamp>' format
                    // 
                    // BAM auxiliary tags format: <TAG>:<TYPE>:<VALUE>
                    // - st: two-character tag name for "start time"
                    // - Z: type code for null-terminated string
                    // - <timestamp>: the actual timestamp value (format varies by basecaller)
                    //
                    // Common timestamp formats encountered:
                    // - ISO 8601: st:Z:2021-02-15T14:30:05.123Z (most human-readable)
                    // - Unix timestamp: st:Z:1613397005.123 (seconds since Unix epoch)
                    // - MinKNOW format: st:Z:2021-02-15 14:30:05.123+00:00
                    // - Numeric variants: st:i:1613397005 or st:f:1613397005.123
                    //
                    // This implementation preserves the original format for downstream analysis
                    let timestamp = {
                        // Create the 'st' tag identifier for start time
                        // BAM tags are identified by 2-byte arrays: [b's', b't'] = "st"
                        let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
                        
                        // Try to get the 'st' tag from the BAM record's auxiliary data
                        // The data() method returns a collection of auxiliary fields
                        match bam_record.data().get(&st_tag) {
                            Some(Ok(field_value)) => {
                                // Successfully found and parsed the 'st' tag
                                match field_value {
                                    sam::alignment::record::data::field::Value::String(ts_bytes) => {
                                        // Tag contains string data (Z type) - this is what we expect
                                        // Convert raw bytes to UTF-8 string representation
                                        std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                                    }
                                    sam::alignment::record::data::field::Value::Int8(val) => {
                                        // 8-bit signed integer timestamp (rare)
                                        format!("INT8_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt8(val) => {
                                        // 8-bit unsigned integer timestamp (rare)  
                                        format!("UINT8_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Int16(val) => {
                                        // 16-bit signed integer timestamp (rare)
                                        format!("INT16_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt16(val) => {
                                        // 16-bit unsigned integer timestamp (rare)
                                        format!("UINT16_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Int32(val) => {
                                        // 32-bit signed integer - might be Unix timestamp
                                        format!("INT32_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::UInt32(val) => {
                                        // 32-bit unsigned integer - might be Unix timestamp
                                        format!("UINT32_{}", val)
                                    }
                                    sam::alignment::record::data::field::Value::Float(val) => {
                                        // Float timestamp - might be Unix timestamp with fractional seconds
                                        format!("FLOAT_{}", val)
                                    }
                                    _ => {
                                        // Tag exists but contains unexpected data type
                                        // Could be array, hex string, etc.
                                        format!("UNKNOWN_TYPE_{:?}", field_value)
                                    }
                                }
                            }
                            Some(Err(e)) => {
                                // Tag was found but there was an error parsing its value
                                // This could indicate corrupted auxiliary data
                                format!("PARSE_ERROR_{:?}", e)
                            }
                            None => {
                                // 'st' tag not present in this read's auxiliary data
                                // This is normal for reads that don't have timing information
                                "NA".to_string()
                            }
                        }
                    };

                    // ===== OUTPUT GENERATION =====
                    // Write a tab-separated line to the output file for this read
                    // Format: read_id, timestamp, contains_variant, variant_description, variant_type, region, region_name
                    writeln!(
                        output,
                        "{}\t{}\t{}\t{}\tSNV\t{}\t{}",
                        read_id,              // Unique identifier for this sequencing read
                        timestamp,            // When the read was sequenced (simplified to "NA")
                        contains_variant,     // Boolean: does this read span the variant?
                        variant_description,  // Human-readable variant (e.g., "A>T")
                        // "SNV",             // Variant type (hardcoded as Single Nucleotide Variant)
                        region_string,        // Genomic region (e.g., "chr1:1000-2000")
                        region_name          // BED region name (always "NA" for BED3)
                    )?;
                }
            }
        }
    }

    // ===== SUCCESSFUL COMPLETION =====
    // If we reach this point, all file processing completed successfully
    // The output file now contains a comprehensive report of:
    // - Which sequencing reads overlap with variants in the specified genomic regions
    // - Timing information (when available)
    // - Variant descriptions in human-readable format
    Ok(())
}
