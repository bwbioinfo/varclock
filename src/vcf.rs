//! VCF file processing and variant extraction
//!
//! This module provides functionality for:
//! - Querying indexed VCF files for specific genomic regions
//! - Extracting and classifying variants (SNV, INS, DEL, BND, etc.)
//! - Parsing complex variant types including structural variants

use noodles_core;
use noodles_vcf as vcf;
use std::path::PathBuf;

// Import the AlternateBases trait to use iter() method
use noodles_vcf::variant::record::AlternateBases;

use crate::types::Variant;

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
pub fn query_vcf_variants_for_region(
    vcf_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<Vec<Variant>, Box<dyn std::error::Error + Send + Sync>> {
    // Use indexed access - no fallback to ensure consistent performance
    // This guarantees O(log n + k) complexity where k is the number of overlapping variants
    match query_vcf_variants_indexed(vcf_path, chrom, start_pos, end_pos) {
        Ok(variants) => Ok(variants),
        Err(index_error) => {
            let error_msg = format!(
                "VCF index access failed for {}. Please ensure the VCF file is indexed with 'tabix -p vcf'. Error: {:?}",
                vcf_path.display(),
                index_error
            );
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
    // STEP 1: Diagnostic logging for path debugging
    eprintln!("DEBUG: VCF file path: {}", vcf_path.display());
    eprintln!(
        "DEBUG: Current working directory: {:?}",
        std::env::current_dir()
    );
    eprintln!("DEBUG: VCF file exists: {}", vcf_path.exists());
    eprintln!("DEBUG: VCF file is absolute: {}", vcf_path.is_absolute());

    // Check for index file existence
    let index_path = vcf_path.with_extension("vcf.gz.tbi");
    let alt_index_path = format!("{}.tbi", vcf_path.display());
    eprintln!("DEBUG: Expected index path 1: {}", index_path.display());
    eprintln!("DEBUG: Index path 1 exists: {}", index_path.exists());
    eprintln!("DEBUG: Expected index path 2: {alt_index_path}");
    eprintln!(
        "DEBUG: Index path 2 exists: {}",
        std::path::Path::new(&alt_index_path).exists()
    );

    // Check if file is gzipped and indexed
    let is_gzipped = vcf_path
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext == "gz")
        .unwrap_or(false);

    if !is_gzipped {
        return Err("VCF file must be bgzip-compressed and tabix-indexed for efficient region queries. Please run: bgzip file.vcf && tabix -p vcf file.vcf.gz".into());
    }

    // STEP 2: Initialize indexed VCF reader
    // This automatically looks for and loads the corresponding .tbi index file
    eprintln!("DEBUG: Attempting to build indexed reader...");
    let mut indexed_reader =
        match vcf::io::indexed_reader::Builder::default().build_from_path(vcf_path) {
            Ok(reader) => {
                eprintln!("DEBUG: Successfully created indexed reader");
                reader
            }
            Err(e) => {
                eprintln!("DEBUG: Failed to create indexed reader: {e:?}");
                return Err(format!(
                    "Failed to create indexed VCF reader for {}: {:?}",
                    vcf_path.display(),
                    e
                )
                .into());
            }
        };

    // Read the VCF header containing meta-information and sample data
    eprintln!("DEBUG: Attempting to read VCF header...");
    let header = match indexed_reader.read_header() {
        Ok(h) => {
            eprintln!("DEBUG: Successfully read VCF header");

            // Debug: List available reference sequences
            eprintln!("DEBUG: Available reference sequences in VCF:");
            for (i, (name, _)) in h.contigs().iter().enumerate() {
                eprintln!("  {i}: {name}");
                if i >= 10 {
                    eprintln!("  ... and {} more", h.contigs().len() - 10);
                    break;
                }
            }
            eprintln!("DEBUG: Looking for chromosome: {chrom}");

            h
        }
        Err(e) => {
            eprintln!("DEBUG: Failed to read VCF header: {e:?}");
            return Err(format!("Failed to read VCF header: {e:?}").into());
        }
    };

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
    eprintln!(
        "DEBUG: Attempting to query region: {chrom}:{start_pos}-{end_pos}"
    );
    let query = match indexed_reader.query(&header, &region) {
        Ok(q) => {
            eprintln!("DEBUG: Successfully created region query");
            q
        }
        Err(e) => {
            eprintln!("DEBUG: Failed to create region query: {e:?}");
            return Err(format!(
                "Failed to query region {chrom}:{start_pos}-{end_pos}: {e:?}"
            )
            .into());
        }
    };

    // Container for extracted variant information
    let mut variants = Vec::new();

    // STEP 5: Process each variant returned by the indexed query
    for result in query {
        // Parse the variant record (handle potential I/O errors)
        let record = result?;

        // SUBSTEP 5A: Extract genomic position of the variant
        let variant_pos = record
            .variant_start()
            .map(|p| match p {
                Ok(pos) => usize::from(pos),
                Err(_) => 0,
            })
            .unwrap_or(0);

        // SUBSTEP 5B: Extract reference and ALL alternative alleles
        let ref_bases = record.reference_bases().to_string();

        // Collect all alternate alleles
        let mut alt_alleles = Vec::new();
        let mut descriptions = Vec::new();
        let mut variant_types = Vec::new();

        for alt_allele_result in record.alternate_bases().iter() {
            let alt_bases = match alt_allele_result {
                Ok(alt_allele) => alt_allele.to_string(),
                Err(_) => ".".to_string(),
            };

            // Skip missing alleles
            if alt_bases == "." || alt_bases.is_empty() {
                continue;
            }

            // SUBSTEP 5C: Classify variant type and create human-readable description for each allele
            let (description, variant_type) = if alt_bases.contains('[') || alt_bases.contains(']')
            {
                // Breakend (BND) variants - structural variations
                let bnd_description = parse_bnd_variant(&ref_bases, &alt_bases);
                (bnd_description, "BND".to_string())
            } else if ref_bases.len() == 1 && alt_bases.len() == 1 {
                // Single nucleotide variant (SNV)
                (format!("{ref_bases}>{alt_bases}"), "SNV".to_string())
            } else if ref_bases.len() > alt_bases.len() {
                // Deletion (DEL)
                (format!("{ref_bases}>{alt_bases}"), "DEL".to_string())
            } else if ref_bases.len() < alt_bases.len() {
                // Insertion (INS)
                (format!("{ref_bases}>{alt_bases}"), "INS".to_string())
            } else if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 {
                // Multi-nucleotide variant (MNV)
                (format!("{ref_bases}>{alt_bases}"), "MNV".to_string())
            } else {
                // Other/complex variant types
                (format!("{ref_bases}>{alt_bases}"), "OTHER".to_string())
            };

            alt_alleles.push(alt_bases);
            descriptions.push(description);
            variant_types.push(variant_type);
        }

        // Skip variants with no valid alternate alleles
        if alt_alleles.is_empty() {
            continue;
        }

        // STEP 6: Extract depth and quality information
        let quality = record.quality_score().and_then(|q| q.ok());

        // For multi-allelic variants, we'll need to extract depths for each allele
        let total_depth = None; // TODO: Extract DP from INFO
        let ref_depth = None; // TODO: Extract RO or AD[0]

        // Create placeholder depths for each alternate allele
        let alt_depths = vec![None; alt_alleles.len()]; // TODO: Extract AO or AD[1..] for each alt

        // STEP 7: Store extracted variant information with support for multiple alleles
        let variant = Variant {
            chrom: chrom.to_string(),
            pos: variant_pos,
            ref_allele: ref_bases.clone(),
            alt_alleles,
            descriptions,
            variant_types,
            total_depth,
            ref_depth,
            alt_depths,
            quality,
        };

        variants.push(variant);
    }

    Ok(variants)
}

/// Parse breakend (BND) variant notation
///
/// Breakend variants represent structural variations and use bracket notation
/// to indicate the mate position/orientation.
pub fn parse_bnd_variant(ref_bases: &str, alt_bases: &str) -> String {
    // Extract the breakend information from the alt allele
    if let Some(bracket_start) = alt_bases.find('[') {
        if let Some(bracket_end) = alt_bases.find(']') {
            // Format: s1[p[s2 or s1]p]s2
            let mate_info = &alt_bases[bracket_start + 1..bracket_end];
            return format!("BND_{ref_bases}>{alt_bases}_to_{mate_info}");
        }
    } else if let Some(bracket_start) = alt_bases.find(']')
        && let Some(bracket_end) = alt_bases.rfind('[')
    {
        // Format: ]p]s2 or s1]p[
        let mate_info = &alt_bases[bracket_start + 1..bracket_end];
        return format!("BND_{ref_bases}>{alt_bases}_to_{mate_info}");
    }

    // Fallback for complex BND notation
    format!("BND_{ref_bases}>{alt_bases}")
}

/// Debug function to test VCF file access and identify path issues
///
/// This function performs basic checks on VCF file accessibility and indexing
/// to help diagnose issues with absolute paths from different directories.
pub fn debug_vcf_access(
    vcf_path: &PathBuf,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    println!("=== VCF File Access Debug ===");
    println!("VCF path: {}", vcf_path.display());
    println!("Current working directory: {:?}", std::env::current_dir()?);
    println!("VCF file exists: {}", vcf_path.exists());
    println!("VCF file is absolute: {}", vcf_path.is_absolute());

    if let Some(parent) = vcf_path.parent() {
        println!("Parent directory: {}", parent.display());
        println!("Parent directory exists: {}", parent.exists());
    }

    // Check file permissions by trying to open it
    match std::fs::File::open(vcf_path) {
        Ok(_) => println!("✓ VCF file is readable"),
        Err(e) => println!("✗ Cannot read VCF file: {e}"),
    }

    // Check for index files in multiple possible locations
    let possible_index_paths = [format!("{}.tbi", vcf_path.display()),
        vcf_path.with_extension("vcf.gz.tbi").display().to_string(),
        vcf_path.with_extension("tbi").display().to_string()];

    println!("\nChecking for index files:");
    for (i, index_path) in possible_index_paths.iter().enumerate() {
        let path = std::path::Path::new(index_path);
        println!("  Index option {}: {}", i + 1, index_path);
        println!("    Exists: {}", path.exists());
        if path.exists() {
            match std::fs::File::open(path) {
                Ok(_) => println!("    ✓ Readable"),
                Err(e) => println!("    ✗ Not readable: {e}"),
            }
        }
    }

    // Try to create the indexed reader
    println!("\nTesting indexed reader creation:");
    match vcf::io::indexed_reader::Builder::default().build_from_path(vcf_path) {
        Ok(mut reader) => {
            println!("✓ Successfully created indexed reader");

            // Test header reading and list available chromosomes
            match reader.read_header() {
                Ok(header) => {
                    println!("✓ Successfully read VCF header");
                    println!("\nAvailable reference sequences:");
                    for (i, (name, _)) in header.contigs().iter().enumerate() {
                        println!("  {i}: {name}");
                        if i >= 20 {
                            println!("  ... and {} more", header.contigs().len() - 20);
                            break;
                        }
                    }
                }
                Err(e) => println!("✗ Failed to read VCF header: {e:?}"),
            }
        }
        Err(e) => {
            println!("✗ Failed to create indexed reader: {e:?}");
            return Err(e.into());
        }
    }

    println!("=== Debug Complete ===\n");
    Ok(())
}

/// Compare VCF access from current directory vs absolute path
///
/// This function helps diagnose differences between accessing VCF files
/// from the same directory vs using absolute paths from different directories.
pub fn compare_vcf_access_methods(
    vcf_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    println!("=== VCF Access Method Comparison ===");

    let original_dir = std::env::current_dir()?;
    println!("Original directory: {}", original_dir.display());

    // Method 1: Absolute path from current directory
    println!("\n--- Method 1: Absolute Path from Current Directory ---");
    let abs_result = test_vcf_query(vcf_path, chrom, start_pos, end_pos);
    match abs_result {
        Ok(count) => println!("✓ Found {count} variants using absolute path"),
        Err(e) => println!("✗ Failed with absolute path: {e}"),
    }

    // Method 2: From VCF's directory with relative path
    if let Some(vcf_dir) = vcf_path.parent()
        && let Some(filename) = vcf_path.file_name() {
            println!("\n--- Method 2: Relative Path from VCF Directory ---");
            println!("Changing to directory: {}", vcf_dir.display());

            match std::env::set_current_dir(vcf_dir) {
                Ok(_) => {
                    let relative_path = PathBuf::from(filename);
                    let rel_result = test_vcf_query(&relative_path, chrom, start_pos, end_pos);
                    match rel_result {
                        Ok(count) => println!("✓ Found {count} variants using relative path"),
                        Err(e) => println!("✗ Failed with relative path: {e}"),
                    }

                    // Restore original directory
                    if let Err(e) = std::env::set_current_dir(&original_dir) {
                        println!("WARNING: Failed to restore directory: {e}");
                    }
                }
                Err(e) => println!("✗ Failed to change to VCF directory: {e}"),
            }
        }

    println!("\n=== Comparison Complete ===");
    Ok(())
}

/// Test VCF query without detailed logging
fn test_vcf_query(
    vcf_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
) -> Result<usize, Box<dyn std::error::Error + Send + Sync>> {
    let mut indexed_reader =
        vcf::io::indexed_reader::Builder::default().build_from_path(vcf_path)?;

    let header = indexed_reader.read_header()?;

    let start_position = noodles_core::Position::try_from(start_pos)?;
    let end_position = noodles_core::Position::try_from(end_pos)?;
    let interval = start_position..=end_position;
    let region = noodles_core::Region::new(chrom, interval);

    let query = indexed_reader.query(&header, &region)?;

    let mut count = 0;
    for result in query {
        let _record = result?;
        count += 1;
    }

    Ok(count)
}
