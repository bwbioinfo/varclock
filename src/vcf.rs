//! VCF file processing and variant extraction  
//!
//! This module provides functionality for:
//! - Querying indexed VCF files for specific genomic regions
//! - Extracting and classifying variants (SNV, INS, DEL, BND, etc.)
//! - Parsing complex variant types including structural variants

use std::path::PathBuf;
use noodles_vcf as vcf;
use noodles_core;

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
        Ok(variants) => {
            Ok(variants)
        }
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
        
        // STEP 6: Extract depth and quality information
        let quality = record.quality_score().and_then(|q| q.ok());
        
        // For now, let's add placeholder depth information and implement a simpler approach
        // We'll extract basic INFO field data that we can access reliably
        let total_depth = None;  // TODO: Extract DP from INFO
        let ref_depth = None;    // TODO: Extract RO or AD[0] 
        let alt_depth = None;    // TODO: Extract AO or AD[1]
        
        // STEP 7: Store extracted variant information with depth placeholders
        let variant = Variant {
            chrom: chrom.to_string(),
            pos: variant_pos,
            ref_allele: ref_bases.clone(),
            alt_allele: alt_bases.clone(),
            description,
            variant_type,
            total_depth,
            ref_depth,
            alt_depth,
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
        && let Some(bracket_end) = alt_bases.rfind('[') {
            // Format: ]p]s2 or s1]p[
            let mate_info = &alt_bases[bracket_start + 1..bracket_end];
            return format!("BND_{ref_bases}>{alt_bases}_to_{mate_info}");
        }
    
    // Fallback for complex BND notation
    format!("BND_{ref_bases}>{alt_bases}")
}
