// VarClock Library - Modular bioinformatics tool for variant analysis
//
// This library provides functionality for analyzing genomic variants by:
// 1. Reading genomic regions from BED files
// 2. Extracting variants from indexed VCF files
// 3. Analyzing sequencing reads from indexed BAM files
// 4. Performing sequence-based variant detection
// 5. Outputting compressed, indexed results

pub mod bam; // BAM file processing and read analysis
pub mod bed; // BED file parsing and region management
pub mod output; // Output formatting and compression
pub mod types;
pub mod variant; // Variant analysis and sequence comparison
pub mod vcf; // VCF file processing and variant extraction // Common data structures and types

// Re-export commonly used types for convenience
pub use bam::query_bam_records_for_region;
pub use bed::load_bed_regions;
pub use output::BgzOutput;
pub use types::{BedRegion, Variant};
pub use variant::{
    AlleleMatch, MultiAlleleMatch, analyze_read_allele_content_detailed,
    analyze_read_multiallelic_content, analyze_read_variant_content, classify_variant_type,
    extract_read_bases_at_position, parse_breakend_mate_position, read_spans_variant_position,
};
pub use vcf::{parse_bnd_variant, query_vcf_variants_for_region};
