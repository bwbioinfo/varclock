// Common data structures and types used throughout VarClock

/// Represents a genomic region from a BED file
#[derive(Debug, Clone)]
pub struct BedRegion {
    pub chrom: String,
    pub start_pos: usize,
    pub end_pos: usize,
    pub region_name: String,
    pub region_string: String,
}

/// Represents a genetic variant from a VCF file
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub description: String,
    pub variant_type: String,
    pub total_depth: Option<u32>,      // DP from INFO or FORMAT
    pub ref_depth: Option<u32>,        // Reference allele depth (RO or from AD)
    pub alt_depth: Option<u32>,        // Alternative allele depth (AO or from AD)
    pub quality: Option<f32>,          // QUAL field
}
