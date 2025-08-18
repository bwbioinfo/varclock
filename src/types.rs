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
/// Now supports multiple alternate alleles (multi-allelic variants)
#[derive(Debug, Clone)]
pub struct Variant {
    pub chrom: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,     // Multiple alternate alleles
    pub descriptions: Vec<String>,    // Description for each alternate allele
    pub variant_types: Vec<String>,   // Type for each alternate allele (SNV, INS, DEL, etc.)
    pub total_depth: Option<u32>,     // DP from INFO or FORMAT
    pub ref_depth: Option<u32>,       // Reference allele depth (RO or from AD)
    pub alt_depths: Vec<Option<u32>>, // Alternative allele depths for each alt allele
    pub quality: Option<f32>,         // QUAL field
}

impl Variant {
    /// Get the number of alternate alleles
    pub fn num_alts(&self) -> usize {
        self.alt_alleles.len()
    }

    /// Check if this is a multi-allelic variant
    pub fn is_multiallelic(&self) -> bool {
        self.alt_alleles.len() > 1
    }

    /// Get allele info for a specific alternate allele index
    pub fn get_allele_info(&self, index: usize) -> Option<(&str, &str, &str, Option<u32>)> {
        if index < self.alt_alleles.len() {
            Some((
                &self.alt_alleles[index],
                &self.descriptions[index],
                &self.variant_types[index],
                self.alt_depths.get(index).copied().flatten(),
            ))
        } else {
            None
        }
    }
}
