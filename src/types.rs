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

    /// Generate a unique variant group identifier for multi-allelic variants
    /// Format: chrom:pos:ref>alt1,alt2,alt3
    pub fn variant_group_id(&self) -> String {
        format!(
            "{}:{}:{}>{}",
            self.chrom,
            self.pos,
            self.ref_allele,
            self.alt_alleles.join(",")
        )
    }

    /// Generate a short variant group identifier (hash-based)
    /// Useful for compact identifiers while preserving uniqueness
    pub fn variant_group_hash(&self) -> String {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut hasher = DefaultHasher::new();
        self.chrom.hash(&mut hasher);
        self.pos.hash(&mut hasher);
        self.ref_allele.hash(&mut hasher);
        for alt in &self.alt_alleles {
            alt.hash(&mut hasher);
        }

        format!("VAR_{:08x}", hasher.finish() as u32)
    }

    /// Get a human-readable variant summary
    pub fn variant_summary(&self) -> String {
        if self.is_multiallelic() {
            format!(
                "{}:{} {}-allelic {}>{}",
                self.chrom,
                self.pos,
                self.num_alts(),
                self.ref_allele,
                self.alt_alleles.join(",")
            )
        } else {
            format!(
                "{}:{} {}>{}",
                self.chrom,
                self.pos,
                self.ref_allele,
                self.alt_alleles.join(",")
            )
        }
    }
}
