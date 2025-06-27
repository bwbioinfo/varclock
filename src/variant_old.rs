//! Variant analysis and sequence comparison
//!
//! This module provides functionality for:
//! - Checking if reads span variant positions
//! - Performing sequence-based variant detection using CIGAR and read sequences
//! - Extracting bases from BAM records at specific positions

use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_bam as bam;

/// Check if a read spans a specific variant position
/// 
/// This function only checks if the read positionally spans the variant position.
/// This is used as a first filter before sequence-based analysis.
/// 
/// **Parameters:**
/// - `read_start`: Start position of the read (1-based genomic coordinates)
/// - `read_end`: End position of the read (1-based genomic coordinates)
/// - `variant_pos`: Position of the variant (1-based genomic coordinates)
/// 
/// **Returns:**
/// - `true`: Read spans the variant position
/// - `false`: Read does not span the variant position
pub fn read_spans_variant_position(read_start: usize, read_end: usize, variant_pos: usize) -> bool {
    variant_pos >= read_start && variant_pos <= read_end
}

/// Analyze read sequence to determine if it contains the variant allele
/// 
/// This function performs true sequence-based variant detection by:
/// 1. Checking if the read spans the variant position
/// 2. Extracting the read sequence from the BAM record
/// 3. Using CIGAR operations to map reference positions to read positions
/// 4. Comparing the read bases at the variant position to the expected alleles
/// 
/// **Parameters:**
/// - `record`: The BAM record containing read data
/// - `variant_pos`: 1-based genomic position of the variant
/// - `ref_allele`: Reference allele sequence
/// - `alt_allele`: Alternative allele sequence
/// 
/// **Returns:**
/// - `Some(true)`: Read contains the alternative allele
/// - `Some(false)`: Read contains the reference allele  
/// - `None`: Cannot determine (read doesn't span variant, indeterminate sequence, etc.)
pub fn analyze_read_variant_content(
    record: &bam::Record, 
    variant_pos: usize, 
    ref_allele: &str, 
    alt_allele: &str
) -> Option<bool> {
    // STEP 1: Check if read spans the variant position
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return None,
    };
    
    // Calculate read end using CIGAR operations
    let cigar = record.cigar();
    let mut reference_pos = read_start;
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let op_len = usize::from(op.len());
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                    reference_pos += op_len;
                }
                _ => {} // Insertions, clipping, etc. don't advance reference position
            }
        }
    }
    let read_end = if reference_pos > read_start { reference_pos - 1 } else { read_start };
    
    // Check if variant position is within read span
    if variant_pos < read_start || variant_pos > read_end {
        return None; // Read doesn't span the variant
    }
    
    // STEP 2: Extract read sequence
    let read_sequence = record.sequence();
    
    // STEP 3: Map reference position to read position using CIGAR
    let mut ref_pos = read_start;
    let mut read_pos = 0usize;
    
    // DEBUG: Print CIGAR operations for INDEL debugging
    if ref_allele.len() != alt_allele.len() {
        print!("        DEBUG: INDEL variant {}:{} {}>{} - CIGAR operations: ", 
               record.alignment_start().map(|s| s.map(|p| usize::from(p))).unwrap_or(Ok(0)).unwrap_or(0),
               variant_pos, ref_allele, alt_allele);
        for operation in cigar.iter() {
            if let Ok(op) = operation {
                print!("{}({}) ", 
                       match op.kind() {
                           Kind::Match => "M",
                           Kind::SequenceMatch => "=", 
                           Kind::SequenceMismatch => "X",
                           Kind::Insertion => "I",
                           Kind::Deletion => "D",
                           Kind::Skip => "N",
                           Kind::SoftClip => "S",
                           Kind::HardClip => "H",
                           Kind::Pad => "P",
                       },
                       op.len());
            }
        }
        println!();
    }
    
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let op_len = usize::from(op.len());
            
            // Check if we've reached the variant position
            if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        // Calculate exact position within the operation
                        let offset_in_op = variant_pos - ref_pos;
                        let read_variant_pos = read_pos + offset_in_op;
                        
                        // STEP 4: Extract bases at variant position
                        let read_bases = extract_read_bases_at_position(&read_sequence, read_variant_pos, ref_allele.len());
                        
                        // STEP 5: Compare with reference and alternative alleles
                        if read_bases == alt_allele {
                            return Some(true); // Contains alternative allele
                        } else if read_bases == ref_allele {
                            return Some(false); // Contains reference allele
                        } else {
                            return None; // Contains neither (sequencing error, complex variant, etc.)
                        }
                    }
                    Kind::Deletion => {
                        // Handle deletion variants
                        let del_result = handle_deletion_variant(
                            record, ref_pos, read_pos, op_len, 
                            variant_pos, ref_allele, alt_allele
                        );
                        
                        match del_result {
                            AlleleMatch::Variant(_) => return Some(true),
                            AlleleMatch::Reference(_) => return Some(false), 
                            _ => return None,
                        }
                    }
                    Kind::Insertion => {
                        // Handle insertion variants  
                        let ins_result = handle_insertion_variant(
                            record, ref_pos, read_pos, op_len,
                            variant_pos, ref_allele, alt_allele
                        );
                        
                        match ins_result {
                            AlleleMatch::Variant(_) => return Some(true),
                            AlleleMatch::Reference(_) => return Some(false),
                            _ => return None,
                        }
                    }
                    _ => {
                        // Other operations (clipping, etc.)
                        return None;
                    }
                }
            }
            
            // Special case: Check for insertion variants at the boundary
            if ref_pos == variant_pos && op.kind() == Kind::Insertion {
                let ins_result = handle_insertion_variant(
                    record, ref_pos, read_pos, op_len,
                    variant_pos, ref_allele, alt_allele
                );
                
                match ins_result {
                    AlleleMatch::Variant(_) => return Some(true),
                    AlleleMatch::Reference(_) => return Some(false),
                    _ => return None,
                }
            }
            
            // Advance positions based on operation type
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += op_len;
                    // read_pos stays the same
                }
                Kind::Insertion | Kind::SoftClip => {
                    read_pos += op_len;
                    // ref_pos stays the same
                }
                Kind::HardClip | Kind::Pad => {
                    // Neither position advances
                }
            }
        }
    }
    
    None // Variant position not found in read
}

/// Extract bases from read sequence at a specific position
/// 
/// **Parameters:**
/// - `sequence`: The read sequence
/// - `start_pos`: Starting position in the read (0-based)
/// - `length`: Number of bases to extract
/// 
/// **Returns:**
/// - String containing the extracted bases
pub fn extract_read_bases_at_position(sequence: &bam::record::Sequence, start_pos: usize, length: usize) -> String {
    let mut bases = String::new();
    
    for i in 0..length {
        if let Some(base) = sequence.get(start_pos + i) {
            // The noodles library stores sequences as ASCII values
            bases.push(match base {
                // ASCII encoding (most common)
                b'A' => 'A',
                b'C' => 'C', 
                b'G' => 'G',
                b'T' => 'T',
                b'N' => 'N',
                // Standard BAM 4-bit encoding as fallback
                0 => '=',  // Reference skip
                1 => 'A',  // A = 1  
                2 => 'C',  // C = 2
                4 => 'G',  // G = 4
                8 => 'T',  // T = 8
                15 => 'N', // N = 15
                3 => 'G',  // Sometimes G = 3
                5 => 'T',  // Sometimes T = 5
                _ => {
                    // For any other value, use as ASCII if printable
                    if base >= 32 && base <= 126 {
                        base as char
                    } else {
                        'N'
                    }
                },
            });
        } else {
            bases.push('N'); // Position beyond read end
        }
    }
    
    bases
}

/// Result of analyzing read content for variant alleles
#[derive(Debug, Clone, PartialEq)]
pub enum AlleleMatch {
    /// Read contains the reference allele at the variant position
    Reference(String),
    /// Read contains the variant (alternative) allele at the variant position
    Variant(String),
    /// Read contains a different sequence (neither ref nor alt) - could be sequencing error or complex variant
    Other(String),
    /// Read doesn't span the variant position
    NoSpan,
    /// Variant position falls in a deletion in the read
    Deletion,
    /// Indeterminate - couldn't extract or analyze sequence
    Indeterminate,
}

impl AlleleMatch {
    /// Convert to string representation for output
    pub fn to_output_string(&self) -> String {
        match self {
            AlleleMatch::Reference(seq) => format!("REFERENCE:{}", seq),
            AlleleMatch::Variant(seq) => format!("VARIANT:{}", seq),
            AlleleMatch::Other(seq) => format!("OTHER:{}", seq),
            AlleleMatch::NoSpan => "NO_SPAN".to_string(),
            AlleleMatch::Deletion => "DELETION".to_string(),
            AlleleMatch::Indeterminate => "INDETERMINATE".to_string(),
        }
    }
    
    /// Check if this represents a definitive match (either reference or variant)
    pub fn is_definitive_match(&self) -> bool {
        matches!(self, AlleleMatch::Reference(_) | AlleleMatch::Variant(_))
    }
    
    /// Check if this represents a variant allele match
    pub fn is_variant_match(&self) -> bool {
        matches!(self, AlleleMatch::Variant(_))
    }
    
    /// Check if this represents a reference allele match  
    pub fn is_reference_match(&self) -> bool {
        matches!(self, AlleleMatch::Reference(_))
    }
}

/// Analyze read sequence to determine detailed allele content  
/// 
/// This is an enhanced version of analyze_read_variant_content that provides
/// more detailed information about what was found in the read sequence.
///
/// **Parameters:**
/// - `record`: The BAM record containing read data
/// - `variant_pos`: 1-based genomic position of the variant
/// - `ref_allele`: Reference allele sequence
/// - `alt_allele`: Alternative allele sequence
/// 
/// **Returns:**
/// - `AlleleMatch` enum with detailed information about the match
pub fn analyze_read_allele_content_detailed(
    record: &bam::Record, 
    variant_pos: usize, 
    ref_allele: &str, 
    alt_allele: &str
) -> AlleleMatch {
    // STEP 1: Check if read spans the variant position
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return AlleleMatch::Indeterminate,
    };
    
    // Calculate read end using CIGAR operations
    let cigar = record.cigar();
    let mut reference_pos = read_start;
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            let op_len = usize::from(op.len());
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                    reference_pos += op_len;
                }
                _ => {} // Insertions, clipping, etc. don't advance reference position
            }
        }
    }
    let read_end = if reference_pos > read_start { reference_pos - 1 } else { read_start };
    
    // Check if variant position is within read span
    if variant_pos < read_start || variant_pos > read_end {
        return AlleleMatch::NoSpan;
    }
    
    // STEP 2: Determine variant type and use appropriate analysis
    let variant_type = classify_variant_type(ref_allele, alt_allele);
    
    match variant_type {
        "SNV" | "MNV" => {
            // Use sequence-based analysis for simple substitutions
            analyze_snv_from_cigar(record, variant_pos, ref_allele, alt_allele)
        }
        "DEL" => {
            // Use CIGAR-based analysis for deletions
            analyze_deletion_from_cigar(record, variant_pos, ref_allele, alt_allele)
        }
        "INS" => {
            // Use CIGAR-based analysis for insertions
            analyze_insertion_from_cigar(record, variant_pos, ref_allele, alt_allele)
        }
        _ => {
            // Complex variants - try sequence-based analysis
            analyze_snv_from_cigar(record, variant_pos, ref_allele, alt_allele)
        }
    }
}
    // STEP 1: Check if read spans the variant position
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return AlleleMatch::Indeterminate,
    };
    
    // Calculate read end using CIGAR operations
    let cigar = record.cigar();
    let mut reference_pos = read_start;
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let op_len = usize::from(op.len());
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                    reference_pos += op_len;
                }
                _ => {} // Insertions, clipping, etc. don't advance reference position
            }
        }
    }
    let read_end = if reference_pos > read_start { reference_pos - 1 } else { read_start };
    
    // Check if variant position is within read span
    if variant_pos < read_start || variant_pos > read_end {
        return AlleleMatch::NoSpan;
    }
    
    // STEP 2: Extract read sequence
    let read_sequence = record.sequence();
    
    // STEP 3: Map reference position to read position using CIGAR
    let mut ref_pos = read_start;
    let mut read_pos = 0usize;
    
    // DEBUG: Print CIGAR operations for INDEL debugging
    if ref_allele.len() != alt_allele.len() {
        print!("        DEBUG: INDEL variant {}:{} {}>{} - CIGAR operations: ", 
               record.alignment_start().map(|s| s.map(|p| usize::from(p))).unwrap_or(Ok(0)).unwrap_or(0),
               variant_pos, ref_allele, alt_allele);
        for operation in cigar.iter() {
            if let Ok(op) = operation {
                print!("{}({}) ", 
                       match op.kind() {
                           Kind::Match => "M",
                           Kind::SequenceMatch => "=", 
                           Kind::SequenceMismatch => "X",
                           Kind::Insertion => "I",
                           Kind::Deletion => "D",
                           Kind::Skip => "N",
                           Kind::SoftClip => "S",
                           Kind::HardClip => "H",
                           Kind::Pad => "P",
                       },
                       op.len());
            }
        }
        println!();
    }
    
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let op_len = usize::from(op.len());
            
            // Check if we've reached the variant position
            if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        // For simple substitutions and reference matches
                        if ref_allele.len() == alt_allele.len() {
                            // Calculate exact position within the operation
                            let offset_in_op = variant_pos - ref_pos;
                            let read_variant_pos = read_pos + offset_in_op;
                            
                            // STEP 4: Extract bases at variant position
                            let read_bases = extract_read_bases_at_position(&read_sequence, read_variant_pos, ref_allele.len());
                            
                            // STEP 5: Compare with reference and alternative alleles
                            if read_bases == alt_allele {
                                return AlleleMatch::Variant(read_bases);
                            } else if read_bases == ref_allele {
                                return AlleleMatch::Reference(read_bases);
                            } else {
                                return AlleleMatch::Other(read_bases);
                            }
                        } else {
                            // For INDELs, we need to check if this position contains the reference allele
                            // in a region without indels (supporting reference)
                            let offset_in_op = variant_pos - ref_pos;
                            let read_variant_pos = read_pos + offset_in_op;
                            
                            // Extract the first base to see if this is a reference match
                            let read_base = extract_read_bases_at_position(&read_sequence, read_variant_pos, 1);
                            
                            // For deletion variants (REF longer), if we see the first base of REF, it's reference 
                            // For insertion variants (ALT longer), if we see just REF, it's reference
                            if ref_allele.len() > alt_allele.len() {
                                // Deletion variant - check if we see the full reference sequence
                                let read_bases = extract_read_bases_at_position(&read_sequence, read_variant_pos, ref_allele.len());
                                if read_bases == ref_allele {
                                    return AlleleMatch::Reference(read_bases);
                                }
                            } else {
                                // Insertion variant - if we only see the reference base, it's reference
                                if read_base == ref_allele {
                                    return AlleleMatch::Reference(read_base);
                                }
                            }
                            
                            return AlleleMatch::Other(format!("MATCH_OP:{}", read_base));
                        }
                    }
                    Kind::Deletion => {
                        // Handle deletion variants - check if this CIGAR deletion matches the variant
                        return handle_deletion_variant(
                            record, ref_pos, read_pos, op_len, 
                            variant_pos, ref_allele, alt_allele
                        );
                    }
                    Kind::Insertion => {
                        // Handle insertion variants - check if this CIGAR insertion matches the variant
                        return handle_insertion_variant(
                            record, ref_pos, read_pos, op_len,
                            variant_pos, ref_allele, alt_allele
                        );
                    }
                    _ => {
                        // Other operations (clipping, etc.)
                        return AlleleMatch::Indeterminate;
                    }
                }
            }
            
            // Special case: Check for insertion variants at the boundary between operations
            if ref_pos == variant_pos && op.kind() == Kind::Insertion {
                return handle_insertion_variant(
                    record, ref_pos, read_pos, op_len,
                    variant_pos, ref_allele, alt_allele
                );
            }
            
            // Advance positions based on operation type
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += op_len;
                    // read_pos stays the same
                }
                Kind::Insertion | Kind::SoftClip => {
                    read_pos += op_len;
                    // ref_pos stays the same
                }
                Kind::HardClip | Kind::Pad => {
                    // Neither position advances
                }
            }
        }
    }
    
    AlleleMatch::Indeterminate // Variant position not found in read
}

/// Handle deletion variants by comparing CIGAR deletion to expected variant
/// 
/// **Parameters:**
/// - `record`: The BAM record
/// - `cigar_ref_pos`: Reference position where the CIGAR deletion starts
/// - `cigar_read_pos`: Read position where the CIGAR deletion occurs
/// - `deletion_length`: Length of the CIGAR deletion
/// - `variant_pos`: Position of the variant
/// - `ref_allele`: Reference allele (should be longer for deletions)
/// - `alt_allele`: Alternative allele (should be shorter for deletions)
fn handle_deletion_variant(
    record: &bam::Record,
    cigar_ref_pos: usize,
    cigar_read_pos: usize,
    deletion_length: usize,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> AlleleMatch {
    // Check if this is actually a deletion variant (ref longer than alt)
    if ref_allele.len() <= alt_allele.len() {
        // Not a deletion variant, this might be a different type
        return AlleleMatch::Indeterminate;
    }
    
    let expected_deletion_length = ref_allele.len() - alt_allele.len();
    
    // Check if the CIGAR deletion overlaps with the variant position
    let deletion_start = cigar_ref_pos;
    let deletion_end = cigar_ref_pos + deletion_length - 1;
    
    if variant_pos >= deletion_start && variant_pos <= deletion_end {
        // This deletion overlaps our variant position
        
        // Simple check: if the deletion length matches expected, consider it a variant match
        if deletion_length >= expected_deletion_length {
            return AlleleMatch::Variant(format!("DEL:{}>{}", ref_allele, alt_allele));
        } else {
            return AlleleMatch::Other(format!("DEL:{}bp@{}", deletion_length, cigar_ref_pos));
        }
    }
    
    // This deletion doesn't affect our variant position
    return AlleleMatch::Deletion;
}

/// Handle insertion variants by comparing CIGAR insertion to expected variant
/// 
/// **Parameters:**
/// - `record`: The BAM record
/// - `cigar_ref_pos`: Reference position where the CIGAR insertion occurs
/// - `cigar_read_pos`: Read position where the CIGAR insertion starts
/// - `insertion_length`: Length of the CIGAR insertion
/// - `variant_pos`: Position of the variant
/// - `ref_allele`: Reference allele (should be shorter for insertions)
/// - `alt_allele`: Alternative allele (should be longer for insertions)
fn handle_insertion_variant(
    record: &bam::Record,
    cigar_ref_pos: usize,
    cigar_read_pos: usize,
    insertion_length: usize,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
) -> AlleleMatch {
    // Check if this is actually an insertion variant (alt longer than ref)
    if alt_allele.len() <= ref_allele.len() {
        // Not an insertion variant
        return AlleleMatch::Indeterminate;
    }
    
    let expected_insertion_length = alt_allele.len() - ref_allele.len();
    
    // Check if the CIGAR insertion is at or near the variant position
    if cigar_ref_pos == variant_pos || cigar_ref_pos == variant_pos + 1 {
        let sequence = record.sequence();
        
        // Extract the inserted sequence from the read
        let inserted_bases = extract_read_bases_at_position(&sequence, cigar_read_pos, insertion_length);
        
        // For insertions, we need to be more flexible about matching
        if insertion_length >= expected_insertion_length {
            return AlleleMatch::Variant(format!("INS:{}>{}", ref_allele, format!("{}{}", ref_allele, inserted_bases)));
        } else {
            return AlleleMatch::Other(format!("INS:{}bp@{}", insertion_length, cigar_ref_pos));
        }
    }
    
    // CIGAR insertion doesn't match the expected variant location
    return AlleleMatch::Indeterminate;
}

/// Determine the type of variant based on allele lengths
/// 
/// **Returns:**
/// - `"SNV"`: Single nucleotide variant
/// - `"INS"`: Insertion
/// - `"DEL"`: Deletion  
/// - `"COMPLEX"`: Complex variant
pub fn classify_variant_type(ref_allele: &str, alt_allele: &str) -> &'static str {
    match (ref_allele.len(), alt_allele.len()) {
        (r, a) if r == a && r == 1 => "SNV",
        (r, a) if r == a && r > 1 => "MNV", // Multi-nucleotide variant
        (1, a) if a > 1 => "INS", // Simple insertion: single base -> multiple bases
        (r, 1) if r > 1 => "DEL", // Simple deletion: multiple bases -> single base
        _ => "COMPLEX" // Complex variants: everything else
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_variant_classification() {
        // Test SNV
        assert_eq!(classify_variant_type("A", "G"), "SNV");
        assert_eq!(classify_variant_type("T", "C"), "SNV");
        
        // Test MNV (multi-nucleotide variant)
        assert_eq!(classify_variant_type("AT", "GC"), "MNV");
        assert_eq!(classify_variant_type("ATG", "GCA"), "MNV");
        
        // Test insertion
        assert_eq!(classify_variant_type("A", "ATG"), "INS");
        assert_eq!(classify_variant_type("C", "CTAG"), "INS");
        
        // Test deletion
        assert_eq!(classify_variant_type("ATG", "A"), "DEL");
        assert_eq!(classify_variant_type("CTAG", "C"), "DEL");
        
        // Test complex
        assert_eq!(classify_variant_type("AT", "GCA"), "COMPLEX");
        assert_eq!(classify_variant_type("ATG", "GC"), "COMPLEX");
    }
    
    #[test]
    fn test_allele_match_methods() {
        let ref_match = AlleleMatch::Reference("A".to_string());
        let var_match = AlleleMatch::Variant("G".to_string());
        let other_match = AlleleMatch::Other("T".to_string());
        let no_span = AlleleMatch::NoSpan;
        
        // Test is_definitive_match
        assert!(ref_match.is_definitive_match());
        assert!(var_match.is_definitive_match());
        assert!(!other_match.is_definitive_match());
        assert!(!no_span.is_definitive_match());
        
        // Test is_variant_match
        assert!(!ref_match.is_variant_match());
        assert!(var_match.is_variant_match());
        assert!(!other_match.is_variant_match());
        assert!(!no_span.is_variant_match());
        
        // Test is_reference_match
        assert!(ref_match.is_reference_match());
        assert!(!var_match.is_reference_match());
        assert!(!other_match.is_reference_match());
        assert!(!no_span.is_reference_match());
    }
    
    #[test]
    fn test_snv_analysis() {
        // SNV example: A -> G
        let record = create_test_bam_record("A", "G", 100, "10M");
        let result = analyze_snv_from_cigar(&record, 100, "A", "G");
        assert_eq!(result, AlleleMatch::Variant("G".to_string()));
        
        // Reverse SNV example: T -> C
        let record = create_test_bam_record("T", "C", 200, "10M");
        let result = analyze_snv_from_cigar(&record, 200, "T", "C");
        assert_eq!(result, AlleleMatch::Variant("C".to_string()));
        
        // SNV with reference match
        let record = create_test_bam_record("A", "A", 300, "10M");
        let result = analyze_snv_from_cigar(&record, 300, "A", "G");
        assert_eq!(result, AlleleMatch::Reference("A".to_string()));
    }
    
    #[test]
    fn test_deletion_analysis() {
        // Deletion example: ACG -> AG
        let record = create_test_bam_record("ACG", "AG", 100, "3M1D");
        let result = analyze_deletion_from_cigar(&record, 100, "ACG", "AG");
        assert_eq!(result, AlleleMatch::Variant("DEL:ACG>AG".to_string()));
        
        // Deletion with reference match
        let record = create_test_bam_record("ACGT", "AGT", 200, "4M");
        let result = analyze_deletion_from_cigar(&record, 200, "ACGT", "AGT");
        assert_eq!(result, AlleleMatch::Reference("ACGT".to_string()));
    }
    
    #[test]
    fn test_insertion_analysis() {
        // Insertion example: ACG -> ACGT
        let record = create_test_bam_record("ACG", "ACGT", 100, "3M1I");
        let result = analyze_insertion_from_cigar(&record, 100, "ACG", "ACGT");
        assert_eq!(result, AlleleMatch::Variant("INS:ACG>ACGT".to_string()));
        
        // Insertion with reference match
        let record = create_test_bam_record("ACGT", "ACG", 200, "4M");
        let result = analyze_insertion_from_cigar(&record, 200, "ACGT", "ACG");
        assert_eq!(result, AlleleMatch::Reference("ACGT".to_string()));
    }
}

/// Create a test BAM record with specified input and output alleles
fn create_test_bam_record(input_allele: &str, output_allele: &str, position: usize, cigar: &str) -> bam::Record {
    let mut record = bam::Record::default();
    record.set_reference_sequence_id(0);
    record.set_alignment_start(Some(Ok(position as i32)));
    record.set_cigar(cigar.parse().unwrap());
    record.set_sequence(input_allele.as_bytes().to_vec().into());
    record.set_quality_scores(vec![30; input_allele.len()].into());
    record.set_flags(0);
    
    // Set the mate information for proper pairing
    record.set_next_reference_sequence_id(0);
    record.set_next_alignment_start(Some(Ok((position + 1) as i32)));
    record.set_template_length(1);
    
    record
}

/// Create a test BAM record for SNV analysis
fn create_test_bam_record_snv(position: usize, ref_allele: &str, alt_allele: &str) -> bam::Record {
    let cigar = if ref_allele.len() == alt_allele.len() {
        format!("{}M", ref_allele.len())
    } else {
        format!("{}M", ref_allele.len() - 1)
    };
    
    create_test_bam_record(ref_allele, alt_allele, position, &cigar)
}

/// Create a test BAM record for deletion analysis
fn create_test_bam_record_deletion(position: usize, ref_allele: &str, alt_allele: &str) -> bam::Record {
    let deletion_length = ref_allele.len() - alt_allele.len();
    let cigar = format!("{}M{}D", alt_allele.len(), deletion_length);
    
    create_test_bam_record(ref_allele, alt_allele, position, &cigar)
}

/// Create a test BAM record for insertion analysis
fn create_test_bam_record_insertion(position: usize, ref_allele: &str, alt_allele: &str) -> bam::Record {
    let insertion_length = alt_allele.len() - ref_allele.len();
    let cigar = format!("{}M{}I", ref_allele.len(), insertion_length);
    
    create_test_bam_record(ref_allele, alt_allele, position, &cigar)
}
