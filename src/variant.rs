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

/// Analyze read sequence to determine detailed allele content  
/// 
/// This is an enhanced version that provides detailed information about what was found in the read sequence.
///
/// **Parameters:**
/// - `record`: The BAM record containing read data
/// - `variant_pos`: 1-based genomic position of the variant
/// - `ref_allele`: Reference allele sequence
/// - `alt_allele`: Alternative allele sequence
/// - `debug`: Whether to print debug information
/// 
/// **Returns:**
/// - `AlleleMatch` enum with detailed information about the match
pub fn analyze_read_allele_content_detailed(
    record: &bam::Record, 
    variant_pos: usize, 
    ref_allele: &str, 
    alt_allele: &str,
    debug: bool,
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
            analyze_snv_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        "DEL" => {
            // Use CIGAR-based analysis for deletions
            analyze_deletion_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        "INS" => {
            // Use CIGAR-based analysis for insertions
            analyze_insertion_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        _ => {
            // Complex variants - try sequence-based analysis
            analyze_snv_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
    }
}

/// Analyze SNV/MNV variants using CIGAR-based position mapping
fn analyze_snv_from_cigar(
    record: &bam::Record,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    _debug: bool,
) -> AlleleMatch {
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return AlleleMatch::Indeterminate,
    };
    
    let cigar = record.cigar();
    let read_sequence = record.sequence();
    
    // Map reference position to read position using CIGAR
    let mut ref_pos = read_start;
    let mut read_pos = 0usize;
    
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            let op_len = usize::from(op.len());
            
            // Check if we've reached the variant position
            if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        // Calculate exact position within the operation
                        let offset_in_op = variant_pos - ref_pos;
                        let read_variant_pos = read_pos + offset_in_op;
                        
                        // Extract bases at variant position
                        let read_bases = extract_read_bases_at_position(&read_sequence, read_variant_pos, ref_allele.len());
                        
                        // Compare with reference and alternative alleles
                        if read_bases == alt_allele {
                            return AlleleMatch::Variant(read_bases);
                        } else if read_bases == ref_allele {
                            return AlleleMatch::Reference(read_bases);
                        } else {
                            return AlleleMatch::Other(read_bases);
                        }
                    }
                    Kind::Deletion => {
                        // Variant position falls in a deletion - this read supports the deletion
                        return AlleleMatch::Deletion;
                    }
                    _ => {
                        return AlleleMatch::Indeterminate;
                    }
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
    
    AlleleMatch::Indeterminate
}

/// Analyze deletion variants using CIGAR operations
fn analyze_deletion_from_cigar(
    record: &bam::Record,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    debug: bool,
) -> AlleleMatch {
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return AlleleMatch::Indeterminate,
    };
    
    let cigar = record.cigar();
    let read_sequence = record.sequence();
    
    let expected_deletion_length = ref_allele.len() - alt_allele.len();
    
    // Map reference position to read position using CIGAR
    let mut ref_pos = read_start;
    let mut read_pos = 0usize;
    let mut found_alt_base = false;
    let mut alt_base_read_pos = 0;
    
    if debug {
        println!("        DEBUG: Analyzing deletion {}>{} at pos {} in read starting at {}", 
                 ref_allele, alt_allele, variant_pos, read_start);
    }
    
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            let op_len = usize::from(op.len());
            
            if debug {
                println!("        DEBUG: CIGAR op {:?}({}) at ref_pos={}, read_pos={}", 
                         op.kind(), op_len, ref_pos, read_pos);
            }
            
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    // Check if variant position falls within this match
                    if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                        let offset_in_op = variant_pos - ref_pos;
                        let read_variant_pos = read_pos + offset_in_op;
                        
                        // For deletion like TACAA>T, check if we see the alt allele (T) at variant position
                        let read_base = extract_read_bases_at_position(&read_sequence, read_variant_pos, alt_allele.len());
                        
                        if debug {
                            println!("        DEBUG: Found match at variant pos, read_base='{}', alt_allele='{}'", 
                                     read_base, alt_allele);
                        }
                        
                        if read_base == alt_allele {
                            found_alt_base = true;
                            alt_base_read_pos = read_variant_pos;
                            
                            // Continue to look for deletion after this base
                            // Don't return yet - check subsequent CIGAR operations
                        } else if read_base == ref_allele {
                            // Full reference sequence present - no deletion
                            return AlleleMatch::Reference(read_base);
                        } else {
                            return AlleleMatch::Other(read_base);
                        }
                    }
                    
                    // Advance both positions for match
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Deletion => {
                    // If we previously found the alt base and now see a deletion
                    // at the position right after the variant position
                    if found_alt_base && ref_pos == variant_pos + alt_allele.len() {
                        if debug {
                            println!("        DEBUG: Found deletion after alt base at ref position {}", ref_pos);
                            println!("        DEBUG: Deletion length: {}, expected length: {}", 
                                     op_len, expected_deletion_length);
                        }
                        
                        if op_len >= expected_deletion_length {
                            return AlleleMatch::Variant(format!("DEL:{}>{}", ref_allele, alt_allele));
                        } else {
                            return AlleleMatch::Other(format!("DEL:{}bp@{}_expected:{}bp", op_len, ref_pos, expected_deletion_length));
                        }
                    }
                    
                    // Check if this deletion overlaps with the variant position (alternative detection)
                    let deletion_start = ref_pos;
                    let deletion_end = ref_pos + op_len - 1;
                    
                    if variant_pos >= deletion_start && variant_pos <= deletion_end {
                        // The variant position itself is deleted
                        if op_len >= expected_deletion_length {
                            return AlleleMatch::Variant(format!("DEL:{}>{}", ref_allele, alt_allele));
                        } else {
                            return AlleleMatch::Other(format!("DEL:{}bp@{}_expected:{}bp", op_len, deletion_start, expected_deletion_length));
                        }
                    }
                    
                    // Advance reference position for deletion
                    ref_pos += op_len;
                }
                Kind::Insertion => {
                    // Insertions don't advance reference position
                    read_pos += op_len;
                }
                Kind::Skip => {
                    ref_pos += op_len;
                }
                Kind::SoftClip => {
                    read_pos += op_len;
                }
                Kind::HardClip | Kind::Pad => {
                    // Neither position advances
                }
            }
        }
    }
    
    // If we found the alt base but no deletion, check if this might be a reference match
    if found_alt_base {
        // Extract sequence to see if we have the full reference or just the alt
        let read_base = extract_read_bases_at_position(&read_sequence, alt_base_read_pos, ref_allele.len().min(read_sequence.len() - alt_base_read_pos));
        
        if read_base == ref_allele {
            return AlleleMatch::Reference(read_base);
        } else if read_base.starts_with(alt_allele) && read_base.len() < ref_allele.len() {
            // We have the alt allele but not the full ref, suggesting deletion
            return AlleleMatch::Variant(format!("DEL:{}>{}", ref_allele, alt_allele));
        } else {
            return AlleleMatch::Other(read_base);
        }
    }
    
    AlleleMatch::Indeterminate
}

/// Analyze insertion variants using CIGAR operations
fn analyze_insertion_from_cigar(
    record: &bam::Record,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    debug: bool,
) -> AlleleMatch {
    let read_start = match record.alignment_start() {
        Some(Ok(start)) => usize::from(start),
        _ => return AlleleMatch::Indeterminate,
    };
    
    let cigar = record.cigar();
    let read_sequence = record.sequence();
    
    let expected_insertion_length = alt_allele.len() - ref_allele.len();
    
    // Map reference position to read position using CIGAR
    let mut ref_pos = read_start;
    let mut read_pos = 0usize;
    let mut found_reference_base = false;
    let mut reference_base_read_pos = 0;
    
    if debug {
        println!("        DEBUG: Analyzing insertion {}>{} at pos {} in read starting at {}", 
                 ref_allele, alt_allele, variant_pos, read_start);
    }
    
    for operation in cigar.iter() {
        if let Ok(op) = operation {
            let op_len = usize::from(op.len());
            
            if debug {
                println!("        DEBUG: CIGAR op {:?}({}) at ref_pos={}, read_pos={}", 
                         op.kind(), op_len, ref_pos, read_pos);
            }
            
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    // Check if variant position falls within this match
                    if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                        let offset_in_op = variant_pos - ref_pos;
                        let read_variant_pos = read_pos + offset_in_op;
                        
                        // Check if the reference base matches
                        let read_base = extract_read_bases_at_position(&read_sequence, read_variant_pos, ref_allele.len());
                        
                        if debug {
                            println!("        DEBUG: Found match at variant pos, read_base='{}', ref_allele='{}'", 
                                     read_base, ref_allele);
                        }
                        
                        if read_base == ref_allele {
                            found_reference_base = true;
                            reference_base_read_pos = read_variant_pos;
                            
                            // Continue to look for insertion after this base
                            // Don't return yet - check subsequent CIGAR operations
                        } else {
                            return AlleleMatch::Other(read_base);
                        }
                    }
                    
                    // Advance both positions for match
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Insertion => {
                    // If we previously found the reference base and now see an insertion
                    // at the position right after the variant position
                    if found_reference_base && ref_pos == variant_pos + ref_allele.len() {
                        if debug {
                            println!("        DEBUG: Found insertion after reference base at ref position {}", ref_pos);
                        }
                        
                        // Extract the inserted sequence
                        let inserted_bases = extract_read_bases_at_position(&read_sequence, read_pos, op_len);
                        
                        if debug {
                            println!("        DEBUG: Inserted bases: '{}', expected length: {}", 
                                     inserted_bases, expected_insertion_length);
                        }
                        
                        // For VCF format insertions like G>GATC, the expected insertion is "ATC" (after G)
                        let expected_insertion = if alt_allele.len() > ref_allele.len() {
                            &alt_allele[ref_allele.len()..]
                        } else {
                            ""
                        };
                        
                        if debug {
                            println!("        DEBUG: Expected insertion sequence: '{}'", expected_insertion);
                        }
                        
                        if op_len >= expected_insertion_length && !expected_insertion.is_empty() {
                            // Check if the insertion matches what we expect
                            if inserted_bases.starts_with(expected_insertion) {
                                return AlleleMatch::Variant(format!("INS:{}>{}", ref_allele, alt_allele));
                            } else {
                                return AlleleMatch::Other(format!("INS:{}bp@{}_found:{}", op_len, ref_pos, inserted_bases));
                            }
                        } else if op_len > 0 {
                            // Any insertion at this position that doesn't match expected
                            return AlleleMatch::Other(format!("INS:{}bp@{}_found:{}", op_len, ref_pos, inserted_bases));
                        }
                    }
                    
                    // Advance read position for insertion (ref_pos doesn't advance for insertions)
                    read_pos += op_len;
                }
                Kind::Deletion => {
                    ref_pos += op_len;
                }
                Kind::Skip => {
                    ref_pos += op_len;
                }
                Kind::SoftClip => {
                    read_pos += op_len;
                }
                Kind::HardClip | Kind::Pad => {
                    // Neither position advances
                }
            }
        }
    }
    
    // If we found the reference base but no insertion, it's a reference match
    if found_reference_base {
        let read_base = extract_read_bases_at_position(&read_sequence, reference_base_read_pos, ref_allele.len());
        return AlleleMatch::Reference(read_base);
    }
    
    AlleleMatch::Indeterminate
}

/// Legacy function for backwards compatibility
pub fn analyze_read_variant_content(
    record: &bam::Record, 
    variant_pos: usize, 
    ref_allele: &str, 
    alt_allele: &str
) -> Option<bool> {
    match analyze_read_allele_content_detailed(record, variant_pos, ref_allele, alt_allele, false) {
        AlleleMatch::Variant(_) => Some(true),
        AlleleMatch::Reference(_) => Some(false),
        _ => None,
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
}
