//! Variant analysis and sequence comparison
//!
//! This module provides functionality for:
//! - Checking if reads span variant positions
//! - Performing sequence-based variant detection using CIGAR and read sequences
//! - Extracting bases from BAM records at specific positions

use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;

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
            AlleleMatch::Reference(seq) => format!("REFERENCE:{seq}"),
            AlleleMatch::Variant(seq) => format!("VARIANT:{seq}"),
            AlleleMatch::Other(seq) => format!("OTHER:{seq}"),
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
pub fn extract_read_bases_at_position(
    sequence: &bam::record::Sequence,
    start_pos: usize,
    length: usize,
) -> String {
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
                    if (32..=126).contains(&base) {
                        base as char
                    } else {
                        'N'
                    }
                }
            });
        } else {
            bases.push('N'); // Position beyond read end
        }
    }

    bases
}

/// Determine the type of variant based on allele sequences and patterns
///
/// **Returns:**
/// - `"SNV"`: Single nucleotide variant
/// - `"MNV"`: Multi-nucleotide variant
/// - `"INS"`: Insertion
/// - `"DEL"`: Deletion
/// - `"DUP"`: Duplication
/// - `"BND"`: Breakend/translocation
/// - `"COMPLEX"`: Complex variant
pub fn classify_variant_type(ref_allele: &str, alt_allele: &str) -> &'static str {
    // Check for breakend variants (BND) - characterized by bracket notation
    if alt_allele.contains('[') || alt_allele.contains(']') {
        return "BND";
    }

    // Check for symbolic alleles indicating structural variants
    if alt_allele.starts_with('<') && alt_allele.ends_with('>') {
        match alt_allele {
            "<DUP>" | "<TDUP>" | "<DUP:TANDEM>" => return "DUP",
            "<DEL>" => return "DEL",
            "<INS>" => return "INS",
            "<INV>" => return "INV",
            "<CNV>" => return "CNV",
            _ => return "COMPLEX",
        }
    }

    // Check for complex breakend notation with angle brackets
    if alt_allele.contains('<') && alt_allele.contains('>') && alt_allele.contains(':') {
        return "BND";
    }

    // Check for duplications based on sequence patterns
    // Tandem duplications often show the reference sequence repeated in the alt
    if alt_allele.len() > ref_allele.len() && alt_allele.starts_with(ref_allele) {
        let inserted_part = &alt_allele[ref_allele.len()..];
        // Check if the inserted part is composed of repetitions of the reference
        // For example: A -> AAA (inserted part "AA" is "A" repeated)
        if inserted_part.len() >= ref_allele.len() {
            // Check if inserted part is made of reference repetitions
            let mut is_duplication = true;
            let ref_bytes = ref_allele.as_bytes();
            let inserted_bytes = inserted_part.as_bytes();

            for chunk_start in (0..inserted_bytes.len()).step_by(ref_bytes.len()) {
                let chunk_end = (chunk_start + ref_bytes.len()).min(inserted_bytes.len());
                let chunk = &inserted_bytes[chunk_start..chunk_end];

                if chunk != ref_bytes {
                    is_duplication = false;
                    break;
                }
            }
            if is_duplication {
                return "DUP";
            }
        } else {
            // For shorter insertions, check if they match the beginning of reference
            if ref_allele.starts_with(inserted_part) {
                return "DUP";
            }
        }
    }

    // Standard size-based classification
    match (ref_allele.len(), alt_allele.len()) {
        (r, a) if r == a && r == 1 => "SNV",
        (r, a) if r == a && r > 1 => "MNV", // Multi-nucleotide variant
        (1, a) if a > 1 => "INS",           // Simple insertion: single base -> multiple bases
        (r, 1) if r > 1 => "DEL",           // Simple deletion: multiple bases -> single base
        _ => "COMPLEX",                     // Complex variants: everything else
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
    for op in cigar.iter().flatten() {
        let op_len = op.len();
        match op.kind() {
            Kind::Match
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
            | Kind::Deletion
            | Kind::Skip => {
                reference_pos += op_len;
            }
            _ => {} // Insertions, clipping, etc. don't advance reference position
        }
    }
    let read_end = if reference_pos > read_start {
        reference_pos - 1
    } else {
        read_start
    };

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
        "DUP" => {
            // Use CIGAR-based analysis for duplications
            analyze_duplication_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        "BND" => {
            // Use CIGAR-based analysis for breakends
            analyze_breakend_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        "INV" | "CNV" => {
            // Use breakend analysis for inversions and CNVs (similar complex patterns)
            analyze_breakend_from_cigar(record, variant_pos, ref_allele, alt_allele, debug)
        }
        _ => {
            // Complex variants and other types - try sequence-based analysis
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

    for op in cigar.iter().flatten() {
        let op_len = op.len();

        // Check if we've reached the variant position
        if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    // Calculate exact position within the operation
                    let offset_in_op = variant_pos - ref_pos;
                    let read_variant_pos = read_pos + offset_in_op;

                    // Extract bases at variant position
                    let read_bases = extract_read_bases_at_position(
                        &read_sequence,
                        read_variant_pos,
                        ref_allele.len(),
                    );

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
        println!(
            "        DEBUG: Analyzing deletion {ref_allele}>{alt_allele} at pos {variant_pos} in read starting at {read_start}"
        );
    }

    for op in cigar.iter().flatten() {
        let op_len = op.len();

        if debug {
            println!(
                "        DEBUG: CIGAR op {:?}({}) at ref_pos={}, read_pos={}",
                op.kind(),
                op_len,
                ref_pos,
                read_pos
            );
        }

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Check if variant position falls within this match
                if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                    let offset_in_op = variant_pos - ref_pos;
                    let read_variant_pos = read_pos + offset_in_op;

                    // For deletion like TACAA>T, check if we see the alt allele (T) at variant position
                    let read_base = extract_read_bases_at_position(
                        &read_sequence,
                        read_variant_pos,
                        alt_allele.len(),
                    );

                    if debug {
                        println!(
                            "        DEBUG: Found match at variant pos, read_base='{read_base}', alt_allele='{alt_allele}'"
                        );
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
                        println!(
                            "        DEBUG: Found deletion after alt base at ref position {ref_pos}"
                        );
                        println!(
                            "        DEBUG: Deletion length: {op_len}, expected length: {expected_deletion_length}"
                        );
                    }

                    if op_len >= expected_deletion_length {
                        return AlleleMatch::Variant(format!("DEL:{ref_allele}>{alt_allele}"));
                    } else {
                        return AlleleMatch::Other(format!(
                            "DEL:{op_len}bp@{ref_pos}_expected:{expected_deletion_length}bp"
                        ));
                    }
                }

                // Check if this deletion overlaps with the variant position (alternative detection)
                let deletion_start = ref_pos;
                let deletion_end = ref_pos + op_len - 1;

                if variant_pos >= deletion_start && variant_pos <= deletion_end {
                    // The variant position itself is deleted
                    if op_len >= expected_deletion_length {
                        return AlleleMatch::Variant(format!("DEL:{ref_allele}>{alt_allele}"));
                    } else {
                        return AlleleMatch::Other(format!(
                            "DEL:{op_len}bp@{deletion_start}_expected:{expected_deletion_length}bp"
                        ));
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

    // If we found the alt base but no deletion, check if this might be a reference match
    if found_alt_base {
        // Extract sequence to see if we have the full reference or just the alt
        let read_base = extract_read_bases_at_position(
            &read_sequence,
            alt_base_read_pos,
            ref_allele
                .len()
                .min(read_sequence.len() - alt_base_read_pos),
        );

        if read_base == ref_allele {
            return AlleleMatch::Reference(read_base);
        } else if read_base.starts_with(alt_allele) && read_base.len() < ref_allele.len() {
            // We have the alt allele but not the full ref, suggesting deletion
            return AlleleMatch::Variant(format!("DEL:{ref_allele}>{alt_allele}"));
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
        println!(
            "        DEBUG: Analyzing insertion {ref_allele}>{alt_allele} at pos {variant_pos} in read starting at {read_start}"
        );
    }

    for op in cigar.iter().flatten() {
        let op_len = op.len();

        if debug {
            println!(
                "        DEBUG: CIGAR op {:?}({}) at ref_pos={}, read_pos={}",
                op.kind(),
                op_len,
                ref_pos,
                read_pos
            );
        }

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // Check if variant position falls within this match
                if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                    let offset_in_op = variant_pos - ref_pos;
                    let read_variant_pos = read_pos + offset_in_op;

                    // Check if the reference base matches
                    let read_base = extract_read_bases_at_position(
                        &read_sequence,
                        read_variant_pos,
                        ref_allele.len(),
                    );

                    if debug {
                        println!(
                            "        DEBUG: Found match at variant pos, read_base='{read_base}', ref_allele='{ref_allele}'"
                        );
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
                        println!(
                            "        DEBUG: Found insertion after reference base at ref position {ref_pos}"
                        );
                    }

                    // Extract the inserted sequence
                    let inserted_bases =
                        extract_read_bases_at_position(&read_sequence, read_pos, op_len);

                    if debug {
                        println!(
                            "        DEBUG: Inserted bases: '{inserted_bases}', expected length: {expected_insertion_length}"
                        );
                    }

                    // For VCF format insertions like G>GATC, the expected insertion is "ATC" (after G)
                    let expected_insertion = if alt_allele.len() > ref_allele.len() {
                        &alt_allele[ref_allele.len()..]
                    } else {
                        ""
                    };

                    if debug {
                        println!(
                            "        DEBUG: Expected insertion sequence: '{expected_insertion}'"
                        );
                    }

                    if op_len >= expected_insertion_length && !expected_insertion.is_empty() {
                        // Check if the insertion matches what we expect
                        if inserted_bases.starts_with(expected_insertion) {
                            return AlleleMatch::Variant(format!("INS:{ref_allele}>{alt_allele}"));
                        } else {
                            return AlleleMatch::Other(format!(
                                "INS:{op_len}bp@{ref_pos}_found:{inserted_bases}"
                            ));
                        }
                    } else if op_len > 0 {
                        // Any insertion at this position that doesn't match expected
                        return AlleleMatch::Other(format!(
                            "INS:{op_len}bp@{ref_pos}_found:{inserted_bases}"
                        ));
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

    // If we found the reference base but no insertion, it's a reference match
    if found_reference_base {
        let read_base = extract_read_bases_at_position(
            &read_sequence,
            reference_base_read_pos,
            ref_allele.len(),
        );
        return AlleleMatch::Reference(read_base);
    }

    AlleleMatch::Indeterminate
}

/// Analyze duplication variants using CIGAR operations
fn analyze_duplication_from_cigar(
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

    if debug {
        println!(
            "        DEBUG: Analyzing duplication {ref_allele}>{alt_allele} at pos {variant_pos} in read starting at {read_start}"
        );
    }

    // For duplications, we expect to see the duplicated sequence in the read
    // Handle symbolic duplications like <DUP>
    if alt_allele == "<DUP>" || alt_allele == "<TDUP>" || alt_allele == "<DUP:TANDEM>" {
        // For symbolic duplications, we can't easily verify from sequence alone
        // Use CIGAR to look for insertions or complex patterns
        let mut ref_pos = read_start;
        let mut read_pos = 0usize;

        for op in cigar.iter().flatten() {
            let op_len = op.len();

            if debug {
                println!(
                    "        DEBUG: CIGAR op {:?}({}) at ref_pos={}, read_pos={}",
                    op.kind(),
                    op_len,
                    ref_pos,
                    read_pos
                );
            }

            // Check if we're in the variant region
            if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        // In a duplication, we might see normal matching
                        return AlleleMatch::Reference(ref_allele.to_string());
                    }
                    Kind::Insertion => {
                        // Duplications often manifest as insertions in the read
                        return AlleleMatch::Variant(format!("DUP:{alt_allele}"));
                    }
                    _ => {}
                }
            }

            // Advance positions
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += op_len;
                }
                Kind::Insertion | Kind::SoftClip => {
                    read_pos += op_len;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }

        // For symbolic duplications, default to indeterminate if we can't determine
        return AlleleMatch::Indeterminate;
    }

    // For sequence-based duplications (e.g., A > ATAT), analyze the sequence
    if alt_allele.len() > ref_allele.len() {
        // Map to variant position and extract sequence
        let mut ref_pos = read_start;
        let mut read_pos = 0usize;

        for op in cigar.iter().flatten() {
            let op_len = op.len();

            if ref_pos <= variant_pos && variant_pos < ref_pos + op_len {
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        let offset_in_op = variant_pos - ref_pos;
                        let read_variant_pos = read_pos + offset_in_op;

                        // Extract sequence at variant position
                        let read_bases = extract_read_bases_at_position(
                            &read_sequence,
                            read_variant_pos,
                            alt_allele.len(),
                        );

                        if debug {
                            println!(
                                "        DEBUG: Found duplication candidate, read_bases='{read_bases}', alt_allele='{alt_allele}'"
                            );
                        }

                        if read_bases == alt_allele {
                            return AlleleMatch::Variant(read_bases);
                        } else if read_bases.starts_with(ref_allele) {
                            return AlleleMatch::Reference(ref_allele.to_string());
                        } else {
                            return AlleleMatch::Other(read_bases);
                        }
                    }
                    _ => {}
                }
            }

            // Advance positions
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    ref_pos += op_len;
                    read_pos += op_len;
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += op_len;
                }
                Kind::Insertion | Kind::SoftClip => {
                    read_pos += op_len;
                }
                Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    AlleleMatch::Indeterminate
}

/// Analyze breakend (BND) variants
fn analyze_breakend_from_cigar(
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

    if debug {
        println!(
            "        DEBUG: Analyzing breakend {ref_allele}>{alt_allele} at pos {variant_pos} in read starting at {read_start}"
        );
    }

    // Breakend variants are complex structural variants that involve chromosome rearrangements
    // They are typically represented with special notation like:
    // - t[chr:pos[ (translocation)
    // - ]chr:pos]t (translocation)
    // - t<chr:pos> (other complex rearrangement)

    // For breakends, we rely on supplementary alignment (SA) tag validation
    // which provides the most reliable evidence for structural variants

    // Check if read spans the variant position
    let read_end = {
        let cigar = record.cigar();
        let mut reference_pos = read_start;
        for op in cigar.iter().flatten() {
            let op_len = op.len();
            match op.kind() {
                Kind::Match
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
                | Kind::Deletion
                | Kind::Skip => {
                    reference_pos += op_len;
                }
                _ => {}
            }
        }
        if reference_pos > read_start {
            reference_pos - 1
        } else {
            read_start
        }
    };

    if variant_pos < read_start || variant_pos > read_end {
        return AlleleMatch::NoSpan;
    }

    // Check supplementary alignment (SA) tag for breakend validation
    let sa_tag_matches = check_sa_tag_for_breakend(record, alt_allele, debug);

    // For breakends, rely exclusively on SA tag validation
    if sa_tag_matches {
        if debug {
            println!("        DEBUG: Breakend evidence found (SA_matches=true)");
        }
        return AlleleMatch::Variant(format!("BND:{alt_allele}"));
    }

    // If no breakend evidence, assume reference
    AlleleMatch::Reference(ref_allele.to_string())
}

/// Check if the supplementary alignment (SA) tag matches the breakend variant call
///
/// The SA tag contains information about chimeric alignments where parts of the same read
/// map to different locations. For breakend variants, we check if the SA mapping span
/// contains the breakpoint position.
///
/// SA tag format: "rname,pos,strand,CIGAR,mapQ,NM;"
/// Example: "chr12,11875518,+,1S8285M27D179S,60,192;"
/// Multiple SA entries are separated by semicolons.
fn check_sa_tag_for_breakend(record: &bam::Record, alt_allele: &str, debug: bool) -> bool {
    use noodles_sam::alignment::record::data::field::Tag;

    // Extract the SA tag from the BAM record
    let sa_tag = Tag::from([b'S', b'A']);
    let sa_value = match record.data().get(&sa_tag) {
        Some(Ok(field_value)) => match field_value {
            noodles_sam::alignment::record::data::field::Value::String(sa_bytes) => {
                match std::str::from_utf8(sa_bytes) {
                    Ok(sa_str) => sa_str.to_string(),
                    Err(_) => {
                        if debug {
                            println!("        DEBUG: SA tag contains invalid UTF-8");
                        }
                        return false;
                    }
                }
            }
            _ => {
                if debug {
                    println!("        DEBUG: SA tag is not a string");
                }
                return false;
            }
        },
        Some(Err(_)) => {
            if debug {
                println!("        DEBUG: Error reading SA tag");
            }
            return false;
        }
        None => {
            if debug {
                println!("        DEBUG: No SA tag found");
            }
            return false;
        }
    };

    if debug {
        println!("        DEBUG: SA tag value: {sa_value}");
    }

    // Parse the breakend alt allele to extract mate position and strand
    let expected_mate_info = parse_breakend_mate_position(alt_allele);
    if expected_mate_info.is_none() {
        if debug {
            println!("        DEBUG: Could not parse mate position from alt allele: {alt_allele}");
        }
        return false;
    }

    let (expected_chrom, expected_pos, _expected_strand) = expected_mate_info.unwrap();

    // Parse SA tag entries (multiple entries separated by semicolons)
    for sa_entry in sa_value.split(';') {
        if sa_entry.trim().is_empty() {
            continue;
        }

        let parts: Vec<&str> = sa_entry.split(',').collect();
        if parts.len() < 6 {
            if debug {
                println!("        DEBUG: Invalid SA entry format: {sa_entry}");
            }
            continue;
        }

        let sa_chrom = parts[0].trim();
        let Ok(sa_pos) = parts[1].trim().parse::<usize>() else {
            if debug {
                println!("        DEBUG: Invalid SA position: {}", parts[1]);
            }
            continue;
        };
        let sa_strand = parts[2].trim();
        let sa_cigar = parts[3].trim();

        if debug {
            println!(
                "        DEBUG: SA entry - chrom: {sa_chrom}, pos: {sa_pos}, strand: {sa_strand}, cigar: {sa_cigar}"
            );
        }

        // Check if chromosomes match (case-insensitive comparison)
        if sa_chrom.to_lowercase() != expected_chrom.to_lowercase() {
            if debug {
                println!(
                    "        DEBUG: Chromosome mismatch: '{sa_chrom}' != '{expected_chrom}' (case-insensitive)"
                );
            }
            continue;
        }

        // Calculate the reference span covered by this SA alignment with strand awareness
        if let Some((span_start, span_end)) =
            calculate_reference_span_from_cigar_with_strand(sa_cigar, sa_pos, sa_strand)
        {
            if debug {
                println!(
                    "        DEBUG: SA alignment spans reference positions {span_start}-{span_end} (breakpoint at {expected_pos})"
                );
            }

            // Check if the breakpoint falls within the SA alignment span with tolerance
            let tolerance = 100; // Allow 100bp tolerance around the span
            let tolerant_start = span_start.saturating_sub(tolerance);
            let tolerant_end = span_end + tolerance;

            if expected_pos >= tolerant_start && expected_pos <= tolerant_end {
                if debug {
                    println!(
                        "        DEBUG: MATCH! Breakpoint {expected_pos} falls within SA span {span_start}-{span_end} (with {tolerance}bp tolerance: {tolerant_start}-{tolerant_end})"
                    );
                }
                return true;
            } else if debug {
                println!(
                    "        DEBUG: Breakpoint {expected_pos} outside SA span {span_start}-{span_end} (with tolerance {tolerant_start}-{tolerant_end})"
                );
            }
        } else if debug {
            println!("        DEBUG: Could not parse CIGAR string: {sa_cigar}");
        }
    }

    if debug {
        println!("        DEBUG: SA tag does not match expected breakend mate position");
    }
    false
}

/// Parse the mate position and strand orientation from a breakend alt allele string
///
/// Returns Some((chromosome, position, `expected_strand_compatible`)) if successfully parsed, None otherwise
/// The `expected_strand_compatible` indicates the strand orientation expected at the mate position
pub fn parse_breakend_mate_position(alt_allele: &str) -> Option<(String, usize, char)> {
    // Handle different breakend notation formats with their strand implications:
    // - t[chr:pos[ (piece extending to the right, forward strand at mate)
    // - ]chr:pos]t (piece extending to the left, reverse strand at mate)
    // - t<chr:pos> (other complex rearrangement, forward strand assumed)

    // Look for bracket notation A[chr:pos[
    if let Some(bracket_start) = alt_allele.find('[')
        && let Some(bracket_end) = alt_allele.rfind('[')
        && bracket_start != bracket_end
    {
        let mate_info = &alt_allele[bracket_start + 1..bracket_end];
        if let Some((chrom, pos)) = parse_mate_coordinate(mate_info) {
            // A[chr:pos[ notation - forward strand at mate position
            return Some((chrom, pos, '+'));
        }
    }

    // Look for bracket notation ]chr:pos]T
    if let Some(bracket_start) = alt_allele.find(']')
        && let Some(bracket_end) = alt_allele.rfind(']')
        && bracket_start != bracket_end
    {
        let mate_info = &alt_allele[bracket_start + 1..bracket_end];
        if let Some((chrom, pos)) = parse_mate_coordinate(mate_info) {
            // ]chr:pos]t notation - reverse strand at mate position
            return Some((chrom, pos, '-'));
        }
    }

    // Handle angle bracket notation <chr:pos>
    if let Some(angle_start) = alt_allele.find('<')
        && let Some(angle_end) = alt_allele.find('>')
    {
        let mate_info = &alt_allele[angle_start + 1..angle_end];
        if let Some((chrom, pos)) = parse_mate_coordinate(mate_info) {
            // <chr:pos> notation - assume forward strand
            return Some((chrom, pos, '+'));
        }
    }

    None
}

/// Parse chromosome and position from mate coordinate string (format: "chr:pos")
fn parse_mate_coordinate(mate_info: &str) -> Option<(String, usize)> {
    let parts: Vec<&str> = mate_info.split(':').collect();
    if parts.len() != 2 {
        return None;
    }

    // Handle uppercase chromosome names like CHR5 by converting to lowercase
    let chrom = parts[0].trim().to_lowercase();
    let pos_str = parts[1].trim();

    let Ok(pos) = pos_str.parse::<usize>() else {
        return None;
    };

    Some((chrom, pos))
}

/// Calculate the reference span covered by a CIGAR alignment with strand awareness
///
/// Returns Some((start, end)) where start is the alignment start position
/// and end is the last reference position covered (inclusive).
///
/// NOTE: The current implementation treats forward and reverse strand alignments
/// the same way for span calculation. The strand parameter is preserved for
/// future enhancement when the specific requirements for reverse strand
/// span calculation are clarified.
///
/// Returns None if the CIGAR string is invalid.
fn calculate_reference_span_from_cigar_with_strand(
    cigar_str: &str,
    start_pos: usize,
    strand: &str,
) -> Option<(usize, usize)> {
    if cigar_str.is_empty() {
        return Some((start_pos, start_pos.saturating_sub(1))); // No operations
    }

    let mut reference_consumed = 0;
    let mut i = 0;
    let cigar_bytes = cigar_str.as_bytes();

    while i < cigar_bytes.len() {
        // Parse the number
        let mut j = i;
        while j < cigar_bytes.len() && cigar_bytes[j].is_ascii_digit() {
            j += 1;
        }

        if j == i || j >= cigar_bytes.len() {
            return None; // Invalid format
        }

        let length_str = std::str::from_utf8(&cigar_bytes[i..j]).ok()?;
        let length = length_str.parse::<usize>().ok()?;

        // Parse the operation
        let op = cigar_bytes[j] as char;

        // Check which operations consume reference
        match op {
            'M' | '=' | 'X' | 'D' | 'N' => {
                // Match, sequence match, sequence mismatch, deletion, skip consume reference
                reference_consumed += length;
            }
            'I' | 'S' | 'H' | 'P' => {
                // Insertion, soft clip, hard clip, padding do not consume reference
            }
            _ => {
                return None; // Unknown operation
            }
        }

        i = j + 1;
    }

    if reference_consumed == 0 {
        Some((start_pos, start_pos.saturating_sub(1))) // No reference consumed
    } else {
        // Calculate span based on strand orientation
        if strand == "-" {
            // For negative strand alignments:
            // - The position in SA tag represents the rightmost coordinate of the alignment
            // - We calculate backwards from this position
            // - The span extends from (start_pos - reference_consumed + 1) to start_pos
            let span_start = start_pos.saturating_sub(reference_consumed - 1);
            let span_end = start_pos;
            Some((span_start, span_end))
        } else {
            // For forward strand alignments (default):
            // - The position in SA tag represents the leftmost coordinate
            // - The span extends from start_pos to (start_pos + reference_consumed - 1)
            let span_start = start_pos;
            let span_end = start_pos + reference_consumed - 1;
            Some((span_start, span_end))
        }
    }
}

/// Legacy function for backwards compatibility
pub fn analyze_read_variant_content(
    record: &bam::Record,
    variant_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
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

        // Test duplication
        assert_eq!(classify_variant_type("A", "AAA"), "DUP"); // Tandem duplication pattern
        assert_eq!(classify_variant_type("AT", "ATAT"), "DUP"); // Tandem duplication
        assert_eq!(classify_variant_type("G", "<DUP>"), "DUP"); // Symbolic duplication
        assert_eq!(classify_variant_type("C", "<TDUP>"), "DUP"); // Tandem duplication symbolic

        // Test breakend
        assert_eq!(classify_variant_type("A", "A[chr2:1000["), "BND"); // Breakend notation
        assert_eq!(classify_variant_type("T", "]chr3:5000]T"), "BND"); // Breakend notation
        assert_eq!(classify_variant_type("G", "G<chr4:2000>"), "BND"); // Complex breakend

        // Test symbolic variants
        assert_eq!(classify_variant_type("ATG", "<DEL>"), "DEL"); // Symbolic deletion
        assert_eq!(classify_variant_type("C", "<INS>"), "INS"); // Symbolic insertion
        assert_eq!(classify_variant_type("ATCG", "<INV>"), "INV"); // Inversion
        assert_eq!(classify_variant_type("A", "<CNV>"), "CNV"); // Copy number variant

        // Test complex
        assert_eq!(classify_variant_type("AT", "GCA"), "COMPLEX");
        assert_eq!(classify_variant_type("ATG", "GC"), "COMPLEX");
        assert_eq!(classify_variant_type("G", "<UNKNOWN>"), "COMPLEX"); // Unknown symbolic
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
    fn test_parse_breakend_mate_position() {
        // Test bracket notation with real-world coordinates
        assert_eq!(
            parse_breakend_mate_position("A[chr12:11875518["),
            Some(("chr12".to_string(), 11875518, '+'))
        );
        assert_eq!(
            parse_breakend_mate_position("]chr3:5000]T"),
            Some(("chr3".to_string(), 5000, '-'))
        );

        // Test uppercase chromosome handling
        assert_eq!(
            parse_breakend_mate_position("A]CHR5:10443321]"),
            Some(("chr5".to_string(), 10443321, '-'))
        );
        assert_eq!(
            parse_breakend_mate_position("T[CHR12:11875518["),
            Some(("chr12".to_string(), 11875518, '+'))
        );

        // Test angle bracket notation
        assert_eq!(
            parse_breakend_mate_position("G<chr4:2000000>"),
            Some(("chr4".to_string(), 2000000, '+'))
        );

        // Test invalid formats
        assert_eq!(parse_breakend_mate_position("A>G"), None);
        assert_eq!(parse_breakend_mate_position("INS"), None);
        assert_eq!(parse_breakend_mate_position("A[invalid]"), None);
    }

    #[test]
    fn test_parse_mate_coordinate() {
        // Test valid coordinates with real-world values
        assert_eq!(
            parse_mate_coordinate("chr12:11875518"),
            Some(("chr12".to_string(), 11875518))
        );
        assert_eq!(
            parse_mate_coordinate("chrX:50000000"),
            Some(("chrx".to_string(), 50000000))
        );

        // Test uppercase chromosome names
        assert_eq!(
            parse_mate_coordinate("CHR5:10443321"),
            Some(("chr5".to_string(), 10443321))
        );
        assert_eq!(
            parse_mate_coordinate("CHRX:50000000"),
            Some(("chrx".to_string(), 50000000))
        );

        // Test invalid coordinates
        assert_eq!(parse_mate_coordinate("chr1"), None);
        assert_eq!(parse_mate_coordinate("chr1:invalid"), None);
        assert_eq!(parse_mate_coordinate("invalid:format:here"), None);
    }

    #[test]
    fn test_sa_tag_only_breakend_validation() {
        // Test that breakend analysis uses only SA tag validation
        // This is a conceptual test to verify the simplified logic

        // Test 1: Valid SA tag parsing with real-world coordinates
        let test_cases = vec![
            ("A[chr12:11875518[", "chr12:11875518"),
            ("]chr3:5000]T", "chr3:5000"),
            ("G<chr4:2000000>", "chr4:2000000"),
            ("A]CHR5:10443321]", "chr5:10443321"),
            ("T[CHRX:50000000[", "chrx:50000000"),
        ];

        for (alt_allele, expected) in test_cases {
            let result = parse_breakend_mate_position(alt_allele);
            let expected_parts: Vec<&str> = expected.split(':').collect();
            let expected_chrom = expected_parts[0].to_string();
            let expected_pos = expected_parts[1].parse::<usize>().unwrap();
            // Determine expected strand based on bracket notation
            let expected_strand = match alt_allele {
                s if s.contains("[") && (s.starts_with("[") || s.ends_with("[")) => '+',
                s if s.contains("]") && (s.starts_with("]") || s.ends_with("]")) => '-',
                s if s.contains("<") => '+', // angle brackets default to forward
                _ => '+',
            };

            assert_eq!(
                result,
                Some((expected_chrom, expected_pos, expected_strand))
            );
        }

        // Test 2: Invalid breakend notation
        assert_eq!(parse_breakend_mate_position("A>G"), None);
        assert_eq!(parse_breakend_mate_position("INS"), None);
        assert_eq!(parse_breakend_mate_position("DEL"), None);

        // Test 3: Mate coordinate parsing with real-world coordinates
        assert_eq!(
            parse_mate_coordinate("chr12:11875518"),
            Some(("chr12".to_string(), 11875518))
        );
        assert_eq!(parse_mate_coordinate("invalid"), None);
    }

    #[test]
    fn test_cigar_reference_span_calculation() {
        // Test basic CIGAR operations with strand awareness
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("100M", 1000, "+"),
            Some((1000, 1099))
        );

        // Test complex CIGAR with real data
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("1S8285M27D179S", 11875518, "+"),
            Some((11875518, 11875518 + 8285 + 27 - 1))
        );

        // Test CIGAR with only non-reference consuming operations
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("50S100I50S", 1000, "+"),
            Some((1000, 999)) // No reference consumed, so end < start
        );

        // Test mixed CIGAR operations
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("50M10I30M5D20M", 1000, "+"),
            Some((1000, 1000 + 50 + 30 + 5 + 20 - 1)) // 105 total reference consumed
        );

        // Invalid CIGAR
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("invalid", 1000, "+"),
            None
        );
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("100", 1000, "+"),
            None
        ); // Missing operation
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("M100", 1000, "+"),
            None
        ); // Invalid format

        // Empty CIGAR
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("", 1000, "+"),
            Some((1000, 999)) // No operations, so end < start
        );
    }

    #[test]
    fn test_breakend_span_validation() {
        // Test that breakend validation correctly checks SA alignment spans

        // Case 1: Breakpoint within SA span
        let span = calculate_reference_span_from_cigar_with_strand("100M", 1000, "+");
        assert_eq!(span, Some((1000, 1099)));

        // Breakpoint at 1050 should be within span 1000-1099
        let breakpoint = 1050;
        let (start, end) = span.unwrap();
        let tolerance = 100;
        assert!(breakpoint >= start.saturating_sub(tolerance) && breakpoint <= end + tolerance);

        // Case 2: Breakpoint outside SA span but within tolerance
        let breakpoint = 1150; // 51bp after end of span
        assert!(breakpoint >= start.saturating_sub(tolerance) && breakpoint <= end + tolerance);

        // Case 3: Breakpoint outside SA span and tolerance
        let breakpoint = 1250; // 151bp after end of span
        assert!(!(breakpoint >= start.saturating_sub(tolerance) && breakpoint <= end + tolerance));

        // Case 4: Complex CIGAR example from real data
        let real_span =
            calculate_reference_span_from_cigar_with_strand("1S8285M27D179S", 11875518, "+");
        assert_eq!(real_span, Some((11875518, 11875518 + 8285 + 27 - 1)));

        let (real_start, real_end) = real_span.unwrap();
        // Breakpoint near the beginning of the span
        let breakpoint1 = 11875520;
        assert!(breakpoint1 >= real_start && breakpoint1 <= real_end);

        // Breakpoint near the end of the span
        let breakpoint2 = real_end - 10;
        assert!(breakpoint2 >= real_start && breakpoint2 <= real_end);
    }

    #[test]
    fn test_multiple_sa_alignments_real_example() {
        // Test the exact multi-SA format provided by the user:
        // SA:Z:chr13,98718353,+,9226S9072M112D11S,60,683;chr12,11882001,-,9082S1824M5D7403S,60,34;

        // Test first SA alignment (chr13)
        let span1 =
            calculate_reference_span_from_cigar_with_strand("9226S9072M112D11S", 98718353, "+");
        assert_eq!(span1, Some((98718353, 98718353 + 9072 + 112 - 1))); // 98718353-98727536

        // Test second SA alignment (chr12) - negative strand
        let span2 =
            calculate_reference_span_from_cigar_with_strand("9082S1824M5D7403S", 11882001, "-");
        // For negative strand: start_pos - reference_consumed + 1 = 11882001 - 1829 + 1 = 11880173
        assert_eq!(span2, Some((11880173, 11882001))); // 11880173-11882001

        // Verify breakpoints would fall within these spans
        let (start1, end1) = span1.unwrap();
        let (start2, end2) = span2.unwrap();

        // Breakpoint in first alignment
        assert!(98720000 >= start1 && 98720000 <= end1);

        // Verify spans are correctly sized
        assert_eq!(end1 - start1 + 1, 9072 + 112); // 9184 bp
        assert_eq!(end2 - start2 + 1, 1824 + 5); // 1829 bp
    }

    #[test]
    fn test_strand_aware_breakend_parsing() {
        // Test that breakend notation correctly identifies strand requirements

        // Forward strand notation tests (A[chr:pos[ format)
        assert_eq!(
            parse_breakend_mate_position("A[chr12:11875518["),
            Some(("chr12".to_string(), 11875518, '+'))
        );
        assert_eq!(
            parse_breakend_mate_position("TTGTAGTA[chr1:100000["),
            Some(("chr1".to_string(), 100000, '+'))
        );

        // Reverse strand notation tests (]chr:pos]T format)
        assert_eq!(
            parse_breakend_mate_position("]chr3:5000]T"),
            Some(("chr3".to_string(), 5000, '-'))
        );
        assert_eq!(
            parse_breakend_mate_position("]chr12:11875518]ACGT"),
            Some(("chr12".to_string(), 11875518, '-'))
        );

        // Test different bracket combinations
        assert_eq!(
            parse_breakend_mate_position("G]chr4:2000000]"),
            Some(("chr4".to_string(), 2000000, '-'))
        );
        assert_eq!(
            parse_breakend_mate_position("[chr5:3000000[C"),
            Some(("chr5".to_string(), 3000000, '+'))
        );

        // Test angle bracket notation (defaults to forward)
        assert_eq!(
            parse_breakend_mate_position("N<chr6:4000000>"),
            Some(("chr6".to_string(), 4000000, '+'))
        );

        // Test case insensitive chromosome handling
        assert_eq!(
            parse_breakend_mate_position("A[CHR12:11875518["),
            Some(("chr12".to_string(), 11875518, '+'))
        );
        assert_eq!(
            parse_breakend_mate_position("]CHRX:50000000]T"),
            Some(("chrx".to_string(), 50000000, '-'))
        );
    }

    #[test]
    fn test_strand_aware_span_calculation() {
        // Test that strand information is properly considered in span calculation

        // Forward strand alignment
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("100M", 1000, "+"),
            Some((1000, 1099))
        );

        // Reverse strand alignment - calculates span in reverse
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("100M", 1000, "-"),
            Some((901, 1000)) // 1000 - 100 + 1 = 901
        );

        // Complex CIGAR with forward strand
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("1S8285M27D179S", 11875518, "+"),
            Some((11875518, 11875518 + 8285 + 27 - 1))
        );

        // Same CIGAR with reverse strand - calculates span in reverse
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("1S8285M27D179S", 11875518, "-"),
            Some((11867207, 11875518)) // 11875518 - 8312 + 1 = 11867207
        );

        // Test with insertions and deletions on both strands
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("50M10I30M5D20M", 1000, "+"),
            Some((1000, 1000 + 50 + 30 + 5 + 20 - 1))
        );
        // Reverse strand calculates backwards
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("50M10I30M5D20M", 1000, "-"),
            Some((896, 1000)) // 1000 - 105 + 1 = 896
        );

        // Test simpler CIGAR
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("10S50M5D", 1000, "+"),
            Some((1000, 1000 + 50 + 5 - 1))
        );
        // Reverse strand calculates backwards
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("10S50M5D", 1000, "-"),
            Some((946, 1000)) // 1000 - 55 + 1 = 946
        );
    }

    #[test]
    fn test_comprehensive_strand_span_calculation() {
        // This test documents the current behavior of span calculation for both strands
        // and validates that the implementation correctly handles various CIGAR patterns

        // Test 1: Simple match operations
        // For both strands, position is leftmost coordinate, span extends rightward
        let test_cases = vec![
            // (cigar, start_pos, strand, expected_start, expected_end, description)
            (
                "100M",
                1000,
                "+",
                1000,
                1099,
                "Simple match on forward strand",
            ),
            (
                "100M",
                1000,
                "-",
                901,
                1000,
                "Simple match on reverse strand - span calculated in reverse",
            ),
            // Test 2: Soft clips (don't consume reference)
            (
                "10S90M",
                1000,
                "+",
                1000,
                1089,
                "Forward strand with 5' soft clip",
            ),
            (
                "10S90M",
                1000,
                "-",
                911,
                1000,
                "Reverse strand with soft clip",
            ),
            (
                "90M10S",
                1000,
                "+",
                1000,
                1089,
                "Forward strand with 3' soft clip",
            ),
            (
                "90M10S",
                1000,
                "-",
                911,
                1000,
                "Reverse strand with soft clip",
            ),
            // Test 3: Hard clips (don't consume reference)
            (
                "10H90M",
                1000,
                "+",
                1000,
                1089,
                "Forward strand with hard clip",
            ),
            (
                "10H90M",
                1000,
                "-",
                911,
                1000,
                "Reverse strand with hard clip",
            ),
            // Test 4: Deletions (consume reference)
            (
                "50M10D40M",
                1000,
                "+",
                1000,
                1099,
                "Forward strand with deletion",
            ),
            (
                "50M10D40M",
                1000,
                "-",
                901,
                1000,
                "Reverse strand with deletion",
            ),
            // Test 5: Insertions (don't consume reference)
            (
                "50M10I40M",
                1000,
                "+",
                1000,
                1089,
                "Forward strand with insertion",
            ),
            (
                "50M10I40M",
                1000,
                "-",
                911,
                1000,
                "Reverse strand with insertion",
            ),
            // Test 6: Complex real-world CIGAR
            (
                "1S8285M27D179S",
                11875518,
                "+",
                11875518,
                11883829,
                "Real SA tag CIGAR forward",
            ),
            (
                "1S8285M27D179S",
                11875518,
                "-",
                11867207,
                11875518,
                "Real SA tag CIGAR reverse",
            ),
            // Test 7: Multiple operations
            (
                "9226S9072M112D11S",
                98718353,
                "+",
                98718353,
                98727536,
                "Complex multi-op forward",
            ),
            (
                "9226S9072M112D11S",
                98718353,
                "-",
                98709170,
                98718353,
                "Complex multi-op reverse",
            ),
            // Test 8: Edge cases
            (
                "100S",
                1000,
                "+",
                1000,
                999,
                "Only soft clips - no ref consumed",
            ),
            ("100S", 1000, "-", 1000, 999, "Only soft clips reverse"),
            (
                "100I",
                1000,
                "+",
                1000,
                999,
                "Only insertion - no ref consumed",
            ),
            ("100I", 1000, "-", 1000, 999, "Only insertion reverse"),
        ];

        for (cigar, start_pos, strand, expected_start, expected_end, description) in test_cases {
            let result = calculate_reference_span_from_cigar_with_strand(cigar, start_pos, strand);
            assert_eq!(
                result,
                Some((expected_start, expected_end)),
                "Failed for {}: CIGAR={}, pos={}, strand={}",
                description,
                cigar,
                start_pos,
                strand
            );
        }

        // Test error cases
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("invalid", 1000, "+"),
            None,
            "Invalid CIGAR should return None"
        );
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("100", 1000, "-"),
            None,
            "CIGAR without operation should return None"
        );
        assert_eq!(
            calculate_reference_span_from_cigar_with_strand("M100", 1000, "+"),
            None,
            "CIGAR with operation before number should return None"
        );
    }

    #[test]
    fn test_breakpoint_span_validation_with_strands() {
        // Test that breakpoint validation correctly uses span calculations
        // for both forward and reverse strand alignments

        // Case 1: Forward strand SA alignment
        let sa_cigar = "100M";
        let sa_pos = 1000;
        let sa_strand = "+";
        let breakpoint = 1050;

        let span = calculate_reference_span_from_cigar_with_strand(sa_cigar, sa_pos, sa_strand);
        assert_eq!(span, Some((1000, 1099)));

        // Breakpoint at 1050 is within span 1000-1099
        let (start, end) = span.unwrap();
        assert!(breakpoint >= start && breakpoint <= end);

        // Case 2: Reverse strand SA alignment
        let sa_strand = "-";
        let span = calculate_reference_span_from_cigar_with_strand(sa_cigar, sa_pos, sa_strand);
        assert_eq!(span, Some((901, 1000))); // Reverse calculation

        // Breakpoint at 1050 is now outside the reverse strand span
        let (start, end) = span.unwrap();
        assert!(!(breakpoint >= start && breakpoint <= end)); // Should be outside

        // Case 3: Real-world example with complex CIGAR
        let sa_cigar = "1S8285M27D179S";
        let sa_pos = 11875518;
        let breakpoint = 11880000;

        // Forward strand
        let span_fwd = calculate_reference_span_from_cigar_with_strand(sa_cigar, sa_pos, "+");
        assert_eq!(span_fwd, Some((11875518, 11883829)));
        let (start, end) = span_fwd.unwrap();
        assert!(breakpoint >= start && breakpoint <= end);

        // Reverse strand - different span calculation
        let span_rev = calculate_reference_span_from_cigar_with_strand(sa_cigar, sa_pos, "-");
        assert_eq!(span_rev, Some((11867207, 11875518))); // Calculated in reverse
        assert_ne!(span_rev, span_fwd); // Should be different from forward strand
    }
}
