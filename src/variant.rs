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

    // For breakends, we primarily look for soft/hard clipping or split reads
    // which indicate the read spans a breakpoint

    let cigar = record.cigar();
    let mut has_clipping = false;
    let mut has_complex_pattern = false;

    for op in cigar.iter().flatten() {
        let op_len = op.len();

        match op.kind() {
            Kind::SoftClip | Kind::HardClip => {
                has_clipping = true;
                if debug {
                    println!(
                        "        DEBUG: Found clipping {:?}({}) - potential breakend support",
                        op.kind(),
                        op_len
                    );
                }
            }
            Kind::Skip => {
                // Large skips might indicate structural variants
                if op_len > 1000 {
                    has_complex_pattern = true;
                    if debug {
                        println!(
                            "        DEBUG: Found large skip {op_len}bp - potential breakend support"
                        );
                    }
                }
            }
            _ => {}
        }
    }

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

    // For breakends, the presence of clipping or complex patterns near the variant
    // position suggests the read might support the breakend
    if has_clipping || has_complex_pattern || sa_tag_matches {
        if debug {
            println!(
                "        DEBUG: Breakend evidence found (clipping={has_clipping}, complex={has_complex_pattern}, SA_matches={sa_tag_matches})"
            );
        }
        return AlleleMatch::Variant(format!("BND:{alt_allele}"));
    }

    // If no breakend evidence, assume reference
    AlleleMatch::Reference(ref_allele.to_string())
}

/// Check if the supplementary alignment (SA) tag matches the breakend variant call
///
/// The SA tag contains information about chimeric alignments where parts of the same read
/// map to different locations. For breakend variants, this indicates the read spans a breakpoint.
///
/// SA tag format: "rname,pos,strand,CIGAR,mapQ,NM"
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

    // Parse the breakend alt allele to extract mate position
    let expected_mate_pos = parse_breakend_mate_position(alt_allele);
    if expected_mate_pos.is_none() {
        if debug {
            println!("        DEBUG: Could not parse mate position from alt allele: {alt_allele}");
        }
        return false;
    }

    let (expected_chrom, expected_pos) = expected_mate_pos.unwrap();

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

        let sa_chrom = parts[0];
        let sa_pos = match parts[1].parse::<usize>() {
            Ok(pos) => pos,
            Err(_) => {
                if debug {
                    println!("        DEBUG: Invalid SA position: {}", parts[1]);
                }
                continue;
            }
        };

        if debug {
            println!("        DEBUG: SA entry - chrom: {sa_chrom}, pos: {sa_pos}");
        }

        // Check if SA tag mate position matches the expected breakend location
        if sa_chrom == expected_chrom && sa_pos == expected_pos {
            if debug {
                println!("        DEBUG: SA tag matches breakend mate position!");
            }
            return true;
        }

        // Allow for some tolerance in position matching (within 10bp)
        if sa_chrom == expected_chrom && sa_pos.abs_diff(expected_pos) <= 10 {
            if debug {
                println!("        DEBUG: SA tag matches breakend mate position within tolerance!");
            }
            return true;
        }
    }

    if debug {
        println!("        DEBUG: SA tag does not match expected breakend mate position");
    }
    false
}

/// Parse the mate position from a breakend alt allele string
///
/// Returns Some((chromosome, position)) if successfully parsed, None otherwise
fn parse_breakend_mate_position(alt_allele: &str) -> Option<(String, usize)> {
    // Handle different breakend notation formats:
    // - t[chr:pos[ (translocation)
    // - ]chr:pos]t (translocation)
    // - t<chr:pos> (other complex rearrangement)

    // Look for bracket notation A[chr:pos[
    if let Some(bracket_start) = alt_allele.find('[')
        && let Some(bracket_end) = alt_allele.rfind('[')
        && bracket_start != bracket_end
    {
        let mate_info = &alt_allele[bracket_start + 1..bracket_end];
        return parse_mate_coordinate(mate_info);
    }

    // Look for bracket notation ]chr:pos]T
    if let Some(bracket_start) = alt_allele.find(']')
        && let Some(bracket_end) = alt_allele.rfind(']')
        && bracket_start != bracket_end
    {
        let mate_info = &alt_allele[bracket_start + 1..bracket_end];
        return parse_mate_coordinate(mate_info);
    }

    // Handle angle bracket notation <chr:pos>
    if let Some(angle_start) = alt_allele.find('<')
        && let Some(angle_end) = alt_allele.find('>')
    {
        let mate_info = &alt_allele[angle_start + 1..angle_end];
        return parse_mate_coordinate(mate_info);
    }

    None
}

/// Parse chromosome and position from mate coordinate string (format: "chr:pos")
fn parse_mate_coordinate(mate_info: &str) -> Option<(String, usize)> {
    let parts: Vec<&str> = mate_info.split(':').collect();
    if parts.len() != 2 {
        return None;
    }

    let chrom = parts[0].to_string();
    let pos = match parts[1].parse::<usize>() {
        Ok(p) => p,
        Err(_) => return None,
    };

    Some((chrom, pos))
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
        // Test bracket notation
        assert_eq!(
            parse_breakend_mate_position("A[chr2:1000["),
            Some(("chr2".to_string(), 1000))
        );
        assert_eq!(
            parse_breakend_mate_position("]chr3:5000]T"),
            Some(("chr3".to_string(), 5000))
        );

        // Test angle bracket notation
        assert_eq!(
            parse_breakend_mate_position("G<chr4:2000>"),
            Some(("chr4".to_string(), 2000))
        );

        // Test invalid formats
        assert_eq!(parse_breakend_mate_position("A>G"), None);
        assert_eq!(parse_breakend_mate_position("INS"), None);
        assert_eq!(parse_breakend_mate_position("A[invalid]"), None);
    }

    #[test]
    fn test_parse_mate_coordinate() {
        // Test valid coordinates
        assert_eq!(
            parse_mate_coordinate("chr1:1000"),
            Some(("chr1".to_string(), 1000))
        );
        assert_eq!(
            parse_mate_coordinate("chrX:50000"),
            Some(("chrX".to_string(), 50000))
        );

        // Test invalid coordinates
        assert_eq!(parse_mate_coordinate("chr1"), None);
        assert_eq!(parse_mate_coordinate("chr1:invalid"), None);
        assert_eq!(parse_mate_coordinate("invalid:format:here"), None);
    }
}
