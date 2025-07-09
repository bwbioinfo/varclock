//! BAM file processing and read analysis
//!
//! This module provides functionality for:
//! - Querying indexed BAM files for specific genomic regions
//! - Extracting sequencing reads with metadata (timestamps, mapping quality)
//! - Processing read alignments and calculating genomic coordinates

use std::path::PathBuf;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_core;

/// Query BAM file for reads overlapping a genomic region
/// 
/// Returns a vector of tuples containing:
/// - BAM record
/// - Read ID
/// - Read start position (1-based)
/// - Read end position (1-based)
/// - Timestamp (from st tag)
/// - Mapping quality
pub fn query_bam_records_for_region(
    bam_path: &PathBuf,
    chrom: &str,
    start_pos: usize,
    end_pos: usize,
    debug: bool,
) -> Result<Vec<(bam::Record, String, usize, usize, String, u8)>, Box<dyn std::error::Error + Send + Sync>> {
    if debug {
        println!("    DEBUG: BAM query - file: {bam_path:?}, region: {chrom}:{start_pos}-{end_pos}");
    }
    // STEP 1: Initialize indexed BAM reader
    if debug {
        println!("    DEBUG: Initializing BAM reader...");
    }
    let mut indexed_reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    
    if debug {
        println!("    DEBUG: Reading BAM header...");
    }
    let header = indexed_reader.read_header()?;
    
    // STEP 2: Create genomic region specification for query
    let start_position = noodles_core::Position::try_from(start_pos)?;
    let end_position = noodles_core::Position::try_from(end_pos)?;
    let interval = start_position..=end_position;
    let region = noodles_core::Region::new(chrom, interval);
    
    // STEP 3: Execute indexed query
    if debug {
        println!("    DEBUG: Executing BAM query for region {region}...");
    }
    let query = indexed_reader.query(&header, &region)?;
    
    let mut records = Vec::new();
    let mut record_count = 0;
    
    // STEP 4: Process each record and extract both metadata and full record
    for result in query {
        record_count += 1;
        if debug && record_count <= 5 {
            println!("    DEBUG: Processing BAM record #{record_count}");
        }
        
        let record = result?;

        // Extract read identifier
        let read_id = record.name()
            .map(|n| std::str::from_utf8(n).unwrap_or("unknown"))
            .unwrap_or("unknown")
            .to_string();

        // Calculate genomic coordinates
        let (read_start, read_end) = if let Some(Ok(alignment_start)) = record.alignment_start() {
            let start = usize::from(alignment_start);
            let cigar = record.cigar();
            let mut reference_pos = start;

            for op in cigar.iter().flatten() {
                use noodles_sam::alignment::record::cigar::op::Kind;
                let op_len = op.len();

                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip => {
                        reference_pos += op_len;
                    }
                    _ => {}
                }
            }

            let read_end = if reference_pos > start { reference_pos - 1 } else { start };
            (start, read_end)
        } else {
            continue;
        };

        // Extract timestamp
        let timestamp = {
            let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);
            match record.data().get(&st_tag) {
                Some(Ok(field_value)) => {
                    match field_value {
                        sam::alignment::record::data::field::Value::String(ts_bytes) => {
                            std::str::from_utf8(ts_bytes).unwrap_or("INVALID_UTF8").to_string()
                        }
                        sam::alignment::record::data::field::Value::Int8(val) => format!("INT8_{val}"),
                        sam::alignment::record::data::field::Value::UInt8(val) => format!("UINT8_{val}"),
                        sam::alignment::record::data::field::Value::Int16(val) => format!("INT16_{val}"),
                        sam::alignment::record::data::field::Value::UInt16(val) => format!("UINT16_{val}"),
                        sam::alignment::record::data::field::Value::Int32(val) => format!("INT32_{val}"),
                        sam::alignment::record::data::field::Value::UInt32(val) => format!("UINT32_{val}"),
                        sam::alignment::record::data::field::Value::Float(val) => format!("FLOAT_{val}"),
                        _ => format!("UNKNOWN_TYPE_{field_value:?}"),
                    }
                }
                Some(Err(e)) => format!("PARSE_ERROR_{e:?}"),
                None => "NA".to_string(),
            }
        };

        // Extract mapping quality
        let mapping_quality = record.mapping_quality()
            .map(u8::from)
            .unwrap_or(0);

        // Store the full record along with metadata
        if debug && record_count <= 5 {
            // Also show the sequence for debugging
            let sequence_debug = {
                let seq = record.sequence();
                let mut seq_str = String::new();
                let display_len = seq.len().min(30);
                
                for i in 0..display_len {
                    if let Some(base) = seq.get(i) {
                        seq_str.push(match base {
                            // Try ASCII encoding first (most common)
                            b'A' => 'A',
                            b'C' => 'C', 
                            b'G' => 'G',
                            b'T' => 'T',
                            b'N' => 'N',
                            // Also try the standard BAM 4-bit encoding as fallback
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
                                    '?'
                                }
                            },
                        });
                    } else {
                        seq_str.push('?');
                    }
                }
                
                if seq.len() > 30 {
                    format!("{}... (len={})", seq_str, seq.len())
                } else {
                    format!("{} (len={})", seq_str, seq.len())
                }
            };
            
            println!("    DEBUG: Stored record #{record_count}: {read_id} at {chrom}:{read_start}-{read_end} (mapq={mapping_quality}, timestamp={timestamp})");
            println!("    DEBUG: Record #{record_count} sequence: {sequence_debug}");
        }
        records.push((record, read_id, read_start, read_end, timestamp, mapping_quality));
    }

    if debug {
        println!("    DEBUG: BAM query completed - found {} total records for region {}:{}-{}", 
                records.len(), chrom, start_pos, end_pos);
        
        // Show summary of all reads found (not just first 5)
        if records.len() > 5 {
            println!("    DEBUG: Additional {} reads found (showing first 5 above)", records.len() - 5);
            println!("    DEBUG: Read position range: {}:{}-{} to {}:{}-{}", 
                    chrom, 
                    records.iter().map(|(_, _, start, _, _, _)| *start).min().unwrap_or(0),
                    records.iter().map(|(_, _, _, end, _, _)| *end).max().unwrap_or(0),
                    chrom,
                    records.iter().map(|(_, _, start, _, _, _)| *start).max().unwrap_or(0),
                    records.iter().map(|(_, _, _, end, _, _)| *end).max().unwrap_or(0));
        }
    }
    
    Ok(records)
}
