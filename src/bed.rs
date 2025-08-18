// BED file processing and genomic region management

use crate::types::BedRegion;
use noodles_bed as bed;
use std::{fs::File, io::BufReader, path::PathBuf};

/// Load all BED regions from file
pub fn load_bed_regions(
    bed_path: &PathBuf,
    debug: bool,
) -> Result<Vec<BedRegion>, Box<dyn std::error::Error>> {
    let bed_metadata = std::fs::metadata(bed_path)?;
    if bed_metadata.len() == 0 {
        return Err("BED file is empty!".into());
    }

    let bed_file = File::open(bed_path)?;
    let mut bed_reader = bed::io::Reader::<4, _>::new(BufReader::new(bed_file));
    let mut bed_record = bed::Record::<4>::default();
    let mut regions = Vec::new();

    loop {
        match bed_reader.read_record(&mut bed_record) {
            Ok(0) => break, // EOF
            Ok(_) => {
                let chrom = std::str::from_utf8(bed_record.reference_sequence_name())?.to_string();
                let start = bed_record.feature_start()?;
                let end = bed_record.feature_end().ok_or("Missing end position")??;

                let region_name = if let Some(name_bytes) = bed_record.name() {
                    std::str::from_utf8(name_bytes)
                        .unwrap_or("INVALID_UTF8")
                        .to_string()
                } else {
                    "NA".to_string()
                };

                let start_pos = usize::from(start);
                let end_pos = usize::from(end);
                let region_string = format!("{chrom}:{start_pos}-{end_pos}");

                regions.push(BedRegion {
                    chrom: chrom.clone(),
                    start_pos,
                    end_pos,
                    region_name: region_name.clone(),
                    region_string: region_string.clone(),
                });

                if debug {
                    println!(
                        "  Loaded BED region: {region_string} (chromosome: '{chrom}', start: {start_pos}, end: {end_pos}, name: '{region_name}')"
                    );
                }
            }
            Err(e) => return Err(format!("Error reading BED record: {e:?}").into()),
        }
    }

    if debug {
        println!("Loaded {} BED regions", regions.len());
    }
    Ok(regions)
}
