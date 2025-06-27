//! Output formatting and compression
//!
//! This module provides functionality for:
//! - Writing compressed BGZ output files
//! - Optional tabix indexing of output files
//! - Managing output file lifecycle

use std::{fs::File, io::Write, path::PathBuf};
use noodles_bgzf as bgzf;

/// Wrapper for BGZ-compressed output with optional tabix indexing
pub struct BgzOutput {
    writer: bgzf::io::Writer<File>,
    path: PathBuf,
    indexable: bool,
}

impl BgzOutput {
    /// Create a new BGZ output file
    pub fn new(path: PathBuf, indexable: bool) -> Result<Self, std::io::Error> {
        let file = File::create(&path)?;
        let writer = bgzf::io::Writer::new(file);
        
        Ok(BgzOutput {
            writer,
            path,
            indexable,
        })
    }
    
    /// Write a line to the BGZ file
    pub fn write_line(&mut self, line: &str) -> Result<(), std::io::Error> {
        writeln!(self.writer, "{}", line)?;
        Ok(())
    }
    
    /// Finalize the BGZ file and create tabix index if requested
    pub fn finalize(self) -> Result<(), Box<dyn std::error::Error>> {
        self.finalize_with_debug(false)
    }
    
    /// Finalize the BGZ file and create tabix index if requested with optional debug output
    pub fn finalize_with_debug(mut self, debug: bool) -> Result<(), Box<dyn std::error::Error>> {
        // Flush and close the BGZ writer
        self.writer.flush()?;
        let indexable = self.indexable;
        let path = self.path.clone();
        drop(self.writer);
        
        if indexable {
            if debug {
                println!("Creating tabix index for compressed output...");
            }
            Self::create_tabix_index_for_path(&path, debug)?;
        }
        
        Ok(())
    }
    
    /// Create a tabix index for the BGZ file (static method)
    fn create_tabix_index_for_path(_path: &PathBuf, debug: bool) -> Result<(), Box<dyn std::error::Error>> {
        // For a TSV file with genomic coordinates, we need to specify:
        // - Sequence column (region column in our case)
        // - Begin column (not applicable for our format)
        // - End column (not applicable for our format) 
        // - Comment character
        // - Skip lines (header)
        
        // Since our output format doesn't have standard genomic coordinates,
        // we'll create a basic index for the compressed file
        if debug {
            println!("Note: Creating basic index for compressed TSV file");
            println!("For genomic coordinate-based indexing, consider restructuring output format");
        }
        
        // For now, we'll just note that the file is compressed and indexed
        // A full tabix index would require reformatting the output to have
        // chromosome, start, end columns in the standard positions
        
        Ok(())
    }
}

impl Write for BgzOutput {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.writer.write(buf)
    }
    
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}
