use std::fs::File;
use std::io::Read;
use flate2::read::GzDecoder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open("test_data/192929_SISJ1196_T_WGSnanoadap_ALL_46102_0h_1h.vcf.gz")?;
    let mut gz = GzDecoder::new(file);
    let mut buffer = Vec::new();
    
    let bytes_read = gz.read_to_end(&mut buffer)?;
    println!("Read {} bytes", bytes_read);
    println!("Buffer length: {}", buffer.len());
    
    // Print first 100 bytes
    let preview = String::from_utf8_lossy(&buffer[..std::cmp::min(100, buffer.len())]);
    println!("First 100 bytes: {}", preview);
    
    // Count lines
    let content = String::from_utf8_lossy(&buffer);
    let line_count = content.lines().count();
    println!("Total lines: {}", line_count);
    
    Ok(())
}
