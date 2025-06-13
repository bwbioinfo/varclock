# Timestamp Extraction Enhancement

## Overview

The varclock tool has been enhanced with comprehensive timestamp extraction from BAM files, specifically designed to handle nanopore sequencing data that includes timing information in auxiliary tags.

## Implementation Details

### BAM Auxiliary Tag Format
The tool searches for the `st` (start time) auxiliary tag in the format:
```
st:Z:<timestamp>
```

Where:
- `st`: Two-character tag identifier for "start time"  
- `Z`: BAM type code indicating null-terminated string
- `<timestamp>`: The actual timestamp value

### Supported Timestamp Formats

The implementation handles multiple timestamp formats commonly found in nanopore data:

1. **ISO 8601 Format** (Recommended)
   ```
   st:Z:2021-02-15T14:30:05.123Z
   st:Z:2021-02-15T14:30:05.123+00:00
   ```

2. **Unix Timestamp**
   ```
   st:Z:1613397005.123
   st:i:1613397005
   st:f:1613397005.123
   ```

3. **MinKNOW Format**
   ```
   st:Z:2021-02-15 14:30:05.123+00:00
   ```

### Error Handling

The implementation provides detailed error information for debugging:

- `NA`: No timestamp tag found (normal for reads without timing data)
- `INVALID_UTF8`: String tag found but contains invalid UTF-8
- `INT32_<value>`: Timestamp stored as 32-bit integer
- `FLOAT_<value>`: Timestamp stored as floating-point number
- `PARSE_ERROR_<details>`: Error occurred while parsing tag data
- `UNKNOWN_TYPE_<details>`: Tag found but in unexpected format

### Code Structure

The timestamp extraction follows this pattern:

```rust
// Create tag identifier
let st_tag = sam::alignment::record::data::field::Tag::from([b's', b't']);

// Extract from BAM record auxiliary data
match bam_record.data().get(&st_tag) {
    Some(Ok(field_value)) => {
        // Handle different value types (String, Int32, Float, etc.)
    }
    Some(Err(e)) => {
        // Handle parsing errors
    }
    None => {
        // Handle missing tag
    }
}
```

## Usage

The enhanced tool now properly extracts timestamps from nanopore BAM files:

```bash
./varclock --bed regions.bed --vcf variants.vcf --bam reads.bam --output results.tsv
```

Output includes the extracted timestamp in the second column:

```tsv
read_id	timestamp	contains_variant	variant_description	variant_type	region	region_name
read001	2021-02-15T14:30:05.123Z	true	A>G	SNV	chr1:1000-2000	region1
read002	1613397005.123	false	C>T	SNV	chr1:1000-2000	region1
read003	NA	true	G>A	SNV	chr1:1000-2000	region1
```

## Benefits

1. **Comprehensive Format Support**: Handles multiple timestamp formats
2. **Robust Error Handling**: Provides detailed error information for debugging
3. **Performance Optimized**: Efficient tag lookup and extraction
4. **Standards Compliant**: Follows BAM specification for auxiliary data access
5. **Debugging Friendly**: Clear error messages help identify data issues

## Future Enhancements

Potential improvements for future versions:

1. **Timestamp Parsing**: Convert timestamps to standardized format
2. **Time Zone Handling**: Normalize timestamps to UTC
3. **Duration Calculation**: Calculate sequencing duration from start/end times
4. **Statistical Analysis**: Provide timing statistics in output
5. **Multiple Tag Support**: Extract additional timing tags (end time, etc.)
