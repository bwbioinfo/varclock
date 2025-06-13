# VarClock - Bioinformatics Variant Timing Scanner

A Rust-based tool for analyzing genetic variants in genomic regions with timing information from nanopore sequencing data.

## Overview

VarClock integrates data from three key bioinformatics file formats to provide comprehensive variant analysis:

- **BED files**: Define genomic regions of interest
- **VCF files**: Contain genetic variant information  
- **BAM files**: Store aligned sequencing reads with timing data

## Features

- ✅ **Multi-format Integration**: Seamlessly combines BED, VCF, and BAM data
- ✅ **Timestamp Extraction**: Advanced parsing of nanopore timing data from BAM auxiliary tags
- ✅ **Variant Overlap Detection**: Identifies which reads contain specific variants
- ✅ **Comprehensive Output**: TSV format with detailed variant and timing information
- ✅ **Error Handling**: Robust error reporting and graceful degradation
- ✅ **Memory Efficient**: Optimized for large genomic datasets

## Installation

### Prerequisites

- Rust 1.70+ 
- Cargo package manager

### Build from Source

```bash
git clone https://github.com/username/varclock.git
cd varclock
cargo build --release
```

## Usage

### Basic Command

```bash
./target/release/varclock \
    --bed regions.bed \
    --vcf variants.vcf \
    --bam reads.bam \
    --output results.tsv
```

### Command Line Options

```
Options:
  -b, --bed <BED>        Input BED file containing genomic regions
  -v, --vcf <VCF>        Input VCF file containing genetic variants  
  -a, --bam <BAM>        Input BAM file containing sequencing reads
  -o, --output <OUTPUT>  Output TSV file for results
  -h, --help             Print help information
  -V, --version          Print version information
```

## Input File Formats

### BED File (Genomic Regions)
```
chr1    1000    2000
chr2    5000    6000
chr3    10000   11000
```

### VCF File (Variants)
```
##fileformat=VCFv4.2
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    1500    .       A       G       60      PASS    .
chr2    5500    .       C       T       45      PASS    .
```

### BAM File (Sequencing Reads)
Standard BAM format with optional timing information in auxiliary tags:
- `st:Z:<timestamp>` - Start time in various formats

## Output Format

Tab-separated values (TSV) with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| read_id | Unique read identifier | read_001 |
| timestamp | Sequencing start time | 2021-02-15T14:30:05.123Z |
| contains_variant | Boolean overlap indicator | true |
| variant_description | Human-readable variant | A>G |
| variant_type | Type of genetic variant | SNV |
| region | Genomic coordinates | chr1:1000-2000 |
| region_name | BED region identifier | NA |

### Example Output
```tsv
read_id	timestamp	contains_variant	variant_description	variant_type	region	region_name
read001	2021-02-15T14:30:05.123Z	true	A>G	SNV	chr1:1000-2000	NA
read002	1613397005.123	false	C>T	SNV	chr1:1000-2000	NA
read003	NA	true	G>A	SNV	chr2:5000-6000	NA
```

## Timestamp Extraction

VarClock includes advanced timestamp extraction capabilities for nanopore sequencing data:

### Supported Formats
- ISO 8601: `2021-02-15T14:30:05.123Z`
- Unix timestamp: `1613397005.123`
- MinKNOW format: `2021-02-15 14:30:05.123+00:00`
- Numeric variants: Integer and float timestamps

### Error Handling
- `NA`: No timestamp available
- `INVALID_UTF8`: String parsing error
- `INT32_<value>`: Integer timestamp
- `PARSE_ERROR_<details>`: Extraction error

See [TIMESTAMP_EXTRACTION.md](TIMESTAMP_EXTRACTION.md) for detailed documentation.

## Performance Considerations

### Current Implementation
- Reads entire files for each analysis (simple but inefficient)
- Suitable for small to medium datasets
- Memory usage optimized with record reuse

### Future Optimizations
- BAM/VCF indexing for random access
- Parallel processing for large datasets
- Memory mapping for huge files

## Dependencies

Built with the [noodles](https://github.com/zaeleus/noodles) bioinformatics library ecosystem:

- `noodles-bam`: BAM file parsing
- `noodles-bed`: BED file parsing  
- `noodles-vcf`: VCF file parsing
- `noodles-sam`: SAM auxiliary data access
- `clap`: Command-line interface

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use VarClock in your research, please cite:

```
VarClock: A Rust-based tool for variant timing analysis in nanopore sequencing data
```

## Support

For issues, questions, or feature requests, please open an issue on GitHub.