# VarClock - Bioinformatics Variant Timing Scanner

A Rust-based tool for analyzing genetic variants in genomic regions with timing information from nanopore sequencing data.

## Overview

VarClock integrates data from three key bioinformatics file formats to provide comprehensive variant analysis with full support for multi-allelic variants:

- **BED files**: Define genomic regions of interest
- **VCF files**: Contain genetic variant information including multi-allelic sites
- **BAM files**: Store aligned sequencing reads with timing data

## Features

- ✅ **Multi-Allelic Variant Support**: Full detection and analysis of variants with multiple alternate alleles
- ✅ **Enhanced Variant Grouping**: Preserves multi-allelic variant relationships with unique group identifiers
- ✅ **Multi-format Integration**: Seamlessly combines BED, VCF, and BAM data
- ✅ **Timestamp Extraction**: Advanced parsing of nanopore timing data from BAM auxiliary tags
- ✅ **Variant Overlap Detection**: Identifies which reads contain specific variants and which alternate allele
- ✅ **Comprehensive Output**: TSV format with detailed variant and timing information plus grouping metadata
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
| allele_match | Which allele the read contains | VARIANT:ALT1:G |
| variant_chrom | Variant chromosome | chr1 |
| variant_pos | Variant position (1-based) | 1500 |
| variant_ref | Reference allele | A |
| variant_alt | Alternative allele(s) | G,T |
| variant_description | Human-readable variant | A>G |
| variant_type | Type of genetic variant | SNV |
| region | Genomic coordinates | chr1:1000-2000 |
| region_name | BED region identifier | gene_region_1 |
| read_start | Read alignment start | 1450 |
| read_end | Read alignment end | 1550 |
| mapping_quality | Read mapping quality | 60 |
| num_alts | Number of alternate alleles | 2 |

### Example Output

#### Single Allele Variant
```tsv
read_id	timestamp	allele_match	variant_chrom	variant_pos	variant_ref	variant_alt	variant_description	variant_type	region	region_name	read_start	read_end	mapping_quality	num_alts
read001	2021-02-15T14:30:05.123Z	VARIANT:ALT1:G	chr1	1500	A	G	A>G	SNV	chr1:1000-2000	region1	1450	1550	60	1
read002	1613397005.123	REFERENCE	chr1	1500	A	G	A>G	SNV	chr1:1000-2000	region1	1480	1580	55	1
```

#### Multi-Allelic Variant
```tsv
read_id	timestamp	allele_match	variant_chrom	variant_pos	variant_ref	variant_alt	variant_description	variant_type	region	region_name	read_start	read_end	mapping_quality	num_alts
read003	2021-02-15T14:30:06.456Z	VARIANT:ALT1:T	chr2	2000	C	T	C>T	SNV	chr2:1900-2100	region2	1950	2050	62	3
read004	2021-02-15T14:30:07.789Z	VARIANT:ALT2:G	chr2	2000	C	G	C>G	SNV	chr2:1900-2100	region2	1960	2060	58	3
read005	2021-02-15T14:30:08.012Z	VARIANT:ALT3:A	chr2	2000	C	A	C>A	SNV	chr2:1900-2100	region2	1970	2070	61	3
read006	2021-02-15T14:30:09.345Z	REFERENCE	chr2	2000	C	T,G,A	C>T; C>G; C>A	SNV,SNV,SNV	chr2:1900-2100	region2	1980	2080	59	3
```

### Multi-Allelic Support & Grouping

VarClock fully supports multi-allelic variants (positions with multiple alternate alleles) with enhanced grouping features:
- The `allele_match` column indicates which specific alternate allele was detected (e.g., `VARIANT:ALT1:T` for the first alternate)
- The `num_alts` column shows the total number of alternate alleles at that position
- **NEW**: `variant_group_id` provides unique identifiers for variant groups (e.g., `chr1:10030:ATG>A,ATGTG,AT`)
- **NEW**: `variant_summary` offers human-readable descriptions (e.g., `chr1:10030 3-allelic ATG>A,ATGTG,AT`)
- For multi-allelic sites, `variant_alt` may contain comma-separated alleles when showing reference matches
- Each read is analyzed against all possible alleles and reports the best match

#### Enhanced Output Format
The new output includes grouping columns that preserve multi-allelic variant relationships:

```tsv
read_id	timestamp	allele_match	variant_chrom	variant_pos	variant_ref	variant_alt	variant_description	variant_type	variant_group_id	variant_summary	region	region_name	read_start	read_end	mapping_quality	num_alts
read001	2021-02-15T14:30:05.123Z	VARIANT:ALT1:T	chr1	1000	C	T	C>T	SNV	chr1:1000:C>T	chr1:1000 C>T	chr1:900-1100	region1	950	1050	60	1
read002	2021-02-15T14:30:06.456Z	OTHER:N	chr2	2000	ATG	A,ATGTG,AT	ATG>A; ATG>ATGTG; ATG>AT	DEL,INS,DEL	chr2:2000:ATG>A,ATGTG,AT	chr2:2000 3-allelic ATG>A,ATGTG,AT	chr2:1900-2100	region2	1950	2050	58	3
```

#### Grouping Analysis Examples

```bash
# Count reads per variant group
zcat results.tsv.gz | tail -n +2 | cut -f10,11 | sort | uniq -c

# Filter multi-allelic variants only
zcat results.tsv.gz | awk -F'\t' 'NR==1 || $NF > 1' > multiallelic_only.tsv

# Analyze specific variant group
zcat results.tsv.gz | grep "chr1:10030:ATG>A,ATGTG,AT" | cut -f1,3
```

Use the provided analysis script for comprehensive grouping analysis:
```bash
python examples/analyze_multiallelic_groups.py results.tsv.gz analysis_output
```

See [MULTIALLELIC_GROUPING.md](MULTIALLELIC_GROUPING.md) for detailed documentation on variant grouping features.

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