# vcf-toolkit

![CI](https://github.com/NIBIL401/vcf-toolkit/actions/workflows/ci.yml/badge.svg)

A small, **dependency-free** Python toolkit for exploring and filtering variant calls
in [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format, the standard output
of NGS variant callers such as GATK HaplotypeCaller, bcftools and DeepVariant.

## Features

- **Streaming parser** for `.vcf` / `.vcf.gz` with INFO flags/values, FORMAT/sample fields, multi-allelic and symbolic alleles
- **Variant classification**: SNP / INDEL / MNP / OTHER (symbolic SVs like `<DEL>`)
- **Summary statistics**:
  - variant counts by type and chromosome
  - transition/transversion (**Ti/Tv**) ratio, a key call-set quality metric (~2.0–2.1 for WGS, ~3.0 for exomes)
  - substitution spectrum (e.g. `C>T`, `A>G`)
  - QUAL and depth distributions
  - indel length distribution
  - per-sample genotype counts and **het/hom ratio**
- **Filtering** by QUAL, depth, allele frequency, FILTER=PASS, variant type, genomic region(s) and genotype, with a report of how many variants each criterion removed

## Install

```bash
git clone https://github.com/NIBIL401/vcf-toolkit.git
cd vcf-toolkit
pip install -e ".[test]"
```

## Usage

```bash
vcfkit stats examples/sample.vcf
```

```
Total variants     20
  INDEL            4
  MNP              1
  OTHER            1
  SNP              14
Multi-allelic      1
Ti / Tv            10 / 3  (ratio 3.333)
QUAL               min 9.8  median 80.3  mean 102.66  max 310.4
...
Sample TUMOR       {'het': 16, 'hom_alt': 4}  het/hom=4.0
```

```bash
# High-confidence SNPs only
vcfkit filter examples/sample.vcf --pass-only --type snp --min-qual 30 --min-depth 20 -o hq_snps.vcf

# Variants in a region, heterozygous in the first sample, written gzipped
vcfkit filter calls.vcf.gz --region chr2:25,200,000-25,300,000 --genotype het -o dnmt3a.vcf.gz

# Pipe to stdout
vcfkit filter calls.vcf --type indel | grep -v '^##' | head
```

Run `vcfkit filter -h` for every option. Filter statistics go to stderr:

```
kept 12 / 20 variants
  removed (type): 5
  removed (not_pass): 3
```

### As a library

```python
from vcf_toolkit.vcf import VCFReader
from vcf_toolkit.analysis import VariantFilter, VariantStats

flt, stats = VariantFilter(pass_only=True, types={"SNP"}), VariantStats()
with VCFReader("calls.vcf.gz") as reader:
    for v in flt.apply(reader):
        stats.add(v, reader.samples)
print(stats.summary()["ti_tv_ratio"])
```

## Example data

`examples/sample.vcf` is a small **simulated** tumour/normal call set on GRCh38-style
coordinates. It covers SNPs, indels, an MNP, a multi-allelic site, a symbolic deletion,
failed filters and a missing genotype.

## Tests

```bash
pytest -v
```

## License

MIT
