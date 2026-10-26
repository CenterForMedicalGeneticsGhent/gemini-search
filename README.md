# gemini-search

A minimal, pure-Python library for querying and searching GEMINI genetic variation databases (SQLite). Compatible with Python 3.14+.

## Overview

`gemini-search` extracts the core search and query capabilities from the original GEMINI framework into a lightweight, installable Python library. All legacy CLI commands, loading tools, and web UI components have been removed.

Features:
- SQL query execution against GEMINI databases with automatic decompression of genotype BLOBs.
- Named sample genotype access (`gts.SAMPLE`, `gt_types.SAMPLE`, `gt_depths.SAMPLE`).
- Genotype filters (`gt_filter = "gt_types.sample1 == HET and gt_types.sample2 == HOM_REF"`).
- Wildcard filters based on sample pedigree attributes (e.g. `(gt_types).(phenotype == 2).(==HET).(all)`).
- Genomic region filtering (`region = "chr1:100-200"`).
- VCF INFO tag dictionary unpacking (`row.info`).
- Flexible output formatting (`DefaultRowFormat`, `JSONRowFormat`, `VCFRowFormat`, `CarrierSummary`, `TPEDRowFormat`, `SampleDetailRowFormat`).

## Installation

### Using Pixi
```bash
pixi install
```

### Using Pip
```bash
pip install .
```

## Quick Start

```python
from gemini import GeminiQuery, HET, HOM_REF

gq = GeminiQuery("my_gemini.db")

# Basic query
gq.run("SELECT chrom, start, end, ref, alt, gene, info FROM variants WHERE gene = 'SCNN1D'")
for row in gq:
    print(row["chrom"], row["start"], row["gene"], row.info.get("AF"))

# Query with genotype columns and filters
gq.run(
    "SELECT chrom, start, gts.sample1, gt_types.sample2 FROM variants",
    gt_filter="gt_types.sample1 == HET and gt_types.sample2 == HOM_REF",
    region="chr1:1000000-2000000",
)
for row in gq:
    print(row["chrom"], row["start"], row["gts.sample1"])
```

## Development

Run tests:
```bash
pixi run test
```

Run linter:
```bash
pixi run lint
```

## License

MIT License (see [LICENSE](LICENSE)).
