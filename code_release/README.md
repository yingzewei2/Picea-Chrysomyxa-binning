# Host/pathogen transcript binning and host differential expression for *Picea schrenkiana* infected by *Chrysomyxa deformans*

Code accompanying:

> Ying, Z.; Hao, J.; Chen, H.; Shi, Y.; Mamuti, M.; Wang, W.; Zhu, X. Transcriptome and metabolome analysis reveals downregulation of cuticle, photosynthesis and diterpene resin acid pathways in *Picea schrenkiana* needles colonised by *Chrysomyxa deformans*. *Journal of Fungi* (submitted).

Raw RNA-seq reads: NCBI SRA, BioProject [PRJNA1547998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1547998).
Raw LC–MS/MS data: OMIX, BioProject PRJCA077544.

## Scripts

| Script | Purpose | Manuscript |
|---|---|---|
| `scripts/01_binning_host_pathogen.py` | Assigns each unigene of the de novo assembly to host or pathogen and writes separate count matrices | Section 2.5; Supplementary Methods S1; Table S1 |
| `scripts/02_host_DE_KEGG.py` | Re-normalises host counts within the host subset, runs DESeq2 with a paired design and performs KEGG enrichment against a host background | Section 2.6; Tables S3, S5 |

`run_01_binning.bat` and `run_02_host_DE_KEGG.bat` are Windows launchers that run the scripts from the folder they are in.

### Binning rules (`01_binning_host_pathogen.py`)

Applied in a fixed order:

- **A (expression):** a unigene with zero counts in all asymptomatic (H) samples is assigned to the pathogen; otherwise to the host.
- **B (taxonomy):** the genus of the NR best hit is matched against curated plant and fungal genus lists; B overrides A whenever it applies, and discordance is recorded.
- **C (GC content, unigenes without an NR hit only):** GC < 42.0% supports host and GC > 52.0% supports pathogen.

Confidence: high = B available; medium = A and C agree; low = A only.

### Differential expression (`02_host_DE_KEGG.py`)

- Model: `~ tree + stage` (PyDESeq2), using the paired within-tree design.
- Threshold: Benjamini–Hochberg q < 0.05 and |log2FC| ≥ 1.
- KEGG enrichment: hypergeometric test with host unigenes carrying a KO as background, BH-adjusted.

## Inputs

Both scripts read the files delivered by the sequencing provider (OE Biotech), placed in the script folder:

| File | Used by |
|---|---|
| `counts_anno.xls` (tab-delimited counts + annotation table) | 01, 02 |
| `Unigene.fa` (de novo assembly, 59,665 unigenes) | 01 (GC content); 02 (optional, FPKM) |
| `host_counts_integer_for_DESeq2.tsv` (output of 01) | 02 |
| `Unigene.KEGG.gene.anno.xls`, `Unigene.KEGG.classification.xls` | 02 |

Sample names follow the SRA submission: H-1…H-5, S1-1…S1-5 and S2-1…S2-5 (trees T1–T5). Samples H-excl, S1-excl and S2-excl (one tree excluded because its S1 sample failed sequencing QC) were used for the assembly only and are not part of the count tables analysed here.

## Requirements

Python 3.10 with pandas 2.3.3, numpy 2.2.6, scipy 1.15.3 and pydeseq2 0.5.2 (versions used for the manuscript).

```
pip install pandas==2.3.3 numpy==2.2.6 scipy==1.15.3 pydeseq2==0.5.2
```

## Usage

```
python scripts/01_binning_host_pathogen.py      # writes ./分箱结果/
python scripts/02_host_DE_KEGG.py               # writes ./2.gene/
```

Console messages and some output file names are in Chinese.

## License

MIT (see `LICENSE`).
