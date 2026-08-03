# WRATH Codebase Documentation

**WRATH** (WRapped Analysis of Tagged Haplotypes) is a bioinformatics tool for exploring structural variation from haplotagging/linked-read sequencing data. This document provides a detailed description of each program in the codebase and explains how they work together.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Main Orchestrator Script](#main-orchestrator-script)
3. [Barcode Parsing Module](#barcode-parsing-module)
4. [SV Detection Module](#sv-detection-module)
5. [Workflow Diagrams](#workflow-diagrams)
6. [Data Flow](#data-flow)

---

## Architecture Overview

WRATH consists of two main functional modules:

| Module | Language | Purpose |
|--------|----------|---------|
| **Barcode Parsing** | Python | Pre-processing of raw haplotagging data to extract and encode barcodes |
| **SV Detection** | Python + R | Core analysis pipeline for matrix construction, outlier detection, and visualization |

The main entry point is a Bash script (`wrath`) that orchestrates the entire pipeline.

---

## Main Orchestrator Script

### `wrath` (Bash)

**Location:** `wrath` (root directory)

**Author:** Anna Farre Orteu (2022)

**Purpose:** Main CLI entry point that orchestrates the complete WRATH pipeline for structural variant detection from haplotagging data.

#### Command-Line Options

| Flag | Parameter | Description |
|------|-----------|-------------|
| `-g` | FASTAFILE | Reference genome in FASTA format |
| `-c` | CHROMOSOMENAME | Target chromosome/scaffold name |
| `-w` | WINDOWSIZE | Size of genomic windows (default: 50000 bp) |
| `-a` | FILELIST | Text file listing BAM file paths |
| `-t` | THREADS | Number of threads (default: 1) |
| `-p` | - | Skip plotting the heatmap |
| `-v` | - | Verbose mode (for matrix generation) |
| `-x` | STEP | Resume from specific step: `makewindows`, `getbarcodes`, `matrix`, `outliers`, `plot` |
| `-l` | - | Enable automatic SV detection |
| `-s` | START | Start position for subsetting windows |
| `-e` | END | End position for subsetting windows |

#### Pipeline Stages

The script executes the following stages sequentially:

1. **Make Genomic Windows**
   - Uses `samtools faidx` to index the reference genome
   - Extracts chromosome sizes
   - Uses `bedtools makewindows` to create fixed-width genomic windows
   - Output: `wrath_out/beds/windows_{winSize}_{chromosome}_{start}_{end}.bed`

2. **Get Barcodes**
   - Extracts barcodes from BAM files using `samtools view`
   - Filters for mapping quality >= 20
   - Parses BX:Z tags containing molecule barcodes
   - Sorts, compresses, and indexes barcode BED files
   - Output: `wrath_out/beds/barcodes_{chromosome}_{start}_{end}_sorted_{group}.bed.gz`

3. **Generate Similarity Matrix**
   - Calls `jaccard_matrix_simplequeue.py`
   - Computes Jaccard index between all window pairs
   - Output: `wrath_out/matrices/jaccard_matrix_{winSize}_{chromosome}_{start}_{end}_{group}.txt`

4. **Detect Outliers** (if `-l` flag provided)
   - Calls `outlier_detection.R`
   - Fits background model and identifies statistical outliers
   - Output: `wrath_out/outliers/outliers_{winSize}_{chromosome}_{start}_{end}_{group}.csv`

5. **Generate Visualization**
   - Calls `plot_heatmap.py` (basic mode) or `sv_detection_and_heatmap.py` (with SV detection)
   - Output: `wrath_out/plots/heatmap_{winSize}_{chromosome}_{start}_{end}_{group}.png`

---

## Barcode Parsing Module

Located in `barcode_parsing/`

### `parse_haptag_barcodes.py`

**Author:** Simon Martin (2022)

**Purpose:** Parses haplotagging barcodes from Illumina index reads and adds them to sequencing read headers.

#### How It Works

1. **Input Processing**
   - Reads paired FASTQ files: read1, read2, index1, index2
   - Index reads are 13 bp: `CCCCCCNAAAAAA` (I1) and `DDDDDDNBBBBBB` (I2)

2. **Barcode Extraction**
   - Extracts 4 barcode components (A, B, C, D) from index reads
   - A and C from index1, B and D from index2
   - Each barcode is 6 bp

3. **Barcode Matching**
   - Looks up barcode sequences in reference files (`BC_A.txt`, `BC_B.txt`, etc.)
   - Supports exact matching or single-mismatch tolerance
   - Generates correction dictionaries for error-tolerant matching

4. **Output Generation**
   - Constructs BX tag: concatenation of A+C+B+D codes
   - Adds BX:Z, RX:Z (raw sequences), and QX:Z (quality scores) to read headers
   - Optionally demultiplexes by sample using C barcode
   - Uses multiprocessing queues for parallel file writing

#### Key Classes

- **`Read`**: Simple container for FASTQ read data (name, sequence, quality)
- **`missing_dict`**: Dictionary returning a default "missing" value for unknown keys
- **`mirror_dict`**: Dictionary returning the key itself if not found

#### Command-Line Usage

```bash
python parse_haptag_barcodes.py \
  -R read1.fq.gz read2.fq.gz \
  -I index1.fq.gz index2.fq.gz \
  --output_label experiment_name \
  --demult_file demult.txt \
  --output_dir /output/path \
  --count_barcodes
```

---

### `make_demult_file.py`

**Author:** Simon Martin (2022)

**Purpose:** Creates a demultiplexing file that maps C barcodes to sample names.

#### How It Works

1. Reads C barcode reference file (`BC_C.txt`) mapping codes to sequences
2. Reads sample barcode file mapping sequences to sample names
3. Creates lookup table: C-code -> sample name
4. Outputs demultiplexing file for `parse_haptag_barcodes.py`

#### Command-Line Usage

```bash
python make_demult_file.py \
  --C_barcode_file BC_C.txt \
  --sample_barcode_file sample_barcodes.txt \
  -o demult_file.txt
```

---

## SV Detection Module

Located in `sv_detection/`

### `jaccard_matrix_simplequeue.py`

**Author:** Anna Orteu (2023)

**Purpose:** Computes a pairwise Jaccard similarity matrix for barcode sharing between genomic windows.

#### How It Works

1. **Input Parsing**
   - Reads window BED file defining genomic intervals
   - Opens tabix-indexed barcode BED file using `pysam`

2. **Matrix Construction**
   - For each window pair (i, j) where j >= i:
     - Fetches barcodes in window i
     - Fetches barcodes in window j
     - Computes intersection and union of barcode sets
     - Calculates Jaccard index: `J(A,B) = |A ∩ B| / |A ∪ B|`

3. **Parallel Processing Architecture**
   - **Worker Processes** (`freqs_wrapper`): Compute Jaccard values for assigned windows
   - **Sorter Thread** (`sorter`): Maintains order of results despite parallel completion
   - **Writer Thread** (`writer`): Writes sorted results to output file
   - **Checker Thread** (`checkStats`): Reports progress statistics every 10 seconds

4. **Output**
   - CSV matrix where entry (i,j) is the Jaccard index between windows i and j
   - Only upper triangle is computed (symmetric matrix)

#### Command-Line Usage

```bash
python jaccard_matrix_simplequeue.py \
  -w windows.bed \
  -b barcodes.bed.gz \
  -o jaccard_matrix.txt \
  -t 10 \
  --verbose
```

---

### `outlier_detection.R`

**Author:** Anna Orteu

**Purpose:** Identifies statistical outliers in the barcode-sharing matrix using Z-scores and prediction bands from a fitted decay model.

#### Algorithm

1. **Data Preparation**
   - Reads Jaccard matrix
   - Sets lower triangle and diagonal to NA (analyze upper triangle only)
   - Reshapes to long format with row, column, value, and distance-from-diagonal

2. **Z-Score Calculation**
   - Groups values by distance from diagonal
   - Computes Z-score for each value within its distance group
   - Default threshold: |Z| > 2

3. **Model Fitting**
   - Fits double-exponential decay model: `y ~ exp(a + b * exp(-x * c))`
   - Models expected barcode sharing as function of genomic distance
   - Uses `nls()` for non-linear least squares fitting

4. **Prediction Bands**
   - Computes 95% prediction bands using `nlraa::predict2_nls()`
   - Points outside bands are flagged

5. **Outlier Definition**
   - Outliers must satisfy BOTH criteria:
     - |Z-score| > threshold (default: 2)
     - Value outside 95% prediction bands

6. **Output**
   - Scatter plot with fitted model and prediction bands
   - CSV file listing outlier window pairs with statistics

#### Command-Line Usage

```bash
Rscript outlier_detection.R input_matrix.txt output_prefix
```

---

### `plot_heatmap.py`

**Author:** Anna Orteu (2023)

**Purpose:** Generates a heatmap visualization of the Jaccard similarity matrix.

#### How It Works

1. Reads Jaccard matrix and window coordinates
2. Renames axes using genomic positions
3. Log-transforms data: `log(value + 0.0001) * 100`
4. Generates heatmap using seaborn with YlGnBu colormap
5. Configures tick labels based on genomic coordinates
6. Saves high-resolution PNG (30x30 inches)

#### Command-Line Usage

```bash
python plot_heatmap.py \
  -m jaccard_matrix.txt \
  -w windows.bed \
  -o heatmap.png
```

---

### `plot_2matrices_together.py`

**Author:** Anna Orteu (2023)

**Purpose:** Creates a combined heatmap showing two populations' matrices in upper and lower triangles for visual comparison.

#### How It Works

1. Reads two Jaccard matrices and window coordinates
2. Transposes matrix 1 and adds to matrix 2
3. Creates symmetric visualization:
   - Upper triangle: Population 1
   - Lower triangle: Population 2
4. Applies same styling as single-matrix heatmap

#### Use Case

Comparing barcode-sharing patterns between populations that may differ in structural variants (e.g., case vs. control, different morphs).

#### Command-Line Usage

```bash
python plot_2matrices_together.py \
  -m1 matrix_pop1.txt \
  -m2 matrix_pop2.txt \
  -w windows.bed \
  -o comparison_heatmap.png
```

---

### `sv_detection.py`

**Author:** Anna Orteu (2023)

**Purpose:** Clusters detected outliers into candidate structural variant calls.

#### Algorithm

1. **Outlier Clustering**
   - Uses scikit-learn's `AgglomerativeClustering`
   - Single linkage with distance threshold of 3 windows
   - Groups nearby outlier points into SV candidates

2. **Breakpoint Identification**
   - For each cluster, finds:
     - Minimum/maximum column indices
     - Minimum/maximum row indices
   - These define the genomic extent of the candidate SV

3. **SV Output**
   - Calculates SV length: `(maxcol - minrow) * window_size`
   - Sorts by length (descending)
   - Outputs: SV_id, start, end, length

#### Command-Line Usage

```bash
python sv_detection.py \
  -m jaccard_matrix.txt \
  -o outliers.csv \
  -s sv_candidates.txt \
  -f 50000
```

---

### `sv_detection_and_heatmap.py`

**Author:** Anna Orteu (2023)

**Purpose:** Combined script that performs SV detection and generates an annotated heatmap in a single step.

#### How It Works

1. Performs same clustering as `sv_detection.py`
2. Generates heatmap with outlier clusters overlaid:
   - Black dots: cluster start positions (minrow, mincol)
   - Magenta dots: cluster end positions (maxrow+1, maxcol+1)
3. Outputs both visualization and SV candidate table

#### Output

- Heatmap PNG with SV candidates marked
- CSV with columns: SV_id, chromosome, start, end, length

#### Command-Line Usage

```bash
python sv_detection_and_heatmap.py \
  -m jaccard_matrix.txt \
  -o outliers.csv \
  -p heatmap.png \
  -s sv_candidates.txt \
  -f 50000 \
  -c chr1 \
  -w windows.bed
```

---

## Workflow Diagrams

### Complete Pipeline Flow

```
                    RAW HAPLOTAGGING DATA
                            │
                            ▼
    ┌───────────────────────────────────────────┐
    │         BARCODE PARSING MODULE            │
    │                                           │
    │  index1.fq.gz ─┬─► parse_haptag_barcodes.py
    │  index2.fq.gz ─┤         │
    │  read1.fq.gz  ─┤         ▼
    │  read2.fq.gz  ─┘   Barcoded FASTQ files
    │                         │
    │  (demult_file.py)       ▼
    │                    Per-sample FASTQs
    └───────────────────────────────────────────┘
                            │
                    [External: BWA alignment]
                            │
                            ▼
                      BAM FILES (with BX tags)
                            │
                            ▼
    ┌───────────────────────────────────────────┐
    │            WRATH MAIN PIPELINE            │
    │                                           │
    │  1. makewindows (bedtools)                │
    │         │                                 │
    │         ▼                                 │
    │  2. getbarcodes (samtools)                │
    │         │                                 │
    │         ▼                                 │
    │  3. jaccard_matrix_simplequeue.py         │
    │         │                                 │
    │         ▼                                 │
    │  4. outlier_detection.R  ◄── optional     │
    │         │                                 │
    │         ▼                                 │
    │  5. sv_detection_and_heatmap.py           │
    │     or plot_heatmap.py                    │
    └───────────────────────────────────────────┘
                            │
                            ▼
                    ┌───────────────┐
                    │   OUTPUTS     │
                    │ - Heatmaps    │
                    │ - SV table    │
                    │ - Outliers    │
                    │ - Matrices    │
                    └───────────────┘
```

### Python-R Integration

```
jaccard_matrix_simplequeue.py
        │
        │ (CSV matrix)
        ▼
  outlier_detection.R
        │
        │ (CSV outliers + PNG plot)
        ▼
sv_detection_and_heatmap.py
        │
        ▼
  Final outputs
```

---

## Data Flow

### File Dependencies

```
Reference Genome (.fa)
        │
        ├──► size.genome
        │
BAM files ──► barcodes_{chr}_{region}_sorted_{group}.bed.gz
        │
        ├──► windows_{winSize}_{chr}_{region}.bed
        │           │
        │           ▼
        └──► jaccard_matrix_{params}.txt
                    │
                    ├──► outliers_{params}.csv
                    │           │
                    │           ▼
                    └──► sv_{params}.txt
                                │
                                ▼
                        heatmap_{params}.png
```

### Key Data Structures

| File Type | Format | Description |
|-----------|--------|-------------|
| Window BED | TSV | `chrom  start  end` |
| Barcode BED | TSV (gzipped) | `chrom  pos  pos  BX:Z:barcode` |
| Jaccard Matrix | CSV | N x N symmetric matrix of Jaccard values |
| Outliers | CSV | `nrow,ncol,value,Estimate,Est.Error,Q2.5,Q97.5,upper,lower,z_score` |
| SV Candidates | CSV | `SV_id,chromosome,start,end,length` |

---

## Dependencies Summary

### Python Packages

| Package | Used By | Purpose |
|---------|---------|---------|
| `numpy` | All Python scripts | Array operations |
| `pandas` | All Python scripts | DataFrame handling |
| `pysam` | jaccard_matrix | Tabix file access |
| `seaborn` | Plotting scripts | Heatmap visualization |
| `matplotlib` | Plotting scripts | Figure rendering |
| `scikit-learn` | SV detection | Agglomerative clustering |

### R Packages

| Package | Used By | Purpose |
|---------|---------|---------|
| `ggplot2` | outlier_detection.R | Plotting |
| `tidyr` | outlier_detection.R | Data reshaping |
| `dplyr` | outlier_detection.R | Data manipulation |
| `nlraa` | outlier_detection.R | Prediction bands |

### Command-Line Tools

| Tool | Used By | Purpose |
|------|---------|---------|
| `samtools` | wrath script | BAM processing, FASTA indexing |
| `bedtools` | wrath script | Window generation |
| `tabix/bgzip` | wrath script | Barcode file indexing |

---

## Example Run

### SLURM Array Job

The example in `example_run/example_wrath_slurm_array.sh` demonstrates running WRATH across multiple chromosomes in parallel:

```bash
#!/bin/bash
#SBATCH --array=1-195

chr=$(sed -n "$SLURM_ARRAY_TASK_ID"p Hera_chr)
winSize=50000
threads=10
start=700000
end=850000
group=malleti.txt

wrath -g ~/genomes/Hmel/Hmel2.5.fa \
      -c ${chr} \
      -w ${winSize} \
      -a ${group} \
      -t ${threads} \
      -l \
      -s ${start} \
      -e ${end}
```

---

*Documentation generated March 2026*
