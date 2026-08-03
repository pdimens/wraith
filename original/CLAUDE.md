# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**WRATH** (WRapped Analysis of Tagged Haplotypes) is a bioinformatics tool for exploring structural variation from haplotagging / linked-read sequencing data.

In practice, WRATH takes aligned reads carrying molecule barcodes, partitions chromosomes or scaffolds into genomic windows, quantifies barcode sharing between windows, builds barcode-sharing matrices, renders heatmaps, and identifies candidate structural-variant signals as statistical outliers relative to the expected distance-decay background.

The project focuses on:
- visualization of barcode-sharing structure from linked-read data,
- prioritization of candidate structural variants (SVs) for follow-up,
- robust handling of BAM-based population-scale haplotagging datasets,
- reproducible CLI-driven analyses suitable for genomics workflows,
- scientifically correct extension of the current methods without obscuring their assumptions and limitations.

WRATH is not just a generic SV caller. It is a linked-read exploration and prioritization system whose value comes from the combination of:
- barcode-aware matrix construction,
- heatmap-based visual interpretation,
- outlier detection against a fitted background model,
- targeted re-analysis of loci of interest,
- comparison of candidate regions across populations.

# WRATH – Technical Specification

## 1. Overview

**WRATH** is a predominantly Python-based structural-variation analysis toolkit with some R-based plotting / downstream analysis support.

The primary goals are:
- parse or consume haplotagged linked-read alignments that preserve molecule barcode information,
- compute windowed barcode-sharing summaries across chromosomes or user-specified intervals,
- generate interpretable heatmaps and related plots for manual SV exploration,
- identify putative SV-associated outlier window pairs using a statistical background model,
- support both genome-wide scans and fine-scale local re-analysis,
- keep the method transparent about where signals are strong, weak, or ambiguous,
- make it straightforward to extend the existing codebase for new signal types, better QC, and improved workflow integration.

This repository should be treated as scientific software. Changes must preserve methodological traceability, reproducibility, and interpretability.

---

## 2. Scientific / Algorithmic Context

WRATH operates on the core observation that reads derived from the same long DNA molecule share a barcode and, absent rearrangement, should map near one another in the reference genome.

At a high level, WRATH works as follows:
1. Use aligned reads with barcode / molecule tags.
2. Split each chromosome or scaffold into genomic windows.
3. Collect the set of barcodes observed in each window.
4. Compute pairwise barcode-sharing scores between windows using the Jaccard index.
5. Store those values in an `n x n` matrix.
6. Visualize the matrix as a heatmap.
7. Model expected barcode sharing as a function of distance from the diagonal.
8. Flag outlier window pairs as candidate SV signals.

Important biological / technical implications:
- Windows near one another are expected to share more barcodes than distant windows.
- Large inversions, duplications, and translocation-like events can create excess barcode sharing away from the diagonal.
- Deletions can appear as depleted regions and, when large enough, may also produce flanking signatures.
- Some heatmap patterns are non-unique: multiple biological events can produce similar visual signatures.
- WRATH prioritizes candidate SVs; it does **not** fully classify all events automatically.
- Very small SVs relative to molecule length are often not detectable with this approach.

This makes WRATH closer to a scientific exploration platform than a fully automated black-box caller.

---

## 3. Target Runtime Environment

* **Primary OS**: Linux (preferred for bioinformatics workflows and large BAM processing).
* **Also supported**: macOS and WSL2 on Windows when toolchain dependencies are available.
* **Language**: Python 3.11+.
* **Secondary language**: R for legacy plotting / downstream analysis where needed.
* **Dependency Management**: `uv` for Python; use standard project-local dependency management for R if R code remains in scope.
* **Interface style**: command-line first. Any GUI / notebook layer is secondary.
* **Containerization**: Docker is encouraged for reproducibility, especially for multi-user or server deployment.

**CRITICAL**: This project uses `uv` exclusively for Python package management.

**Always use:**
```bash
uv add package_name
uv sync
uv run pytest
uv run python -m wrath.cli
```

**Never use:**
```bash
pip install xxx
python script.py
python -m pytest
```

For R code:
- keep R dependencies explicit,
- isolate R-facing functionality behind clear interfaces,
- do not silently move core scientific logic into R if it weakens testability or portability.

---

## 4. Core Concepts

### 4.1 Haplotagging / Linked-Read Data
Reads originating from the same long molecule carry a shared barcode. WRATH assumes barcode information is present and trustworthy enough to support molecule-level inference.

### 4.2 BX / Molecule Tags
Barcode / molecule identity must be accessible from aligned reads. If preprocessing is required to convert raw index information into BAM-friendly tags, that step must be explicit and testable.

### 4.3 Genomic Windows
Chromosomes / scaffolds are partitioned into fixed-width windows. Window size is a scientifically meaningful parameter that affects sensitivity, precision, runtime, and noise.

### 4.4 Barcode Sets per Window
Each window is represented by the set of barcodes observed among reads mapping into that window.

### 4.5 Barcode-Sharing Matrix
For each pair of windows, WRATH computes a similarity score, currently the Jaccard index:

```text
J(A, B) = |A ∩ B| / |A ∪ B|
```

### 4.6 Background Distance-Decay Model
Expected barcode sharing decays with genomic distance from the diagonal. Candidate signals are defined relative to that background, not by raw sharing alone.

### 4.7 Candidate SVs
Outliers are prioritization targets, not final truth calls. Any workflow or UI must keep this distinction explicit.

### 4.8 Manual Interpretation
Heatmaps are a first-class output. Visual inspection is part of the intended analysis model, not an afterthought.

### 4.9 Methodological Limits
Molecule length and window size impose hard detection limits. Very small SVs or ambiguous repetitive regions may require orthogonal methods.

---

## 5. Functional Requirements

### 5.1 Input Data
WRATH should support the following primary inputs:
- coordinate-sorted BAM files,
- indexed reference genome coordinates / scaffold definitions,
- reads carrying barcode / molecule information,
- user-specified chromosome, scaffold, or genomic interval selection,
- optional groupings of samples by population / phenotype / condition.

Optional supporting inputs may include:
- BED files for intervals of interest,
- chromosome / scaffold allow-lists,
- sample manifests,
- previously generated candidate-SV BED-like outputs for overlap analysis,
- metadata for comparative analyses.

### 5.2 Input Validation
Validation must run before heavy computation.

#### 5.2.1 File-Level Validation
- BAM exists and is readable.
- BAM index exists when required.
- Reference coordinate definitions are available.
- Requested contigs / intervals exist.
- Output directory is writable.

#### 5.2.2 Content-Level Validation
- Reads contain the expected barcode / molecule tags or equivalent fields.
- Mapping quality thresholds are valid.
- Window size is positive and biologically sensible.
- Requested interval size is compatible with selected window size.
- Thread count is valid.

#### 5.2.3 Scientific Validation
- Warn when window size is too large relative to the region of interest.
- Warn when window size is too small for available coverage.
- Warn when the requested analysis targets SV sizes below the linked-read method's likely detection range.
- Warn when low mapping quality or sparse barcode density will likely make the matrix uninterpretable.

Validation behavior:
- Blocking issues stop execution.
- Scientific-quality warnings are shown prominently but do not always block execution.
- Error messages must tell the user exactly what parameter, file, contig, or interval is invalid.

### 5.3 Core Analysis Workflow
The main workflow should remain decomposed into explicit stages:
1. Parse / read aligned barcode-tagged reads.
2. Build per-window barcode sets.
3. Compute pairwise window similarity matrix.
4. Fit background distance-decay model.
5. Compute outlier statistics and apply filters.
6. Generate plots and candidate-SV outputs.
7. Optionally compare outputs across groups.

Each stage should be callable independently where practical.

### 5.4 Barcode Parsing
Support or preserve logic for:
- extracting molecule / barcode identifiers from alignments,
- filtering low-quality reads,
- handling duplicate or malformed barcode tags,
- optionally restricting to target intervals before building barcode sets.

### 5.5 Matrix Construction
- Compute barcode-sharing values per chromosome / scaffold / interval.
- Current default similarity metric is Jaccard index.
- Matrix computation must support chunking or streaming if memory becomes limiting.
- Symmetry should be preserved and tested.
- Diagonal handling should be explicit.

### 5.6 Statistical Outlier Detection
WRATH currently assumes barcode sharing decays with distance from the diagonal and fits a double-exponential background model.

Required behaviors:
- fit the background model reproducibly,
- expose parameters and thresholds clearly,
- identify points outside prediction bands as candidate outliers,
- support configurable alpha / significance settings,
- support configurable Z-score filtering near the diagonal,
- keep the distinction between “candidate signal” and “called SV” explicit.

### 5.7 Plotting and Visualization
Plotting is a core feature, not optional polish.

Required outputs:
- heatmaps of barcode-sharing matrices,
- optional overlays of genes / annotations / loci of interest,
- model-fit plots for outlier detection,
- candidate-SV summaries,
- interval-focused re-plots at finer resolution.

Plots must be reproducible, scriptable, and savable without interactive sessions.

### 5.8 Comparative Analyses
The codebase should support separate analysis of different sample groups and easy downstream comparison.

Examples:
- population A vs population B,
- phenotype-defined subsets,
- case vs control,
- per-sample exploratory runs.

Comparison outputs should be easy to export in BED-like coordinate formats for downstream overlap analysis.

### 5.9 Outputs
#### 5.9.1 Primary Outputs
- heatmap images,
- candidate-SV tables with genomic coordinates,
- model-fit / outlier plots,
- logs containing key parameter settings,
- optional intermediate matrices or serialized summaries.

#### 5.9.2 Recommended Output Formats
- `.png` / `.pdf` for plots,
- `.tsv` / `.bed` / `.bedpe`-like text outputs for candidate regions,
- optional compressed binary formats for cached matrices,
- structured logs in plain text.

#### 5.9.3 Naming
Output files should encode:
- sample or cohort name,
- chromosome / interval,
- window size,
- filter settings,
- timestamp or run identifier.

Example:
```text
wrath_Hmisippus_Mminus_chr26_6731000_6744000_win100_model_outliers.tsv
```

### 5.10 Performance
- Parallel execution must remain supported.
- Thread count should be configurable.
- Expensive operations should be profiled before major refactors.
- Avoid premature optimization that obscures scientific logic.

### 5.11 No Silent Method Drift
Any change that alters:
- barcode parsing behavior,
- similarity calculation,
- model fitting,
- outlier thresholds,
- filtering logic,
- window semantics,
- output interpretation,

must be treated as a method change and documented clearly.

---

## 6. Architecture & Code Organization

### 6.1 Project Structure

Example layout:

```text
wrath-plus/
  pyproject.toml
  README.md
  CLAUDE.md
  src/
    wrath/
      __init__.py
      cli.py                 # command-line entrypoints
      config.py              # defaults for window size, mapq, thresholds
      barcode.py             # barcode extraction / tag parsing logic
      intervals.py           # chromosome / interval handling
      windows.py             # genomic window generation
      matrix.py              # Jaccard matrix construction
      model.py               # background fitting and outlier logic
      candidates.py          # candidate-SV extraction and export
      plotting.py            # heatmaps and model-fit plots
      io.py                  # BAM / TSV / BED / config I/O
      compare.py             # overlap and group comparison helpers
      qc.py                  # coverage / barcode-density / diagnostic checks
      utils.py               # small shared helpers only
  r/
    plots.R                 # legacy or optional R plotting support
    helpers.R
  tests/
    test_barcode.py
    test_windows.py
    test_matrix.py
    test_model.py
    test_candidates.py
    test_compare.py
  docs/
    method.md
    examples/
  docker/
    Dockerfile
  .gitignore
```

### 6.2 Data Models
Prefer simple typed structures unless classes materially improve clarity.

Potential models:
- `Window`
  - contig, start, end, index
- `WindowBarcodeSet`
  - window, barcodes, read_count, barcode_count
- `MatrixBuildResult`
  - contig, window_size, matrix, window_metadata
- `DecayModelFit`
  - parameters, prediction_band, fit_diagnostics
- `CandidateSV`
  - contig_a, start_a, end_a, contig_b, start_b, end_b, score, z_score, p_value, metadata
- `RunConfig`
  - bam_paths, interval, window_size, mapq_threshold, threads, output_dir, mode

### 6.3 Modules

* `barcode.py`
  - extract / normalize barcode tags
  - validate barcode presence

* `windows.py`
  - generate genomic windows
  - map reads to windows

* `matrix.py`
  - construct per-window barcode sets
  - compute Jaccard matrix

* `model.py`
  - fit distance-decay background model
  - compute prediction bands
  - compute Z-scores / outlier flags

* `candidates.py`
  - convert outliers into exported candidate records
  - merge / filter nearby signals when appropriate

* `plotting.py`
  - render heatmaps
  - render model-fit plots
  - render region-specific summaries

* `compare.py`
  - compare group-specific candidate sets
  - export coordinate overlaps for downstream tools

* `cli.py`
  - provide subcommands such as:
    - `wrath scan`
    - `wrath plot`
    - `wrath candidates`
    - `wrath compare`

---

## 7. Non-Functional Requirements

### 7.1 Robustness
- Must handle large BAMs and long scaffolds without brittle failures.
- Must fail clearly when barcode tags are absent or malformed.
- Must not silently drop contigs, reads, or windows without logging.

### 7.2 Reproducibility
- Every run should record its full parameterization.
- Version strings and git commit hashes should be included when possible.
- Method-changing defaults must be documented.

### 7.3 Testability
- Core mathematical and transformation logic must be unit tested.
- BAM-heavy and plotting-heavy workflows should have integration tests with small fixtures.
- R-dependent functionality must be testable in isolation or treated as optional.

### 7.4 Interpretability
- Heatmaps and candidate tables should remain human-auditable.
- The code must preserve the conceptual mapping between mathematical outputs and biological interpretation.

### 7.5 Performance
- Typical targeted-locus analyses should complete quickly.
- Whole-genome scans should scale with threads and avoid unnecessary repeated BAM passes.

### 7.6 Backward Compatibility
- Preserve existing CLI behavior where feasible.
- If breaking output schemas or flags, document them explicitly and add migration notes.

---

## 8. Testing & QA Plan

### 8.1 Unit Tests
Test pure logic thoroughly:
- barcode normalization,
- window generation,
- Jaccard calculation,
- matrix symmetry,
- decay-model input preparation,
- outlier filtering,
- coordinate export formatting.

### 8.2 Integration Tests
Use small fixture BAMs / synthetic window-barcode summaries to test:
- end-to-end scan on a toy contig,
- interval-restricted plotting,
- candidate table generation,
- comparative analysis between two groups.

### 8.3 Regression Tests
Maintain stable fixtures for:
- known bowtie-like inversion signal,
- known deletion-like depletion signal,
- no-signal negative control,
- parameter sensitivity around window-size changes.

### 8.4 Scientific QA
For any method-heavy change, verify:
- heatmap patterns remain interpretable,
- model-fit behavior is sensible,
- candidate counts do not shift unexpectedly without explanation,
- known benchmark examples still produce expected qualitative signals.

### 8.5 Performance QA
Benchmark:
- matrix build time,
- model-fit time,
- memory footprint,
- scaling with thread count.

---

## 9. Development Priorities for an Extensive Refactor / Expansion

When making substantial changes to WRATH, prioritize the following order:

1. **Preserve scientific behavior first**
   - Do not rewrite core logic before capturing current outputs with fixtures.

2. **Separate concerns**
   - Untangle barcode parsing, matrix construction, modeling, plotting, and CLI layers.

3. **Make intermediate data explicit**
   - Prefer named data structures over hidden state.

4. **Strengthen observability**
   - Add logging, summaries, and diagnostics before changing algorithms.

5. **Add tests before changing algorithms**
   - Especially for matrix construction and outlier detection.

6. **Avoid hidden methodological changes**
   - A refactor that alters results is not “just cleanup.”

7. **Keep R optional unless scientifically necessary**
   - Favor Python for core logic and keep cross-language boundaries narrow.

8. **Design for workflow integration**
   - Output formats should plug cleanly into BEDTools, plotting pipelines, and batch execution.

---

## 10. Working with Users: Core Principles

### Before starting, always read the project plan and inspect the relevant code paths.

### 1. Establish Context First
When a user asks for help:
- Ask what biological / analytical goal they are trying to accomplish.
- Ask what dataset, interval, or population comparison they are working on.
- Ask what exact command, parameters, and error / behavior they observed.

### 2. Diagnose Before Fixing (MOST IMPORTANT)

**DO NOT jump to conclusions and write lots of code before understanding the problem.**

Common mistakes to avoid:
- changing window-size logic before confirming the issue is actually parameterization,
- adding defensive exception handling that hides malformed BAM / barcode problems,
- changing statistical thresholds to “fix” noisy output without understanding coverage or molecule size,
- bundling CLI, modeling, and plotting changes into one patch.

**Correct process:**
1. Reproduce the issue.
2. Identify the failing stage.
3. Determine whether it is a software defect, data-quality problem, or method-limit problem.
4. Explain the root cause.
5. Make the smallest defensible change.
6. Re-run the relevant tests / fixture analyses.

### 3. Common Root Causes to Check First
Before editing code, check:
- missing or malformed barcode tags,
- BAM not indexed or not coordinate-sorted,
- wrong contig names,
- interval too small or too large for selected window size,
- mapQ filtering too aggressive,
- insufficient barcode density,
- method limitations caused by molecule size vs target SV size,
- plotting failures caused by empty or degenerate matrices.

### 4. Help Users Help Themselves
Encourage users to:
- inspect logs and intermediate outputs,
- test on a single locus before running whole genomes,
- compare outputs under only one parameter change at a time,
- distinguish biological ambiguity from software failure.

---

## 11. Common Issues and Troubleshooting

### Issue 1: No Candidate Signals in a Region Expected to Contain an SV
**Possible causes:**
- window size too coarse,
- SV smaller than what linked-read molecule length can resolve,
- sparse barcode coverage,
- interval mis-specified,
- mapQ filter removed too many reads.

**Diagnosis:**
1. Confirm contig and coordinates.
2. Inspect barcode density and read coverage.
3. Try a smaller window size only if coverage supports it.
4. Check whether the expected event is below the likely detection regime.

### Issue 2: Excessive False Positives Near the Diagonal
**Possible causes:**
- noisy local barcode sharing,
- insufficient Z-score filtering,
- low-complexity or repetitive sequence,
- model fit degraded by sparse data.

**Diagnosis:**
1. Inspect raw matrix.
2. Review Z-score threshold and fit diagnostics.
3. Check local coverage and repetitiveness.

### Issue 3: Heatmap Looks Empty or Uniform
**Possible causes:**
- barcode parsing failed,
- BAM lacks expected tags,
- interval has little usable data,
- plotting scale unsuitable for the matrix range.

### Issue 4: Runtime or Memory Explosion on Large Contigs
**Possible causes:**
- window size too small for chromosome-wide scan,
- matrix materialization too eager,
- repeated BAM traversal.

**Diagnosis:**
1. Start with a targeted region.
2. Increase window size for exploratory whole-contig work.
3. Profile matrix construction before rewriting algorithms.

---

## 12. For Claude Code (AI Assistant)

When helping with WRATH:

0. **Prepare**
   - Read the codebase, CLI entrypoints, and any method notes before proposing changes.

1. **Respect the scientific method**
   - Treat output changes as methodological changes until proven otherwise.

2. **Use domain vocabulary**
   - Say barcode, molecule, window, Jaccard, heatmap, interval, contig, candidate SV, outlier, diagonal, mapQ.

3. **Prefer minimal, testable changes**
   - Especially in matrix construction and model fitting.

4. **Do not erase interpretability**
   - Avoid replacing transparent logic with opaque abstractions unless there is a strong reason.

5. **Preserve the CLI-first workflow**
   - Do not center the project around notebooks or GUIs unless explicitly requested.

6. **Document method assumptions**
   - Especially around window-size choice, molecule size, and candidate interpretation.

7. **Separate refactor from method change**
   - If both are needed, do them in separate commits / PRs.

8. **Keep it honest**
   - If the data cannot support an inference, say so.

---

*This guide was created to help AI assistants effectively support development of the WRATH project. Last updated: March 2026.*

# Claude Code Guidelines by Sabrina Ramonov

## Implementation Best Practices

### 0 — Purpose

These rules ensure maintainability, safety, and developer velocity.
**MUST** rules are enforced by CI; **SHOULD** rules are strongly recommended.

---

### 1 — Before Coding

- **BP-1 (MUST)** Ask the user clarifying questions for ambiguous or high-impact work.
- **BP-2 (SHOULD)** Draft and confirm an approach for complex work.
- **BP-3 (SHOULD)** If two or more approaches exist, list clear pros and cons.

---

### 2 — While Coding

- **C-1 (MUST)** Name functions with existing domain vocabulary for consistency.
- **C-2 (SHOULD NOT)** Introduce classes when small testable functions suffice.
- **C-3 (SHOULD)** Prefer simple, composable, testable functions.
- **C-4 (SHOULD NOT)** Add comments except for critical scientific caveats; rely on self-explanatory code where possible.
- **C-5 (SHOULD NOT)** Extract a new function unless it will be reused, materially improves testability, or significantly improves readability.

---

### 3 — Testing

- **T-1 (MUST)** Separate pure-logic unit tests from file / BAM / plotting integration tests.
- **T-2 (SHOULD)** Prefer integration tests over heavy mocking.
- **T-3 (SHOULD)** Unit-test complex algorithms thoroughly.
- **T-4 (SHOULD)** Test the entire returned structure in one strong assertion when practical.

---

### 4 — Database

- **D-1** No guidelines yet.

---

### 5 — Code Organization

- **O-1** No guidelines yet.

---

### 6 — Tooling Gates

- **G-1** No guidelines yet.

---

### 7 - Git

- **GH-1 (MUST)** Use Conventional Commits format when writing commit messages.
- **GH-2 (SHOULD NOT)** Refer to Claude or Anthropic in commit messages.

---

## Writing Functions Best Practices

When evaluating whether a function you implemented is good or not, use this checklist:

1. Can you read the function and honestly follow what it's doing?
2. Does the function have very high cyclomatic complexity?
3. Are there clearer data structures or algorithms that would simplify it?
4. Are there unused parameters?
5. Are there unnecessary type casts that belong at the boundary instead?
6. Is the function testable without excessive mocking?
7. Does it hide non-trivial dependencies that should become parameters?
8. Brainstorm 3 better function names and confirm the current one is best.

IMPORTANT: you SHOULD NOT refactor out a separate function unless there is a compelling need, such as:
- the refactored function is used in more than one place,
- the refactored function becomes easily unit testable while the original function is not,
- the original function is genuinely too hard to follow.

## Writing Tests Best Practices

When evaluating whether a test you've implemented is good or not, use this checklist:

1. SHOULD parameterize inputs; do not hide important assumptions in magic literals.
2. SHOULD NOT add a test unless it can fail for a real defect.
3. SHOULD ensure the test description states exactly what the final assertion verifies.
4. SHOULD compare results to independent expectations or domain properties.
5. SHOULD follow the same lint, type-safety, and style rules as production code.
6. SHOULD express invariants or axioms when practical.
7. Unit tests for a function should be grouped under `describe(functionName, () => ...)`.
8. Use strong assertions.
9. SHOULD test edge cases, realistic input, unexpected input, and value boundaries.
10. SHOULD NOT test conditions that are already guaranteed by the type checker.

## Remember Shortcuts

### QNEW
When I type `qnew`, this means:

```text
Understand all best practices listed in CLAUDE.md.
Your code should always follow these best practices.
```

### QPLAN
When I type `qplan`, this means:

```text
Analyze similar parts of the codebase and determine whether your plan:
- is consistent with the rest of the codebase,
- introduces minimal changes,
- reuses existing code.
```

### QCODE
When I type `qcode`, this means:

```text
Implement your plan and make sure your new tests pass.
Always run tests to make sure you didn't break anything else.
```

### QCHECK
When I type `qcheck`, this means:

```text
You are a skeptical senior software engineer.
Perform this analysis for every major code change you introduced:
1. CLAUDE.md checklist Writing Functions Best Practices.
2. CLAUDE.md checklist Writing Tests Best Practices.
3. CLAUDE.md checklist Implementation Best Practices.
```

### QCHECKF
When I type `qcheckf`, this means:

```text
You are a skeptical senior software engineer.
Perform this analysis for every major function you added or edited:
1. CLAUDE.md checklist Writing Functions Best Practices.
```

### QCHECKT
When I type `qcheckt`, this means:

```text
You are a skeptical senior software engineer.
Perform this analysis for every major test you added or edited:
1. CLAUDE.md checklist Writing Tests Best Practices.
```

### QUX
When I type `qux`, this means:

```text
Imagine you are a human UX tester of the feature you implemented.
Output a comprehensive list of scenarios you would test, sorted by highest priority.
```

### QGIT
When I type `qgit`, this means:

```text
Add all changes to staging, create a commit, and push to remote.
Use Conventional Commits format and do not mention Claude or Anthropic.
```
