# Windows Executable

Latest release: in_silico_project_ver3.5.1.exe

※ Download from Releases → Assets

Note: Do not use “Code → Download ZIP”. That only contains the source code and may include LFS pointer files instead of the real executable or large binary assets.

# NGS Dataset Downloader (NDD) & In Silico Sequence Mining (ISSM)

A graphical Python tool for mining, matching, and downloading NGS data (SRA/ENA) with probe-based sequence detection, IUPAC-aware matching, sampling support, and optimized exact/biological matching for FASTQ/FASTA datasets.

---

## Project Overview

This tool was developed to assist molecular biologists and bioinformaticians in mining NGS datasets and verifying probe-target matching in silico.  
It provides an all-in-one graphical interface for downloading and analyzing sequencing data using custom probe sequences, enabling faster and more reproducible in silico screening workflows.

ISSM supports rapid preliminary detection of target nucleotide sequences using existing PCR primers, degenerate primers, hybridization probes, or custom-designed probe sequences.  
Matched reads or sequences can be exported for downstream validation using external alignment tools and genomic-context assessment.

---

## Features

- **In Silico Sequence Matching**:
  - Probe-based sequence detection using user-defined FASTA-format probes
  - Forward and reverse-complement probe matching
  - Support for IUPAC degenerate nucleotide codes in probe and input sequences
  - FASTQ and FASTA input support
  - gzip-compressed and plain-text input support
  - Sampling and batch detection modes
  - Matched read/sequence export in FASTA format

- **Matching Methods in ISSM**:
  - **Exact matching** is automatically used when threshold is 100%
  - **Fast biological matching** is available when threshold is below 100%
    - Uses probe-length percent identity
    - Allows mismatch-tolerant detection based on probe length and threshold
    - Suitable for primer/probe validation and variant-tolerant detection
  - **Legacy fuzzy matching** is available for compatibility
    - Uses the original RapidFuzz partial-ratio based fuzzy matching
    - Suitable for reproducing previous ISSM fuzzy-matching results

- **NGS Downloader GUI**:
  - Search and fetch data from **NCBI SRA** or **EBI ENA**
  - Automated SRA Toolkit download & setup
  - Real-time download status and logging
  - Organism-level search result statistics and filtering

- **Performance**:
  - Rust-based binary matcher for fast exact and biological matching
  - SeqKit-assisted strict sampling
  - Automatic CPU/thread configuration
  - File-level and matcher-level parallel behavior optimized internally
  - Progress reporting based on processed records and file byte progress

---

## GUI Architecture Overview

This project contains **two independent modules**, both built using **PySide6 / Qt Designer**.

### In Silico Sequence Mining (ISSM)
- Built with **PySide6 / Qt Designer**
- Modern interface with probe editor, file selection panels, matching method controls, sampling controls, and batch status
- Modules:
  - `in_silico_sequence_mining.py`
- External tools used by ISSM:
  - `tools/seqkit.exe`
  - `tools/issm_matcher.exe`

### NGS Dataset Downloader (NDD)
- Built with **PySide6 / Qt Designer**
- Lightweight UI for searching, filtering, and downloading datasets from SRA/ENA
- Modules:
  - `NGS_data_downloader.py`

---

## Usage Example

### ISSM:
1. Select the input format:
   - FASTQ reads
   - FASTA sequences
2. Select an input folder containing supported FASTQ or FASTA files.
3. Paste probes in FASTA format:
    ```
    >probe1
    AGTCAGTC
    >probe2
    TTAGGCCA
    ```
4. Set threshold and sampling.
5. If threshold is below 100%, choose the matching method:
   - Fast biological matching
   - Legacy fuzzy matching
6. Choose output path and click “Start Analysis”.

### Supported ISSM Input Formats

| Input type | Supported extensions |
|-----------|----------------------|
| FASTQ | `.fastq.gz`, `.fastq`, `.fq.gz`, `.fq` |
| FASTA | `.fasta.gz`, `.fasta`, `.fa.gz`, `.fa`, `.fna.gz`, `.fna`, `.fas.gz`, `.fas` |

File extension matching is case-insensitive.

Note: FASTQ and FASTA files should be analyzed as separate runs. A single ISSM run uses either FASTQ mode or FASTA mode to avoid ambiguity in sampling and result interpretation.

### ISSM Threshold Interpretation

| Condition | Matching behavior | Threshold meaning |
|----------|------------------|-------------------|
| Threshold = 100% | Exact matching | Exact probe match only |
| Threshold < 100%, Fast biological matching | Rust-based biological matching | Probe-length percent identity |
| Threshold < 100%, Legacy fuzzy matching | Original Python/RapidFuzz matching | RapidFuzz partial-ratio score |

For Fast biological matching, allowed mismatches are calculated based on probe length:

```txt
allowed_mismatches = floor(probe_length × (100 - threshold) / 100)
```

Example:

```txt
20 bp probe, threshold 95% = up to 1 mismatch
20 bp probe, threshold 90% = up to 2 mismatches
20 bp probe, threshold 97% = 0 mismatches
```

### ISSM IUPAC Handling Policy

ISSM supports IUPAC degenerate nucleotide codes in probe sequences and input FASTQ/FASTA sequences, but the interpretation differs by source.

| Sequence source | IUPAC handling | Note |
|----------------|----------------|------|
| Probe sequence | Full IUPAC support | Degenerate bases such as `R`, `Y`, `K`, `M`, `S`, `W`, `B`, `D`, `H`, `V`, and `N` are treated as intended probe ambiguity. |
| Input FASTQ/FASTA sequence | IUPAC support except `N`/`I` wildcard matching | Ambiguous bases such as `R`, `Y`, `K`, `M`, `S`, `W`, `B`, `D`, `H`, and `V` are interpreted as ambiguity codes. |
| Input `N` or `I` | Treated as no-call bases | Input `N`/`I` are not treated as wildcards and do not automatically match probe bases. |

### NDD:
1. Select database (SRA or ENA)
2. Enter query and email
3. Search and select accessions
4. Choose output format and click “Download Start”

### Detailed instructions and troubleshooting
- For detailed instructions and troubleshooting, please refer to the **In silico project(ver3.5.1)_User_manual.pdf** included in the repository or release package.

---

## Tested Environment

- OS: Windows 10, Windows 11
- RAM: 8 GB or higher recommended
- Disk: >5 GB free space for downloaded FASTQ files
- Internet: Required for SRA/ENA data retrieval

---

## File Structure

| File | Description |
|------|-------------|
| `in_silico_project.py` | Main launcher |
| `in_silico_sequence_mining.py` | Probe matching engine and GUI logic (ISSM) |
| `NGS_data_downloader.py` | Dataset downloader GUI (NDD) |
| `ISSM.py`, `resource_rc.py` | PySide6-based UI and resources |
| `tools/seqkit.exe` | External SeqKit binary used for strict sampling and sequence utilities |
| `tools/issm_matcher.exe` | Rust-based ISSM matcher used for exact and fast biological matching |
| `.ico`, `.png` files | Icon & resource files |
| `*_ver*.exe` | Windows packaged executables |

---

## Requirements when you use Python code

```txt
biopython
pyside6
matplotlib
seaborn
pandas
numpy
rapidfuzz
psutil
requests
aiofiles
certifi
beautifulsoup4
pillow
pyinstaller
```

External binaries required for full ISSM functionality:

```txt
tools/seqkit.exe
tools/issm_matcher.exe
```

If building from source, use the same Python environment for installing packages and running PyInstaller:

```bat
python -m pip install -r requirements_build_v8.txt
python -m PyInstaller --clean --noconfirm in_silico_project_onefile.spec
```

---

## License

This project is licensed under the **Modified MIT**.  
See the LICENSE file for details.

---

## Bundled Tools

| Tool | Description | License |
|------|-------------|---------|
| SeqKit | External FASTA/FASTQ processing tool used for strict sampling and sequence utilities | MIT |
| issm_matcher | Rust-based internal matcher used by ISSM for exact and fast biological matching | Same as this project (Modified MIT) |

---

## Third-Party Python Libraries

| Library | License |
|---------|---------|
| Biopython | Biopython License |
| PySide6 | LGPL 3.0 |
| psutil | BSD |
| matplotlib | PSF |
| seaborn | BSD |
| pandas | BSD |
| numpy | BSD |
| rapidfuzz | MIT |
| requests | Apache 2.0 |
| aiofiles | Apache 2.0 |
| certifi | MPL 2.0 |
| BeautifulSoup4 | MIT |
| Pillow (PIL) | HPND |

---

## Third-Party Rust Dependencies

The Rust-based `issm_matcher` component directly uses the following crates.
Dependencies marked as `MIT OR Apache-2.0` are dual-licensed; users may comply with either license option where applicable.

| Crate | License |
|-------|---------|
| anyhow | MIT OR Apache-2.0 |
| clap | MIT OR Apache-2.0 |
| flate2 | MIT OR Apache-2.0 |
| rayon | MIT OR Apache-2.0 |

---

## Maintainer

Developed and maintained by [Hye Kwon Kim, Min Chan Kim, Hye Ji Jung or Viral Infectious Disease Laboratory, Chungbuk National University, Republic of Korea (Professor Hye Kwon Kim's Laboratory)]

For collaboration inquiries, permission requests, or file/access-related requests: **khk1329@chungbuk.ac.kr**

For questions, feedback, or bug reports: **kminchan1010@naver.com**, **ghkrk007512@naver.com**, **khk1329@chungbuk.ac.kr**

---

Made for molecular biology and bioinformatics research.

---

## Release Notes

### v3.5.1
1. ISSM
- Removed redundant reverse-complement calculations during sequence/read matching.
- Added SeqKit path detection and execution wrapper functions.
- Added SeqKit-assisted read/record counting and strict sampling.
- Added sequence text stream handling for compressed and plain-text input files.
- Reduced repeated `probe_seq.upper()` calls during matching.
- Added precomputed probe metadata generation.
- Replaced `fuzz.partial_ratio` with direct substring matching when threshold is 100% in the Python path.
- Increased batch size for faster detection.
- Treated `I` in probe sequences similarly to `N`.
- Generated reverse-complement sequences only for probe sequences.
- Added a Rust-based binary matcher, `tools/issm_matcher.exe`, for fast exact matching and fast biological matching.
- Threshold 100% now uses the Rust matcher instead of the original Python/RapidFuzz path.
- Threshold below 100% supports two matching methods:
  - Fast biological matching: probe-length percent identity and allowed mismatch-based matching.
  - Legacy fuzzy matching: original RapidFuzz partial-ratio based matching for compatibility.
- Added IUPAC-aware matching using a bitmask strategy.
- Applied full IUPAC handling to probe sequences.
- Applied input-side IUPAC handling while treating input `N`/`I` as no-call bases rather than wildcards to reduce false-positive matches from unresolved input bases.
- Added support for gzip-compressed and plain-text FASTQ/FASTA files.
- Added support for `.fastq.gz`, `.fastq`, `.fq.gz`, `.fq`, `.fasta.gz`, `.fasta`, `.fa.gz`, `.fa`, `.fna.gz`, `.fna`, `.fas.gz`, and `.fas` files.
- Added Rust matcher progress reporting using processed record count, read byte count, and total file byte count.
- Added automatic matcher thread control. Fast matcher mode processes one file at a time while using internal matcher threads.
- Removed the basic GUI custom parallel processing option to avoid mode-dependent interpretation and unstable user-defined thread settings.
- Added warning dialog for sampling values between 51% and 99%, because strict sampling in this range may provide limited speed benefit and can be slower than 100% analysis.
- Fixed sampling warning dialog behavior when changing sampling to 50% or 100%.
- Matching method controls are shown only when threshold is below 100%; threshold 100% automatically uses exact matching.
- Fixed GUI spacing when switching threshold from 100% to below 100% and back to 100%.
- Improved completion and cancellation handling so that the final completion message appears only after plot, heatmap, and summary generation are finished.
- Improved cancellation handling so that cleanup is completed before the final cancelled status is shown.
- Improved result and summary integration for Rust matcher outputs.
- Applied ISSM-specific GUI font and widget-size stabilization for packaged executable builds.

### v3.4.1
1. NDD
- Patches to minimize certificate-related issues during SRA-based search and download in NDD have been applied.
- If certificate errors still occur when using version 3.4.1, refer to the Troubleshooting section of the manual.

### v3.4.0
1. ISSM
- Validate required inputs before Run Analysis; launch process_window only after passing checks.
- Capped maximum CPU workers at 61.
- The final completion message is emitted in the process window only when all files reach a terminal state; counts include both Completed + Cancelled.

2. NDD
- Show Organism distribution stats for search results and enable filtering by organism.
- Write start and end timestamps to the log file.
- Migrated the launcher and NGS Dataset Downloader (NDD) GUI from Tkinter to PySide6 for a unified, modern UI.
- Switched NDD internal signaling to thread-safe queues.
- Enabled Select All button in Selected list box.
- The final completion message is emitted in the process window only when all files reach a terminal state; counts include both Completed + Cancelled.
