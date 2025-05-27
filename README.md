# In Silico Sequence Mining Tool

A graphical Python tool for mining, matching, and downloading NGS data (SRA/ENA) with probe-based sequence detection and parallel processing support.

---

## Project Overview

This tool was developed to assist molecular biologists and bioinformaticians in mining NGS datasets and verifying probe-target matching in silico.  
It provides an all-in-one graphical interface for downloading and analyzing sequencing data using custom probe sequences, enabling faster and more reproducible in silico screening workflows.

---

## Features

- **In Silico Sequence Matching**:
  - Fuzzy probe matching (original + reverse complement)
  - Support for IUPAC nucleotide codes
  - FASTA and FASTQ input support
  - Sampling and batch detection modes

- **NGS Downloader GUI**:
  - Search and fetch data from **NCBI SRA** or **EBI ENA**
  - Automated SRA Toolkit download & setup
  - Real-time download status and logging

- **Performance**:
  - Multiprocessing and threading enabled
  - Dynamic memory and CPU usage control
  - Auto batch size and worker tuning

---

## GUI Architecture Overview

This project contains **two independent modules**, each built using a different GUI framework:

### In Silico Sequence Mining (ISSM)
- Built with **PySide6 / Qt Designer**
- Modern interface with probe editor, file selection panels, and batch status
- Modules:
  - `in_silico_sequence_mining.py`

### NGS Dataset Downloader (NDD)
- Built with **Tkinter**
- Lightweight UI for searching, filtering, and downloading datasets from SRA/ENA
- Modules:
  - `NGS_data_downloader.py`

---

## GUI Preview

*(Add screenshots here of both ISSM and NDD GUIs)*

---

## Windows Executable
- `in_silico_project_ver3.2.4.exe`: General user version  

---

## Usage Example

### ISSM:
1. Select input folder containing `.fastq` or `.fasta` files
2. Paste probes in FASTA format:
    ```
    >probe1
    AGTCAGTC
    >probe2
    TTAGGCCA
    ```
3. Set threshold and sampling
4. Choose output path and click “Start Analysis”

### NDD:
1. Select database (SRA or ENA)
2. Enter query and email
3. Search and select accessions
4. Choose output format and click “Download Start”

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
| `in_silico_sequence_mining.py` | Probe matching engine (ISSM) |
| `NGS_data_downloader.py` | Dataset downloader GUI (NDD) |
| `ISSM.py`, `resource_rc.py` | PySide6-based UI |
| `.ico`, `.png` files | Icon & resource files |
| `*_ver*.exe` | Windows packaged executables |

---

## Requirements when you used python cord

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
beautifulsoup4
pillow
```

---

## License

This project is licensed under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

> PySide6 is under LGPL 3.0 — dynamic linking is recommended for distribution.

---

## Third-Party Libraries

| Library         | License      |
|----------------|--------------|
| Biopython      | BSD          |
| PySide6        | LGPL 3.0     |
| psutil         | BSD          |
| matplotlib     | PSF          |
| seaborn        | BSD          |
| pandas         | BSD          |
| numpy          | BSD          |
| rapidfuzz      | MIT          |
| requests       | Apache 2.0   |
| BeautifulSoup4 | MIT          |
| Pillow (PIL)   | HPND         |

---

## Maintainer

Developed and maintained by [Hye Kwon Kim, Min Chan Kim, Hye Ji Jung or Viral Infectious Disease Laboratory, Chungbuk National University, Republic of Korea (Professor Hye Kwon Kim's Laboratory)]

For questions, feedback, or bug reports: **kminchan1010@naver.com**  

---

Made for molecular biology and bioinformatics research.
