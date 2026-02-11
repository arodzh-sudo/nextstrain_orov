# Nextstrain repository for Oropouche virus

> **Note**: This repository has been modified by Florida's Bureau of Public Health Laboratories to support integration of local sequence data alongside public GenBank sequences.

This repository contains two workflows for the analysis of Oropouche virus data:

- [`ingest/`](./ingest) - Download data from GenBank, clean and curate it
- [`phylogenetic/`](./phylogenetic) - Filter sequences, align, construct phylogeny and export for visualization

Each folder contains a README.md with more information.

## Installation

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Quickstart

### Standard workflow (GenBank data only)

Run the default phylogenetic workflow via:
```bash
cd phylogenetic/
nextstrain build .
nextstrain view .
```

### Adding local sequences from assembly pipeline

This workflow supports adding sequences from the [BPHL's Juno reference-based assembly pipeline](https://github.com/BPHL-Molecular/Juno):

1. **Run the ingest workflow** to download GenBank data:
   ```bash
   cd ingest/
   nextstrain build .
   cd ..
   ```

2. **Create a metadata CSV** for local sequences:
   ```bash
   python scripts/add_user_sequences.py --template my_samples.csv
   ```
   
   Edit the CSV with strain names, accessions, and collection dates.

3. **Import sequences** from Juno output:
   ```bash
   python scripts/add_user_sequences.py \
     --assembly-dir /path/to/consensus \
     --metadata my_samples.csv \
     --rebuild
   ```

4. **View results**:
   ```bash
   cd phylogenetic/
   nextstrain view .
   ```

## Documentation

- [Running a pathogen workflow](https://docs.nextstrain.org/en/latest/tutorials/running-a-workflow.html)
