# Adding User Sequences to Nextstrain Oropouche Workflow After Running Juno

## Method 1: Assembly Pipeline Import (Recommended)

### Prerequisites
Run the Nextstrain ingest workflow first:
   ```bash
   cd ingest/
   nextstrain build .
   cd ..
   ```

### Usage

1. **Create a metadata CSV template** (optional):
   ```bash
   python scripts/add_user_sequences.py --template my_samples.csv
   ```



   **Notes:**
   - Leave accession columns empty for segments that didn't pass QC or aren't available
   - Date formats accepted: `YYYY-MM-DD`, `YYYY-MM`, `YYYY`, or `MM/DD/YYYY` (will be auto-converted)
   - Region, country, division, host, authors, and institution are auto-populated for Florida

2. **Import sequences and rebuild phylogeny**:
   ```bash
   python scripts/add_user_sequences.py \
     --assembly-dir path/to/juno_output \
     --metadata my_samples.csv \
     --rebuild
   ```

3. **View results**:
   ```bash
   cd phylogenetic/
   nextstrain view .
   ```

   Or manually run phylogenetic build::
   ```bash
   cd phylogenetic/
   nextstrain build .
   ```

### What Gets Auto-Populated
- **QC Status**: Automatically set based on `summary_report.tsv` (PASS → "good", PASS_W_HIGH_N_BASES → "mediocre")
- **Segment Lengths**: Read from assembly output
- **Release Date**: Set to today's date (when script runs)
- **Metadata**: Region (North America), Country (USA), Division (Florida), Host (Homo sapiens), Authors (BPHL et al.), Institution (Florida BPHL)

---
