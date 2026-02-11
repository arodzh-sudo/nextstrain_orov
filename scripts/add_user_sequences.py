#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import sys
import shutil
import re
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
import subprocess

class AssemblyImporter:

    # Reference to segment mapping
    REFERENCE_TO_SEGMENT = {
        'PQ064919.1': 'L',
        'PQ064920.1': 'M',
        'PQ064921.1': 'S'
    }

    # Valid QC status values
    VALID_QC_STATUS = ['PASS', 'PASS_W_HIGH_N_BASES']

    # Default metadata for all samples
    DEFAULT_METADATA = {
        'region': 'North America',
        'country': 'USA',
        'division': 'Florida',
        'host': 'Homo sapiens',
        'authors': 'BPHL et al.',
        'institution': 'Florida BPHL'
    }

    def __init__(self, assembly_dir):
        self.assembly_dir = Path(assembly_dir)
        self.summary_file = self.assembly_dir / "summary_report.tsv"
        self.consensus_dir = self.assembly_dir / "consensus"

        # Validate paths
        if not self.assembly_dir.exists():
            raise FileNotFoundError(f"Assembly directory not found: {self.assembly_dir}")
        if not self.summary_file.exists():
            raise FileNotFoundError(f"Summary report not found: {self.summary_file}")
        if not self.consensus_dir.exists():
            raise FileNotFoundError(f"Consensus directory not found: {self.consensus_dir}")

    def parse_summary_report(self):

        df = pd.read_csv(self.summary_file, sep='\t')

        # Filter QC status
        good_quality = df[df['qc_pass_fail'].isin(self.VALID_QC_STATUS)].copy()

        # Add segment column
        good_quality['segment'] = good_quality['reference'].map(self.REFERENCE_TO_SEGMENT)

        # Extract strain name (remove reference suffix)
        good_quality['strain'] = good_quality['sampleID'].str.replace(
            r'_PQ064(919|920|921)\.1$', '', regex=True
        )

        return good_quality

    def group_by_strain(self, qc_df):

        strain_data = {}

        for _, row in qc_df.iterrows():
            strain = row['strain']
            segment = row['segment']

            if strain not in strain_data:
                strain_data[strain] = {}

            strain_data[strain][segment] = row

        return strain_data

    def validate_date_format(self, date_str):
        date_str = str(date_str).strip()

        # Already in correct format
        if re.match(r'^\d{4}(-\d{2}(-\d{2})?)?$', date_str):
            # Validate it's a real date
            try:
                if len(date_str) == 10:
                    datetime.strptime(date_str, '%Y-%m-%d')
                elif len(date_str) == 7:
                    datetime.strptime(date_str, '%Y-%m')
                elif len(date_str) == 4:
                    int(date_str)
                return date_str
            except ValueError:
                raise ValueError(f"Invalid date: {date_str}")

        # Try common formats and convert
        formats = [
            ('%m/%d/%Y', '%Y-%m-%d'),
            ('%m/%d/%y', '%Y-%m-%d'),
            ('%Y/%m/%d', '%Y-%m-%d'),
        ]

        for in_fmt, out_fmt in formats:
            try:
                dt = datetime.strptime(date_str, in_fmt)
                return dt.strftime(out_fmt)
            except ValueError:
                continue

        raise ValueError(
            f"Unrecognized date format: {date_str}. Use YYYY-MM-DD, YYYY-MM, or YYYY"
        )

    def validate_user_metadata(self, metadata_csv, strain_data):

        user_df = pd.read_csv(metadata_csv)

        # Check required columns
        required = ['strain', 'date']
        missing = [col for col in required if col not in user_df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Check each strain
        for idx, row in user_df.iterrows():
            strain = row['strain']

            if strain not in strain_data:
                raise ValueError(
                    f"Row {idx}: Strain '{strain}' not found in QC-passed sequences\n"
                    f"Available strains: {list(strain_data.keys())}"
                )

            # Validate date
            try:
                normalized_date = self.validate_date_format(row['date'])
                user_df.at[idx, 'date'] = normalized_date
            except ValueError as e:
                raise ValueError(f"Row {idx}: {e}")

            # Check accessions match available segments
            available_segments = strain_data[strain].keys()
            for seg in ['L', 'M', 'S']:
                accession_col = f'accession_{seg}'
                has_accession = (accession_col in row and
                               pd.notna(row[accession_col]) and
                               str(row[accession_col]).strip() != '')
                has_segment = seg in available_segments

                if has_accession and not has_segment:
                    print(f"⚠️  Warning row {idx}: {strain} has accession for {seg} "
                          f"but segment did not pass QC (will be skipped)")
                if has_segment and not has_accession:
                    print(f"⚠️  Warning row {idx}: {strain} has {seg} segment that passed QC "
                          f"but no accession provided (will be skipped)")

        return user_df

    def read_consensus_sequence(self, sampleID, reference):

        fasta_file = self.consensus_dir / f"{sampleID}.consensus.fasta"

        if not fasta_file.exists():
            raise FileNotFoundError(f"Consensus file not found: {fasta_file}")

        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if len(sequences) == 0:
            raise ValueError(f"No sequences in {fasta_file}")

        seq_str = str(sequences[0].seq)
        return seq_str, len(seq_str)

    def create_nextstrain_metadata_row(self, user_row, strain_data_dict, strain_qc_df):

        strain = user_row['strain']
        collection_date = user_row['date']
        release_date = datetime.now().strftime('%Y-%m-%d')

        # Count segments that have both QC pass AND accession
        n_segments = 0
        for seg in ['L', 'M', 'S']:
            has_segment = seg in strain_data_dict
            accession_col = f'accession_{seg}'
            has_accession = (accession_col in user_row and
                           pd.notna(user_row[accession_col]) and
                           str(user_row[accession_col]).strip() != '')
            if has_segment and has_accession:
                n_segments += 1

        # Base metadata (auto-populated)
        metadata = {
            'strain': strain,
            'n_segments': n_segments,
            'date': collection_date,
            'region': self.DEFAULT_METADATA['region'],
            'country': self.DEFAULT_METADATA['country'],
            'division': self.DEFAULT_METADATA['division'],
            'location': '',
            'host': self.DEFAULT_METADATA['host'],
            'authors': self.DEFAULT_METADATA['authors'],
            'abbr_authors': self.DEFAULT_METADATA['authors'],
            'institution': self.DEFAULT_METADATA['institution'],
        }

        # Initialize all segment columns as empty/zero
        # Use pd.NA for numeric columns to avoid export issues
        for seg in ['L', 'M', 'S']:
            metadata.update({
                f'accession_{seg}': '',
                f'accession_version_{seg}': '',
                f'length_{seg}': pd.NA,
                f'date_released_{seg}': '',
                f'date_updated_{seg}': '',
                f'sra_accessions_{seg}': '',
                f'segment_{seg}': '0',
                f'qc_{seg}': ''
            })

        # Fill segment-specific data
        for segment, qc_row in strain_data_dict.items():
            accession_col = f'accession_{segment}'
            accession = user_row.get(accession_col, '')

            # Only add if user provided accession
            if accession and str(accession).strip():
                accession = str(accession).strip()
                # Add version if not already present
                accession_version = accession if '.' in accession else f"{accession}.1"

                # Map QC status
                qc_status = 'good' if qc_row['qc_pass_fail'] == 'PASS' else 'mediocre'

                metadata.update({
                    f'accession_{segment}': accession.split('.')[0],
                    f'accession_version_{segment}': accession_version,
                    f'length_{segment}': int(qc_row['assembly_length']),
                    f'date_released_{segment}': release_date,
                    f'date_updated_{segment}': release_date,
                    f'segment_{segment}': '1',
                    f'qc_{segment}': qc_status
                })

        return metadata

class UserSequenceProcessor:
    def __init__(self, base_dir="."):
        self.base_dir = Path(base_dir)
        self.ingest_dir = self.base_dir / "ingest"
        self.phylo_dir = self.base_dir / "phylogenetic"
        self.results_dir = self.ingest_dir / "results"

    def create_assembly_template(self, output_file):
        """Create a simplified CSV template for assembly pipeline import"""
        template_data = {
            'strain': ['Sample_001', 'Sample_002', 'Sample_003'],
            'accession_L': ['ACC001', 'ACC004', 'ACC007'],
            'accession_M': ['ACC002', 'ACC005', 'ACC008'],
            'accession_S': ['ACC003', 'ACC006', 'ACC009'],
            'date': ['2026-01-15', '2026-01-20', '2026-01-25']
        }

        df = pd.DataFrame(template_data)
        df.to_csv(output_file, index=False)

        print(f"✓ Template created: {output_file}")
        print("\nInstructions:")
        print("1. Replace example data with your actual strain names and accessions")
        print("2. Strain names must match those in your assembly pipeline output")
        print("3. Leave accession columns empty for segments not available")
        print("4. Date format: YYYY-MM-DD, YYYY-MM, YYYY, or MM/DD/YYYY")
        print("5. Save and run:")
        print(f"   python scripts/add_user_sequences.py --assembly-dir <path> --metadata {output_file}")
        print("\nNote: Region, country, division, host, authors, and institution are auto-populated.")

        return output_file

    def create_template(self, output_file):
        """Create a CSV template for manual user input (legacy)"""
        template_data = {
            'strain': ['MyStrain_001', 'MyStrain_002'],
            'segment': ['S', 'S'],
            'sequence_file': ['sequences/strain001.fasta', 'sequences/strain002.fasta'],
            'date': ['2024-03-XX', '2024-04-XX'],
            'region': ['North America', 'North America'],
            'country': ['USA', 'USA'],
            'division': ['Florida', 'Florida'],
            'location': ['Miami', 'Tampa'],
            'host': ['Homo sapiens', 'Homo sapiens'],
            'authors': ['John Doe et al.', 'John Doe et al.'],
            'institution': ['University of Florida', 'University of Florida'],
            'notes': ['', '']
        }

        df = pd.DataFrame(template_data)
        df.to_csv(output_file, index=False)

        print(f"Template created: {output_file}")
        print("\nInstructions:")
        print("1. Open the CSV file in Excel or any spreadsheet software")
        print("2. Replace the example data with your actual sequence information")
        print("3. Make sure sequence files exist in the specified paths")
        print("4. Save as CSV and run: python scripts/add_user_sequences.py --input your_file.csv")

        return output_file

    def validate_csv(self, csv_file):
        """Validate the user CSV file"""
        print(f"Validating {csv_file}...")

        # Required columns
        required_cols = ['strain', 'segment', 'sequence_file', 'date', 'country', 'authors', 'institution']

        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            print(f"Error reading CSV file: {e}")
            return False

        # Check required columns
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Missing required columns: {missing_cols}")
            return False

        # Check for empty required fields
        for col in required_cols:
            if df[col].isna().any() or (df[col] == '').any():
                empty_rows = df[df[col].isna() | (df[col] == '')].index.tolist()
                print(f"Empty values in required column '{col}' at rows: {empty_rows}")
                return False

        # Check sequence files exist
        for idx, row in df.iterrows():
            seq_file = Path(row['sequence_file'])
            if not seq_file.exists():
                print(f"Sequence file not found: {seq_file} (row {idx})")
                return False

            # Validate FASTA format
            try:
                sequences = list(SeqIO.parse(seq_file, "fasta"))
                if len(sequences) == 0:
                    print(f"No sequences found in {seq_file} (row {idx})")
                    return False
                elif len(sequences) > 1:
                    print(f"Multiple sequences in {seq_file} (row {idx}) - will use first sequence")

                # Check minimum length for S segment
                seq_len = len(sequences[0].seq)
                if row['segment'] == 'S' and seq_len < 700:
                    print(f"Sequence in {seq_file} is {seq_len}bp, minimum for S segment is 700bp")

            except Exception as e:
                print(f"Error reading FASTA file {seq_file}: {e}")
                return False

        print(f"Validation passed! Found {len(df)} sequences to process")
        return True

    def process_sequences(self, csv_file, dry_run=False):
        """Process user sequences and add to Nextstrain data"""
        if not self.validate_csv(csv_file):
            return False

        print(f"Processing sequences from {csv_file}...")

        # Load user data
        user_df = pd.read_csv(csv_file)

        # Load existing metadata
        existing_metadata = self.results_dir / "metadata.tsv"
        if existing_metadata.exists():
            nextstrain_df = pd.read_csv(existing_metadata, sep='\t')
            print(f"Found {len(nextstrain_df)} existing strains in metadata")
        else:
            print("No existing metadata found. Run ingest workflow first.")
            return False

        # Process each user sequence
        new_rows = []
        sequences_to_add = {}  # {segment: [(strain, sequence), ...]}

        for idx, row in user_df.iterrows():
            print(f"Processing {row['strain']} ({row['segment']} segment)...")

            # Read sequence
            seq_file = Path(row['sequence_file'])
            sequences = list(SeqIO.parse(seq_file, "fasta"))
            sequence = sequences[0]

            # Create metadata row in Nextstrain format
            metadata_row = self.create_metadata_row(row, sequence)
            new_rows.append(metadata_row)

            # Store sequence for adding to FASTA files
            segment = row['segment']
            if segment not in sequences_to_add:
                sequences_to_add[segment] = []
            sequences_to_add[segment].append((row['strain'], str(sequence.seq)))

        if dry_run:
            print("🔍 DRY RUN - No files modified")
            print(f"Would add {len(new_rows)} new strains:")
            for row in new_rows:
                print(f"  - {row['strain']} ({row.get('segment_S', 0) + row.get('segment_M', 0) + row.get('segment_L', 0)} segments)")
            return True

        # Add metadata rows
        updated_df = pd.concat([nextstrain_df, pd.DataFrame(new_rows)], ignore_index=True)
        updated_df.to_csv(existing_metadata, sep='\t', index=False)
        print(f"Updated metadata with {len(new_rows)} new strains")

        # Add sequences to FASTA files
        for segment, sequences in sequences_to_add.items():
            fasta_file = self.results_dir / segment / "sequences.fasta"
            with open(fasta_file, 'a') as f:
                for strain, sequence in sequences:
                    f.write(f">{strain}\n{sequence}\n")
            print(f"Added {len(sequences)} sequences to {segment} segment")

        print("Successfully added all user sequences!")
        return True

    def process_assembly_output(self, assembly_dir, metadata_csv, dry_run=False):
        """
        Process sequences from assembly pipeline output.

        Args:
            assembly_dir: Path to assembly output directory
            metadata_csv: User's simplified metadata CSV (strain, accessions, date)
            dry_run: Preview without making changes
        """
        print(f"Processing assembly output from {assembly_dir}...\n")

        try:
            # Initialize importer
            importer = AssemblyImporter(assembly_dir)

            # Parse QC summary
            qc_df = importer.parse_summary_report()
            print(f"✓ Found {len(qc_df)} sequences that passed QC")

            # Group by strain
            strain_data = importer.group_by_strain(qc_df)
            print(f"✓ Found {len(strain_data)} strains with passing segments:")
            for strain, segments in strain_data.items():
                seg_list = ', '.join(sorted(segments.keys()))
                print(f"    {strain}: {seg_list}")

            print()

            # Validate user metadata
            user_df = importer.validate_user_metadata(metadata_csv, strain_data)
            print(f"✓ Validated user metadata: {len(user_df)} strains\n")

        except (FileNotFoundError, ValueError) as e:
            print(f"❌ Error: {e}")
            return False

        # Load existing Nextstrain metadata
        existing_metadata = self.results_dir / "metadata.tsv"
        if not existing_metadata.exists():
            print("❌ Error: No existing metadata found.")
            print("   Run the ingest workflow first: cd ingest && nextstrain build .")
            return False

        nextstrain_df = pd.read_csv(existing_metadata, sep='\t')

        # Check for duplicates
        existing_strains = set(nextstrain_df['strain'].values)
        duplicate_strains = [row['strain'] for _, row in user_df.iterrows()
                             if row['strain'] in existing_strains]

        if duplicate_strains:
            print(f"❌ Error: Duplicate strains found: {duplicate_strains}")
            print("   These strains already exist in Nextstrain metadata.")
            return False

        # Process each strain
        new_metadata_rows = []
        sequences_to_add = {'L': [], 'M': [], 'S': []}

        print("Processing strains:")
        for _, user_row in user_df.iterrows():
            strain = user_row['strain']
            print(f"\n  {strain}:")

            # Get QC data for this strain
            strain_qc_dict = strain_data[strain]

            # Read consensus sequences for segments with accessions
            strain_sequences = {}
            for segment, qc_row in strain_qc_dict.items():
                accession_col = f'accession_{segment}'
                has_accession = (accession_col in user_row and
                               pd.notna(user_row[accession_col]) and
                               str(user_row[accession_col]).strip() != '')

                if has_accession:
                    try:
                        # Get reference for this segment
                        reference = [k for k, v in importer.REFERENCE_TO_SEGMENT.items()
                                   if v == segment][0]

                        # Read sequence
                        seq, length = importer.read_consensus_sequence(
                            qc_row['sampleID'],
                            reference
                        )
                        strain_sequences[segment] = (seq, length)

                        qc_label = qc_row['qc_pass_fail']
                        print(f"    ✓ {segment} segment: {length}bp ({qc_label})")

                        # Store for adding to FASTA
                        sequences_to_add[segment].append((strain, seq))

                    except Exception as e:
                        print(f"    ❌ Error reading {segment} segment: {e}")
                        return False
                else:
                    print(f"    - {segment} segment: skipped (no accession)")

            # Create metadata row
            try:
                metadata_row = importer.create_nextstrain_metadata_row(
                    user_row,
                    strain_qc_dict,
                    qc_df[qc_df['strain'] == strain]
                )
                new_metadata_rows.append(metadata_row)
            except Exception as e:
                print(f"    ❌ Error creating metadata: {e}")
                return False

        # Dry run summary
        if dry_run:
            print("\n" + "="*60)
            print("🔍 DRY RUN - No files modified")
            print("="*60)
            print(f"\nWould add {len(new_metadata_rows)} strains to Nextstrain:\n")
            for row in new_metadata_rows:
                segments_added = [s for s in ['L', 'M', 'S'] if row[f'segment_{s}'] == '1']
                print(f"  • {row['strain']}")
                print(f"      Segments: {', '.join(segments_added)}")
                print(f"      Collection date: {row['date']}")
            print()
            return True

        # Add to metadata
        try:
            updated_df = pd.concat([nextstrain_df, pd.DataFrame(new_metadata_rows)],
                                  ignore_index=True)
            updated_df.to_csv(existing_metadata, sep='\t', index=False)
            print(f"✅ Updated metadata with {len(new_metadata_rows)} new strains")
        except Exception as e:
            print(f"❌ Error updating metadata: {e}")
            return False

        # Add sequences to FASTA files
        print("\nAdding sequences to FASTA files...")
        try:
            for segment in ['L', 'M', 'S']:
                if sequences_to_add[segment]:
                    fasta_file = self.results_dir / segment / "sequences.fasta"

                    # Ensure directory exists
                    fasta_file.parent.mkdir(parents=True, exist_ok=True)

                    # Check write permissions
                    if fasta_file.exists() and not os.access(fasta_file, os.W_OK):
                        raise PermissionError(f"Cannot write to {fasta_file}")

                    with open(fasta_file, 'a') as f:
                        for strain, sequence in sequences_to_add[segment]:
                            f.write(f">{strain}\n{sequence}\n")

                    print(f"  ✓ Added {len(sequences_to_add[segment])} sequences to {segment} segment")
        except (OSError, PermissionError) as e:
            print(f"❌ Error writing sequence files: {e}")
            return False

        print("\n" + "="*60)
        print("🎉 Successfully imported all assembly sequences!")
        print("="*60)
        return True

    def create_metadata_row(self, user_row, sequence):
        """Create a metadata row in Nextstrain format"""
        # Base metadata that all rows need
        metadata = {
            'strain': user_row['strain'],
            'n_segments': 1,
            'date': user_row['date'],
            'region': user_row.get('region', ''),
            'country': user_row['country'],
            'division': user_row.get('division', ''),
            'location': user_row.get('location', ''),
            'host': user_row.get('host', ''),
            'authors': user_row['authors'],
            'abbr_authors': user_row['authors'],
            'institution': user_row['institution'],
        }

        # Initialize all segment columns as empty
        for seg in ['L', 'M', 'S']:
            metadata.update({
                f'accession_{seg}': '',
                f'accession_version_{seg}': '',
                f'length_{seg}': '',
                f'date_released_{seg}': '',
                f'date_updated_{seg}': '',
                f'sra_accessions_{seg}': '',
                f'segment_{seg}': '0',
                f'qc_{seg}': ''
            })

        # Fill in the specific segment data
        segment = user_row['segment']
        today = datetime.now().strftime('%Y-%m-%d')
        metadata.update({
            f'accession_{segment}': user_row['strain'],
            f'accession_version_{segment}': f"{user_row['strain']}.1",
            f'length_{segment}': len(sequence.seq),
            f'date_released_{segment}': today,
            f'date_updated_{segment}': today,
            f'segment_{segment}': '1',
            f'qc_{segment}': 'good'
        })

        return metadata

    def rebuild_phylogeny(self):
        """Trigger phylogenetic workflow rebuild"""
        print("🔄 Rebuilding phylogenetic workflow...")

        try:
            result = subprocess.run(
                ['nextstrain', 'build', '.', '--forceall'],
                cwd=self.phylo_dir,
                capture_output=True,
                text=True,
                check=True
            )
            print("Phylogenetic rebuild completed successfully!")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Phylogenetic rebuild failed: {e}")
            print(f"Output: {e.stdout}")
            print(f"Error: {e.stderr}")
            return False

def main():
    parser = argparse.ArgumentParser(
        description="Add user sequences to Nextstrain Oropouche workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Import from assembly pipeline (recommended)
  %(prog)s --template my_samples.csv
  %(prog)s --assembly-dir results/assembly_output --metadata my_samples.csv --dry-run
  %(prog)s --assembly-dir results/assembly_output --metadata my_samples.csv --rebuild

  # Manual CSV entry (legacy mode)
  %(prog)s --create-template template.csv
  %(prog)s --input filled_template.csv --validate-only
  %(prog)s --input filled_template.csv --rebuild
        """
    )

    # Mode 1: Assembly import (NEW - recommended)
    assembly_group = parser.add_argument_group('Assembly Pipeline Import (Recommended)')
    assembly_group.add_argument('--template', type=str, metavar='FILE',
                                help='Create simplified metadata CSV template for assembly import')
    assembly_group.add_argument('--assembly-dir', type=str,
                                help='Assembly pipeline output directory')
    assembly_group.add_argument('--metadata', type=str,
                                help='User metadata CSV (strain,accession_L,accession_M,accession_S,date)')

    # Mode 2: Manual CSV (EXISTING - legacy)
    manual_group = parser.add_argument_group('Manual CSV Entry (Legacy)')
    manual_group.add_argument('--input', type=str,
                              help='Input CSV file with manual sequence entries')
    manual_group.add_argument('--create-template', type=str,
                              help='Create CSV template file for manual entry')

    # Common options
    common_group = parser.add_argument_group('Common Options')
    common_group.add_argument('--validate-only', action='store_true',
                              help='Only validate, do not process')
    common_group.add_argument('--dry-run', action='store_true',
                              help='Preview changes without modifying files')
    common_group.add_argument('--rebuild', action='store_true',
                              help='Rebuild phylogenetic workflow after adding')
    common_group.add_argument('--base-dir', type=str, default='.',
                              help='Base directory of Nextstrain repository (default: current directory)')

    args = parser.parse_args()

    processor = UserSequenceProcessor(args.base_dir)

    # Template generation for assembly import
    if args.template:
        processor.create_assembly_template(args.template)
        return

    # Mode 1: Assembly import
    if args.assembly_dir or args.metadata:
        if not args.assembly_dir or not args.metadata:
            parser.error("--assembly-dir and --metadata must be used together")

        success = processor.process_assembly_output(
            args.assembly_dir,
            args.metadata,
            dry_run=args.dry_run
        )

        if success and args.rebuild and not args.dry_run:
            processor.rebuild_phylogeny()

        sys.exit(0 if success else 1)

    # Mode 2: Manual CSV
    if args.create_template:
        processor.create_template(args.create_template)
        return

    if args.input:
        if args.validate_only:
            success = processor.validate_csv(args.input)
            sys.exit(0 if success else 1)

        success = processor.process_sequences(args.input, dry_run=args.dry_run)

        if success and args.rebuild and not args.dry_run:
            processor.rebuild_phylogeny()

        sys.exit(0 if success else 1)

    # No arguments provided
    parser.print_help()

if __name__ == "__main__":
    main()