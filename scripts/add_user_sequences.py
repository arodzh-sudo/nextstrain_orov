#!/usr/bin/env python3
"""
Script to add user sequences to Nextstrain Oropouche workflow.
Supports CSV template for easy metadata entry.
"""

import argparse
import pandas as pd
import os
import sys
import shutil
from pathlib import Path
from Bio import SeqIO
from datetime import datetime
import subprocess

class UserSequenceProcessor:
    def __init__(self, base_dir="."):
        self.base_dir = Path(base_dir)
        self.ingest_dir = self.base_dir / "ingest"
        self.phylo_dir = self.base_dir / "phylogenetic"
        self.results_dir = self.ingest_dir / "results"
        self.backup_dir = self.base_dir / "user_data" / "backup"
        
        # Ensure directories exist
        os.makedirs(self.backup_dir, exist_ok=True)
    
    def create_template(self, output_file):
        """Create a CSV template for user input"""
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
    
    def create_backup(self):
        """Create backup of existing files"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_subdir = self.backup_dir / f"backup_{timestamp}"
        os.makedirs(backup_subdir, exist_ok=True)
        
        # Backup metadata
        metadata_file = self.results_dir / "metadata.tsv"
        if metadata_file.exists():
            shutil.copy2(metadata_file, backup_subdir / "metadata.tsv")
        
        # Backup sequence files
        for segment in ['L', 'M', 'S']:
            seq_file = self.results_dir / segment / "sequences.fasta"
            if seq_file.exists():
                os.makedirs(backup_subdir / segment, exist_ok=True)
                shutil.copy2(seq_file, backup_subdir / segment / "sequences.fasta")
        
        print(f"Backup created: {backup_subdir}")
        return backup_subdir
    
    def process_sequences(self, csv_file, dry_run=False):
        """Process user sequences and add to Nextstrain data"""
        if not self.validate_csv(csv_file):
            return False
        
        print(f"Processing sequences from {csv_file}...")
        
        # Load user data
        user_df = pd.read_csv(csv_file)
        
        if not dry_run:
            # Create backup
            self.create_backup()
        
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
    
    def create_metadata_row(self, user_row, sequence):
        """Create a metadata row in Nextstrain format"""
        # Base metadata that all rows need
        metadata = {
            'strain': user_row['strain'],
            'n_segments': 1,  # Assuming single segment for now
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
    parser = argparse.ArgumentParser(description="Add user sequences to Nextstrain Oropouche workflow")
    parser.add_argument('--create-template', type=str, help='Create CSV template file')
    parser.add_argument('--input', type=str, help='Input CSV file with user sequences')
    parser.add_argument('--validate-only', action='store_true', help='Only validate CSV file, do not process')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without making changes')
    parser.add_argument('--rebuild', action='store_true', help='Rebuild phylogenetic workflow after adding sequences')
    parser.add_argument('--base-dir', type=str, default='.', help='Base directory of Nextstrain repository')
    
    args = parser.parse_args()
    
    processor = UserSequenceProcessor(args.base_dir)
    
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
    
    parser.print_help()

if __name__ == "__main__":
    main()