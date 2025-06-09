# Create template
python scripts/add_user_sequences.py --create-template user_sequences.csv

# Validate user data
python scripts/add_user_sequences.py --input user_sequences.csv --validate-only

# Dry run (preview changes)
python scripts/add_user_sequences.py --input user_sequences.csv --dry-run

# Process sequences and rebuild phylogeny
python scripts/add_user_sequences.py --input user_sequences.csv --rebuild