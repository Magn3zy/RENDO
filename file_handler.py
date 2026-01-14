import os
from pathlib import Path
import shutil
from datetime import datetime

BASE_DIR = os.getcwd()
TEMP_DIR = os.path.join(BASE_DIR, "temp_files", datetime.now().strftime("%Y%m%d_%H%M%S"))
FASTA_OUTPUT_DIR = os.path.join(BASE_DIR, "fasta_outputs")
LOG_DIR = os.path.join(BASE_DIR, "logs")

def initialize_directories():
    for directory in [TEMP_DIR, FASTA_OUTPUT_DIR, LOG_DIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)

def clean_temp_dir():
    if os.path.exists(TEMP_DIR):
        shutil.rmtree(TEMP_DIR)
    os.makedirs(TEMP_DIR)

def validate_sequence(sequence):
    valid_bases = set('ATCGN')
    invalid_chars = set(char.upper() for char in sequence) - valid_bases
    if invalid_chars:
        raise ValueError(f"Invalid characters found in sequence: {invalid_chars}")
    return True

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement_dict.get(base.upper(), base) for base in reversed(sequence))

def load_fasta_with_orientation(file_path, orientation):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    if orientation not in ['3_5', '5_3']:
        raise ValueError("Orientation must be either '3_5' or '5_3'")

    sequences = {}
    file_name = os.path.basename(file_path)
    name, ext = os.path.splitext(file_name)

    temp_file_name = f"{name}_{orientation}_{datetime.now().strftime('%Y%m%d_%H%M%S')}{ext}"
    temp_file_path = os.path.join("temp_files", temp_file_name)

    try:
        with open(file_path, "r") as file:
            seq_name = None
            seq_data = []
            for line in file:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    if seq_name and seq_data:
                        sequence = ''.join(seq_data)
                        sequences[seq_name] = sequence
                    seq_name = line[1:]
                    seq_data = []
                else:
                    seq_data.append(line.upper())

            if seq_name and seq_data:
                sequence = ''.join(seq_data)
                sequences[seq_name] = sequence

        if orientation == "3_5":
            sequences = {name: reverse_complement(seq) for name, seq in sequences.items()}

        with open(temp_file_path, "w") as temp_file:
            for name, seq in sequences.items():
                temp_file.write(f">{name}\n{seq}\n")

        return sequences, temp_file_path

    except Exception as e:
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)
        raise e

if __name__ == "__main__":
    try:
        file_path = "path_to_your_fasta_file.fasta"
        orientation = "5_3"  # or "3_5"
        sequences, temp_file_path = load_fasta_with_orientation(file_path, orientation)
        print(f"Sequences: {sequences}")
        print(f"Temporary file path: {temp_file_path}")
    except Exception as e:
        print(f"Error: {e}")