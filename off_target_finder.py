import time
from Bio import SeqIO
import gzip
import tempfile
import shutil
import os
import subprocess

# Configuration constants
DEFAULT_MAX_OFF_TARGETS = 4
DEFAULT_PAM_SEQUENCE = "GG"
GENOME_TYPE_FLAG = "G"
CLEANUP_DELAY_SECONDS = 2


def make_temp_file(compressed_file_path):
    with tempfile.NamedTemporaryFile(delete=False, mode='wb', suffix=".fasta") as temp_decompressed_file:
        with gzip.open(compressed_file_path, 'rb') as compressed_input:
            shutil.copyfileobj(compressed_input, temp_decompressed_file)
        decompressed_file_path = temp_decompressed_file.name
        return decompressed_file_path


def prepare_txt_file(genome_file_path, sgrna_input_file_path, max_off_targets, preparation_output_file_path,
                     pam_sequence=DEFAULT_PAM_SEQUENCE):
    temporary_genome_file_path = make_temp_file(genome_file_path)
    sgrna_sequences_content = ""
    pam_pattern = ""
    max_sequence_length = 0

    with open(sgrna_input_file_path, "r") as sgrna_file_handle:
        for i, sequence_record in enumerate(SeqIO.parse(sgrna_file_handle, "fasta")):
            max_sequence_length = len(sequence_record.seq) if max_sequence_length < len(
                sequence_record.seq) else max_sequence_length
            sgrna_sequences_content += f"{sequence_record.seq}\t{max_off_targets}\n"
        pam_pattern = f"N" * (max_sequence_length - len(pam_sequence)) + pam_sequence + "\n"

    with open(preparation_output_file_path, "w") as output_file:
        output_file.write(temporary_genome_file_path + "\n" + pam_pattern + sgrna_sequences_content)
    return temporary_genome_file_path


def find_off_targets(cas_offinder_path, genome_file_path, sgrna_input_file_path, preparation_output_file_path,
                     cas_output_file_path, max_off_targets=DEFAULT_MAX_OFF_TARGETS, pam_sequence=DEFAULT_PAM_SEQUENCE):
    temporary_genome_path = prepare_txt_file(
        genome_file_path,
        sgrna_input_file_path,
        max_off_targets,
        preparation_output_file_path,
        pam_sequence
    )

    subprocess.run([cas_offinder_path, preparation_output_file_path, GENOME_TYPE_FLAG, cas_output_file_path])
    time.sleep(CLEANUP_DELAY_SECONDS)
    os.remove(temporary_genome_path)

    return cas_output_file_path


def run_cas_offinder(cas_offinder_executable, genome_file, sgrna_file, output_dir, gene,conserved_region,max_off_targets=DEFAULT_MAX_OFF_TARGETS):

    cas_output = f"{output_dir}/{gene}_{conserved_region}_offinder_results.txt"
    preparation_output = f"{output_dir}/{gene}_{conserved_region}offinder_preperation.txt"
    result_file = find_off_targets(
        cas_offinder_executable,
        genome_file,
        sgrna_file,
        preparation_output,
        cas_output,
        max_off_targets
    )
    return result_file


