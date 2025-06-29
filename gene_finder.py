#!/usr/bin/env python3
"""
CESA Gene Analysis Pipeline

This script analyzes CESA (Cellulose Synthase) genes in Salvia rosmarinus by:
1. Identifying CESA-related genes from annotations
2. Extracting CDS features from GFF file
3. Loading genomic sequences
4. Building CDS sequences
5. Translating to protein sequences

Dependencies:
- Biopython
- gzip
"""

import csv
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from typing import Dict, List, Tuple, Set
import logging

############################################################################
# Configuration
############################################################################

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


############################################################################
# Annotation Processing Functions
############################################################################

def extract_cesa_gene_ids(annotation_file: str, gene_name, gene_subfamily) -> Set[str]:
    """
    Extract gene IDs related to CESA from an .annotation file.

    Args:
        annotation_file: Path to an .annotation file

    Returns:
        Set of gene IDs related to CESA
    """
    query_ids = set()

    try:
        with gzip.open(annotation_file, "rt") as f:
            # Process header lines and skip comments
            lines = [
                line[1:].strip() if line.startswith("#") and not line.startswith("##")
                else line.strip()
                for line in f if not line.startswith("##")
            ]

        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            if (gene_name.lower() in row["Preferred_name"].lower() or
                    gene_subfamily.lower() in row["Description"].lower()):
                logging.info(f"Found CESA gene: {row['query']} | {row['Preferred_name']}")
                query_ids.add(row["query"])

    except FileNotFoundError:
        logging.error(f"Annotation file not found: {annotation_file}")
    except Exception as e:
        logging.error(f"Error processing annotation file: {e}")

    return query_ids


############################################################################
# GFF Feature Processing
############################################################################

def parse_gff_cds_features(gff_file: str, query_ids: Set[str]) -> Dict[str, List[Tuple]]:
    """
    Parse CDS features from GFF file for specified gene IDs.

    Args:
        gff_file: Path to GFF file
        query_ids: Set of gene IDs to extract CDS features for

    Returns:
        Dictionary mapping gene IDs to lists of CDS coordinates (scaffold, start, end, strand)
    """
    cds_coords = {}

    try:
        with gzip.open(gff_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 9 or parts[2] != "CDS":
                    continue

                # Extract feature information
                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attributes = dict(
                    attr.split("=", 1)
                    for attr in parts[8].split(";")
                    if "=" in attr
                )

                # Get parent ID and store coordinates if matching query
                parent_id = attributes.get("Parent") or attributes.get("ID")
                if parent_id in query_ids:
                    cds_coords.setdefault(parent_id, []).append(
                        (scaffold, start, end, strand)
                    )

    except FileNotFoundError:
        logging.error(f"GFF file not found: {gff_file}")
    except Exception as e:
        logging.error(f"Error processing GFF file: {e}")

    return cds_coords


############################################################################
# Sequence Processing Functions
############################################################################

def load_genome_sequences(genome_file: str) -> Dict[str, Seq]:
    """
    Load genome sequences from FASTA file.

    Args:
        genome_file: Path to genome FASTA file (gzipped)

    Returns:
        Dictionary mapping scaffold IDs to sequences
    """
    genome_dict = {}

    try:
        with gzip.open(genome_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                clean_id = record.id.split()[0].split(".")[0]
                genome_dict[clean_id] = record.seq

        logging.info(f"Loaded {len(genome_dict)} scaffold sequences")

    except FileNotFoundError:
        logging.error(f"Genome file not found: {genome_file}")
    except Exception as e:
        logging.error(f"Error loading genome sequences: {e}")

    return genome_dict


def build_cds_sequences(
        cds_coords: Dict[str, List[Tuple]],
        genome_dict: Dict[str, Seq]
) -> str:
    """
    Build CDS sequences from coordinates and genome sequences.

    Args:
        cds_coords: Dictionary of CDS coordinates by gene ID
        genome_dict: Dictionary of genome sequences by scaffold

    Returns:
        FASTA format string of CDS sequences
    """
    fasta_string = ""

    for gene_id, exon_list in cds_coords.items():
        if not exon_list:
            logging.warning(f"No CDS exons found for {gene_id}")
            continue

        # Get scaffold sequence
        scaffold = exon_list[0][0]
        seq = genome_dict.get(scaffold)
        if not seq:
            scaffold_clean = scaffold.split()[0].split(".")[0]
            seq = genome_dict.get(scaffold_clean)

        if not seq:
            logging.error(f"Scaffold '{scaffold}' not found for gene {gene_id}")
            continue

        # Build CDS sequence
        sorted_exons = sorted(exon_list, key=lambda x: x[1])
        cds_sequence = Seq("")
        strand = sorted_exons[0][3]

        for _, start, end, _ in sorted_exons:
            cds_sequence += seq[start - 1:end]

        # Handle reverse strand
        if strand == "-":
            cds_sequence = cds_sequence.reverse_complement()

        fasta_string += f">{gene_id} {scaffold} stitched CDS\n{cds_sequence}\n"

    return fasta_string


############################################################################
# Conserved Sequence Analysis
############################################################################

def extract_conserved_sequence(
        nucleotide_seq: str,
        protein_sequence: str,
        conserved_residues: str
) -> tuple[str, bool]:
    """
    Extract nucleotide sequence from the closest ATG before conserved_residues motif up to and including conserved_residues.

    Args:
        nucleotide_seq: The nucleotide sequence to analyze
        protein_sequence: The translated protein sequence
        conserved_residues: The conserved amino acid sequence to search for

    Returns:
        Tuple containing:
        - sequence: The extracted sequence from closest ATG to conserved_residues (empty if not found)
        - found: Boolean indicating if conserved sequence was found
    """
    # Find position of conserved_residues motif
    motif_pos = protein_sequence.find(conserved_residues)
    if motif_pos == -1:  # conserved_residues not found
        return "", False

    # Convert protein position to nucleotide position (3 nucleotides per amino acid)
    residue_start = motif_pos * 3  # Start of conserved_residues in nucleotides
    residue_end = residue_start + (len(conserved_residues) * 3)  # End of conserved_residues

    # Find all ATG codons in frame and before conserved_residues
    last_valid_atg = -1
    atg_index = nucleotide_seq.find("ATG")

    while atg_index != -1 and atg_index < residue_start:
        if atg_index % 3 == 0:  # Check if in frame
            last_valid_atg = atg_index
        atg_index = nucleotide_seq.find("ATG", atg_index + 3)

    if last_valid_atg != -1:
        # Extract sequence from last valid ATG before conserved_residues to end of conserved_residues
        valid_seq = nucleotide_seq[last_valid_atg:residue_end]
        return valid_seq, True

    return "", False


############################################################################
# Translation Functions
############################################################################

def translate_cds_to_protein(input_fasta: str, cds_out: str, cnsvrd_orf_out: str, conserved_residues: str):
    """
    Translate CDS sequences to protein sequences and identify conserved residue motifs.

    Args:
        input_fasta: Path to input CDS FASTA file
        cds_out: Path to output protein FASTA file
        cnsvrd_orf_out: Path to output conserved ORF FASTA file
        conserved_residues: The conserved amino acid sequence to search for
    """
    try:
        fasta_string = ""
        with open(cds_out, "w") as out_handle:
            for record in SeqIO.parse(input_fasta, "fasta"):
                # Translate and create new record
                protein_seq = record.seq.translate(to_stop=True)
                protein_record = record[:]
                protein_record.seq = protein_seq
                protein_record.description = "translated protein"

                if conserved_residues not in protein_record.seq.upper():
                    continue

                sequence, found = extract_conserved_sequence(
                    str(record.seq),
                    str(protein_record.seq),
                    conserved_residues
                )
                if found:
                    fasta_string += f">{record.description} {conserved_residues} ORF\n{sequence}\n"
                SeqIO.write(protein_record, out_handle, "fasta")

        with open(cnsvrd_orf_out, 'w') as f:
            f.write(fasta_string)

        logging.info(f"Protein sequences saved to {cds_out}")

    except Exception as e:
        logging.error(f"Error translating sequences: {e}")

############################################################################
# Main Execution
############################################################################

def run_analysis(genome, annotation, gff, output_dir, gene, gene_subfamily, conserved_residues="WPGN"):
    """Main execution function."""
    # File paths
    annotation_file = annotation
    gff_file = gff
    genome_file = genome
    cds_output = f"{output_dir}/{gene}_coding_sequences.fasta"
    protein_output = f"{output_dir}/{gene}_protein_sequences.fasta"
    sgrna_orf_out = f"{output_dir}/{gene}_sgrna_region.fasta"

    # Step 1: Extract CESA gene IDs from annotations
    query_ids = extract_cesa_gene_ids(annotation_file, gene, gene_subfamily)

    # Step 2: Parse CDS features from GFF file
    cds_coords = parse_gff_cds_features(gff_file, query_ids)

    # Step 3: Load genome sequences
    genome_dict = load_genome_sequences(genome_file)

    # Step 4: Build CDS sequences
    fasta_string = build_cds_sequences(cds_coords, genome_dict)

    # Step 5: Write CDS sequences
    with open(cds_output, "w") as f:
        f.write(fasta_string)
    logging.info(f"CDS sequences saved to {cds_output}")

    # Step 6: Translate to protein and identify conserved sequences
    translate_cds_to_protein(cds_output, protein_output, sgrna_orf_out, conserved_residues.upper())