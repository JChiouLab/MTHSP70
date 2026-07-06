#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 13:38:32 2024

@author: user
"""

#%% load functions
import sys
import os
import re
from Bio import SeqIO
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess
from collections import defaultdict

# Ajouter le répertoire du fichier macse_utils.py au chemin d'importation
sys.path.append('/Users/user/Documents/Python_Scripts/')


#from PAML_plant import get_CDS_from_MAS_AA 
#I have to save def in one script that contain only def


def rename_genes_in_fasta(input_fasta_file, output_fasta_file):
    """
    Reads a FASTA file, modifies the gene names according to specific rules,
    and saves the modified sequences to a new FASTA file. Returns the list of modified gene names.
    
    Parameters:
    - input_fasta_file (str): Path to the input FASTA file.
    - output_fasta_file (str): Path to save the output FASTA file with modified gene names.
    
    Returns:
    - list: A list of modified gene names.
    """
    modified_records = []
    gene_names = []

    # Parse the FASTA file, modify gene names, and store the modified SeqRecords
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        gene_name = record.id  # Extract the gene name from the record ID
        # Apply the transformation rule
        if gene_name.startswith("Zosma") and not gene_name.endswith(".1"):
            gene_name += ".1"
        modified_gene_name = (
            gene_name[:-2] if gene_name.endswith('.p') 
            else re.sub(r'_P(\d{3})$', r'_T\1', re.sub(r'(Achr\w+)P', r'\1T', gene_name))
        )
        gene_names.append(modified_gene_name)  # Add modified name to list

        # Create a new SeqRecord with the modified gene name
        modified_record = SeqRecord(
            record.seq[:-1] if record.seq.endswith('*') else record.seq,  # The sequence remains the same
            id=modified_gene_name,  # Use the modified gene name as the ID
            description=""  # Set an empty description or retain the original as needed
        )
        modified_records.append(modified_record)

    # Write the modified records to a new FASTA file
    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "fasta")

    print(f" {len(gene_names)} modified sequences have been saved to {output_fasta_file}")
    return gene_names
def get_CDS_from_AA(genes,orthogroup_fasta):
    # get CDS
    
    # Liste des gènes
    
    # Chemin vers le répertoire contenant les fichiers FASTA
    CDS_dir = "/RAID1/working/R402/hsushen/PlantProject/genome/CDS/"
    
    # Dictionnaire pour stocker les séquences des gènes
    gene_sequences = {}
    all_gene_id = []
    #print(genes)
    # Parcourir chaque fichier dans le répertoire CDS
    for filename in os.listdir(CDS_dir):
        if filename.endswith(".cds.fa"):  # Vérifier que le fichier est un fichier FASTA
            fasta_file = os.path.join(CDS_dir, filename)
            db_name = filename[:-7]
            db_name = db_name.replace('.', '_')
            # Lire les séquences du fichier FASTA
            for record in SeqIO.parse(fasta_file, "fasta"):
                gene_id = record.id.split()[0]  # Extraire l'identifiant du gène
                #print(gene_id)
                gene_id_sc = gene_id[:-5] #assuming that only Scer has -mrna at the end
            
                #no need to define full name with the proteomes db file names when genes are from OG_seq
                full_name   = gene_id #f'{db_name}_protein_{gene_id}'
                full_name_sc =  gene_id_sc #f'{db_name}_{gene_id_sc}'
                
                # passed_seq_id controls whether the full_name or just the gene_id will be passed into the created .cds.fasta file
                passed_seq_id = full_name
                passed_seq_id_sc = full_name_sc
                seq = str(record.seq)#[:-3]) 
                #Supprimer le codon d'arrêt (les 3 derniers nucléotides) de chaque séquence
                if seq.endswith('TAG') or seq.endswith('TAA') or seq.endswith('TGA'):
                    seq = seq[:-3]
                    
                if (full_name in genes or full_name_sc in genes ) and (passed_seq_id not in gene_sequences or passed_seq_id_sc not in gene_sequences):
                  
                    if full_name_sc in genes: 
                        passed_seq_id = passed_seq_id_sc
                        
                        
                    gene_sequences[passed_seq_id] = seq
                    all_gene_id.append(passed_seq_id)
    with open(orthogroup_fasta, "w") as output_file:
        for gene_id, sequence in gene_sequences.items():
            output_file.write(f">{gene_id}\n{sequence}\n")
            
    
    print(f"CDS seq Extraction terminée. {len(gene_sequences)} séquences ont été extraites.")
    return all_gene_id
def insert_gaps_into_cds(aa_aligned, gene_id, cds_sequence):
    """
    Insert gaps from an amino acid alignment sequence into the CDS sequence.
    Skips sequences with length mismatches.

    Parameters:
    - aa_aligned: str, aligned amino acid sequence with gaps ('-')
    - gene_id: str, identifier for the gene
    - cds_sequence: str, the CDS sequence (un-gapped)

    Returns:
    - str, CDS sequence with gaps inserted if matched; None if mismatched.
    - bool, True if mismatched, False otherwise.
    """
    # Count non-gap characters in aa_aligned
    aa_length = len(aa_aligned.replace('-', ''))

    # Check if the CDS sequence is exactly three times the length of the aa sequence
    if len(cds_sequence) != aa_length * 3:
        return None, True  # Return None and indicate a mismatch
    
    cds_index = 0
    cds_with_gaps = ""
    
    # Insert gaps from the aa_aligned sequence into the CDS
    for aa in aa_aligned:
        if aa == '-':
            cds_with_gaps += '---'  # Add a gap of three nucleotides for a codon
        else:
            # Add the next codon from the CDS sequence
            cds_with_gaps += cds_sequence[cds_index:cds_index + 3]
            cds_index += 3

    return cds_with_gaps, False  # Return gapped CDS and indicate no mismatch
def process_alignment_file(alignment_file, cds_file, output_file):
    """
    Process an alignment file and corresponding CDS file to insert gaps.
    Reports mismatches with specific differences.

    Parameters:
    - alignment_file: str, path to the file containing amino acid alignments in fasta format
    - cds_file: str, path to the file containing CDS sequences in fasta format
    - output_file: str, path for saving the output sequences with gaps
    """
    # Load amino acid sequences into a dictionary keyed by sequence ID
    aa_sequences = {record.id: str(record.seq) for record in SeqIO.parse(alignment_file, "fasta")}
    mismatched_genes = []  # List to store mismatched gene details

    with open(output_file, "w") as out_handle:
        # Iterate over each CDS sequence
        for cds_record in SeqIO.parse(cds_file, "fasta"):
            # Find the corresponding amino acid sequence by ID
            aa_aligned = aa_sequences.get(cds_record.id)
            if aa_aligned:
                # Insert gaps into the CDS sequence
                cds_with_gaps, is_mismatched = insert_gaps_into_cds(aa_aligned, cds_record.id, str(cds_record.seq))
                
                if is_mismatched:
                    # Calculate expected vs actual length and store the discrepancy
                    expected_length = len(aa_aligned.replace('-', '')) * 3
                    actual_length = len(cds_record.seq)
                    mismatched_genes.append({
                        "gene_id": cds_record.id,
                        "expected_length": expected_length,
                        "actual_length": actual_length,
                        "discrepancy": actual_length - expected_length
                    })
                    continue  # Skip this sequence
                
                # Write the output sequence with gaps
                out_handle.write(f">{cds_record.id}\n{cds_with_gaps}\n")
            else:
                print(f"Warning: No corresponding amino acid sequence found for CDS ID {cds_record.id}")

    # Print report of mismatched sequences with specific differences
    if mismatched_genes:
        print("Mismatched CDS sequences with length discrepancies:")
        for gene in mismatched_genes:
            print(f"- Gene ID: {gene['gene_id']}")
            print(f"  Expected Length: {gene['expected_length']} bp")
            print(f"  Actual Length: {gene['actual_length']} bp")
            print(f"  Discrepancy: {gene['discrepancy']} bp")
def get_genes(with_p):

    newick_file = "/RAID1/working/R402/hsushen/PlantProject/genome/proteomes/OrthoFinder/Results_Sep30/WorkingDirectory/OrthoFinder/Results_Nov06/Gene_Trees/OG0002368_tree.txt"
    tree = Phylo.read(newick_file, 'newick')
   
    # Extract terminal/leaf names
    genes = [leaf.name for leaf in tree.get_terminals()]
    # Modify the list to remove '.p' where present
    
    if with_p:
        genes = [gene + ".1" if gene.startswith("Zmarina") and not gene.endswith(".1") else gene for gene in genes]
        genes = [gene[:-2] if gene.endswith('.p') else re.sub(r'_P(\d{3})$', r'_T\1',re.sub(r'(Achr\w+)P', r'\1T', gene))for gene in genes]
    
    #genes = [gene[:-2] if (gene.endswith('.p') and with_p)else re.sub(r'_P(\d{3})$', r'_T\1',re.sub(r'Achr(\d{1})P', r'Achr\1T', gene))for gene in genes]
    #genes = [re.sub(r'_P(\d{3})$', r'_T\1', gene) for gene in genes]
    return genes
def remove_short(genes, OG_aa_fasta,cutoff):
    long_genes = []
    for record in SeqIO.parse(OG_aa_fasta, "fasta"):
        gene_name = record.id
        if len(record.seq)>cutoff:
            for gene in genes:
                if gene.endswith(gene_name):
                    long_genes.append(gene) 
    return long_genes
def remove_multiversion(genes):
 # Filter to keep only the .1 version if available or else the first available version
   #genes = [leaf.name for leaf in tree.get_terminals()]
   
   # Dictionary to store the base name of each gene and its version(s)
   gene_dict = defaultdict(list)
   
   for gene in genes:
       if gene.endswith('.p'):
           gene = gene[:-2]  # Remove the '.p' suffix if it exists
           # Handle base name and version extraction for various formats
       if re.search(r'_T\d+$', gene):  # Handle Zmays-style
           base_name = re.sub(r'_T\d+$', '', gene)  # Remove `_TXXX`
           version = re.search(r'_T(\d+)$', gene)
       elif re.search(r'\.t\d+$', gene):  # Handle Carabicageisha-style
           base_name = re.sub(r'\.t\d+$', '', gene)  # Remove `.tX`
           version = re.search(r'\.t(\d+)$', gene)
       elif re.search(r'\.\d+$', gene):  # Handle general `.X` format (e.g., Taestivum, Tcacao)
           base_name = re.sub(r'\.\d+$', '', gene)  # Remove `.X`
           version = re.search(r'\.(\d+)$', gene)
           
       else:  # If no versioning, use the gene as-is
           
           base_name = gene
           version = None
       # Store the gene under its base name with version number if available
       gene_dict[base_name].append((gene, int(version.group(1)) if version else None))
   
   # Filter to keep only the .1 version if available or else the first available version
   processed_genes = []
   for base_name, versions in gene_dict.items():
       # Sort versions to prioritize .1 if it exists
       versions = sorted(versions, key=lambda x: (x[1] != 1, x[1] if x[1] is not None else float('inf')))
       # Append the full gene name with .1 if it exists, otherwise the first available version
       processed_genes.append(versions[0][0])
   return processed_genes

def filter_genes(genes,OG_fa,filtered_OG_fa,lb):

    
    # Chemin vers le répertoire contenant les fichiers FASTA
    CDS_dir = "/RAID1/working/R402/hsushen/PlantProject/genome/CDS/"
    
    # Dictionnaire pour stocker les séquences des gènes
    gene_sequences = {}
    all_gene_id = []
    #print(genes)
    # Parcourir chaque fichier dans le répertoire CDS
    for filename in os.listdir(CDS_dir):
        if filename.endswith(".cds.fa"):  # Vérifier que le fichier est un fichier FASTA
            
            db_name = filename[:-7]
            db_name = db_name.replace('.', '_')
            # Lire les séquences du fichier FASTA
            for record in SeqIO.parse(OG_fa, "fasta"):
                gene_id = record.id.split()[0]  # Extraire l'identifiant du gène
                #print(gene_id)
                gene_id_sc = gene_id #don't need to remove -mrna at the end again
                #define full name with the proteomes db file names
                full_name   = f'{db_name}_protein_{gene_id}'
                full_name_sc =  f'{db_name}_{gene_id_sc}' 
                
                # passed_seq_id controls whether the full_name or just the gene_id will be passed into the created .cds.fasta file
                passed_seq_id = full_name
                passed_seq_id_sc = full_name_sc
                # find matche or matches in genes
                matches = [gene for gene in genes if passed_seq_id in gene]
                if (full_name in genes or full_name_sc in genes ) and (passed_seq_id not in gene_sequences or passed_seq_id_sc not in gene_sequences):
                # if matches:
                #     for one_m in matches: 
                #         if one_m in genes and (passed_seq_id not in gene_sequences or passed_seq_id_sc not in gene_sequences):
                    #print(full_name)
                    # for match in matches:
                    #     print(match)
                    if full_name_sc in genes: 
                       passed_seq_id = passed_seq_id_sc
                    # Supprimer le codon d'arrêt (les 3 derniers nucléotides) de chaque séquence
                    #gene_sequences[passed_seq_id] = str(record.seq[:-3])
                    gene_sequences[passed_seq_id] = str(record.seq)
                    all_gene_id.append(passed_seq_id)
    with open(filtered_OG_fa, "w") as output_file:
        for gene_id, sequence in gene_sequences.items():
            output_file.write(f">{gene_id}\n{sequence}\n")
            
    
    print(f"{lb} seq Extraction terminée. {len(gene_sequences)} séquences ont été extraites.")
    return all_gene_id
def rename_genes_with_ID(genes,orthogroup_fasta,id_fasta_file,gene_id_mapping):
    # rename genes with ID because alignment cuts off the names
    # Créer un fichier de mappage
    mapping = {}
    for i, gene in enumerate(genes):
        mapping[f"ID{i+1}"] = gene
    
    
    # Lire le fichier .phy et remplacer les noms de gènes par leurs IDs
    with open(orthogroup_fasta, 'r') as infile, open(id_fasta_file, 'w') as outfile:
        for line in infile:
            for id, gene in mapping.items():
                #line = re.sub(r'\b' + re.escape(gene) + r'\b', id, line)
                line = re.sub(r'>' + re.escape(gene) + r'\b', '>' + id , line)
            outfile.write(line)
    
    # Écrire le fichier de mappage
    #gene_id_mapping = "/RAID1/working/R402/hsushen/PlantProject/v2_analysis/id_mapping.txt"
    with open(gene_id_mapping, 'w') as mapfile:
        for id, gene in mapping.items():
            mapfile.write(f"{id}\t{gene}\n")
    return mapping

def insert_gaps_into_cds(aa_aligned, gene_id, cds_sequence):
    """
    Insert gaps from an amino acid alignment sequence into the CDS sequence.
    Skips sequences with length mismatches.

    Parameters:
    - aa_aligned: str, aligned amino acid sequence with gaps ('-')
    - gene_id: str, identifier for the gene
    - cds_sequence: str, the CDS sequence (un-gapped)

    Returns:
    - str, CDS sequence with gaps inserted if matched; None if mismatched.
    - bool, True if mismatched, False otherwise.
    """
    # Count non-gap characters in aa_aligned
    aa_length = len(aa_aligned.replace('-', ''))

    # Check if the CDS sequence is exactly three times the length of the aa sequence
    if len(cds_sequence) != aa_length * 3:
        return None, True  # Return None and indicate a mismatch
    
    cds_index = 0
    cds_with_gaps = ""
    
    # Insert gaps from the aa_aligned sequence into the CDS
    for aa in aa_aligned:
        if aa == '-':
            cds_with_gaps += '---'  # Add a gap of three nucleotides for a codon
        else:
            # Add the next codon from the CDS sequence
            cds_with_gaps += cds_sequence[cds_index:cds_index + 3]
            cds_index += 3

    return cds_with_gaps, False  # Return gapped CDS and indicate no mismatch


def process_alignment_file(alignment_file, cds_file, output_file):
    """
    Process an alignment file and corresponding CDS file to insert gaps.
    Reports mismatches with specific differences.

    Parameters:
    - alignment_file: str, path to the file containing amino acid alignments in fasta format
    - cds_file: str, path to the file containing CDS sequences in fasta format
    - output_file: str, path for saving the output sequences with gaps
    """
    # Load amino acid sequences into a dictionary keyed by sequence ID
    aa_sequences = {record.id: str(record.seq) for record in SeqIO.parse(alignment_file, "fasta")}
    mismatched_genes = []  # List to store mismatched gene details

    with open(output_file, "w") as out_handle:
        # Iterate over each CDS sequence
        for cds_record in SeqIO.parse(cds_file, "fasta"):
            # Find the corresponding amino acid sequence by ID
            aa_aligned = aa_sequences.get(cds_record.id)
            if aa_aligned:
                # Insert gaps into the CDS sequence
                cds_with_gaps, is_mismatched = insert_gaps_into_cds(aa_aligned, cds_record.id, str(cds_record.seq))
                
                if is_mismatched:
                    # Calculate expected vs actual length and store the discrepancy
                    expected_length = len(aa_aligned.replace('-', '')) * 3
                    actual_length = len(cds_record.seq)
                    mismatched_genes.append({
                        "gene_id": cds_record.id,
                        "expected_length": expected_length,
                        "actual_length": actual_length,
                        "discrepancy": actual_length - expected_length
                    })
                    continue  # Skip this sequence
                
                # Write the output sequence with gaps
                out_handle.write(f">{cds_record.id}\n{cds_with_gaps}\n")
            else:
                print(f"Warning: No corresponding amino acid sequence found for CDS ID {cds_record.id}")

    # Print report of mismatched sequences with specific differences
    if mismatched_genes:
        print("Mismatched CDS sequences with length discrepancies:")
        for gene in mismatched_genes:
            print(f"- Gene ID: {gene['gene_id']}")
            print(f"  Expected Length: {gene['expected_length']} bp")
            print(f"  Actual Length: {gene['actual_length']} bp")
            print(f"  Discrepancy: {gene['discrepancy']} bp")
            
def correct_phylip_format(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()
        for line in lines[1:]:
            # Utiliser une expression régulière pour extraire les valeurs de ntax et nchar
            match = re.search(r'dimensions\s+ntax=(\d+)\s+nchar=(\d+);', line)
            if match:
                ntax = match.group(1)
                nchar = match.group(2)
                outfile.write(f"{ntax} {nchar}\n")
       
        header = lines[0].strip()
        #outfile.write(header + '\n')
        seq_data = {}
        current_id = None
        
        for line in lines[5:]:

            if line.strip() == ';' :
               print(line.strip())
               break
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) > 1:
                # New sequence entry
                current_id = parts[0]
                #
                if current_id not in seq_data:
                    seq_data[current_id] = ''.join(parts[1:])
                else:
                    seq_data[current_id] += ''.join(parts[1:])
            else:
                #outfile.write(line)
                
                if current_id:
                    seq_data[current_id] += line.strip()

        for seq_id, sequence in seq_data.items():
            outfile.write(f"{seq_id:<10}  {sequence}\n")
            
def remove_sequences_with_n(input_fasta, output_fasta):
    """
    Reads a FASTA file and writes sequences without 'N' to a new FASTA file.
    
    :param input_fasta: Path to the input FASTA file.
    :param output_fasta: Path to the output FASTA file with filtered sequences.
    """
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if 'N' not in record.seq.upper():
                SeqIO.write(record, output_handle, "fasta")
#%% load OG0000105 aa seq
working_dir = "/RAID1/working/R402/hsushen/PlantProject/v2_analysis"
OG_fasta = "/RAID1/working/R402/hsushen/PlantProject/genome/proteomes/OrthoFinder/Results_Sep30/WorkingDirectory/OrthoFinder/Results_Nov06/Orthogroup_Sequences/OG0002368.fa"
OG_aa_fasta = f"{working_dir}/OG0002368_modified.fa"
gene_names = rename_genes_in_fasta(OG_fasta, OG_aa_fasta)

    
#%% get CDS according to OG0002368 aa seq

# Créer le répertoire s'il n'existe pas
if not os.path.exists(working_dir):
    os.makedirs(working_dir)
os.chdir(working_dir)

OG_cds_fasta = f"{working_dir}/OG0002368_cds.fa"
all_gene_ids = get_CDS_from_AA(gene_names, OG_cds_fasta)
#%% compare length 
tbs_file = f"{working_dir}/OG0002368_cds_tbs.fa"             
process_alignment_file(OG_aa_fasta, OG_cds_fasta, tbs_file)
#%% get genes from the tree
genes_filter = get_genes(1)
#%% remove short genes
cutoff = 400
genes_filter = remove_short(genes_filter,OG_aa_fasta,cutoff)
print('Short genes removed')
#%% remove multiple versions of the same gene and keep only .1
genes_filter = remove_multiversion(genes_filter)
print('Extra versions of the same gene removed')
#%% filter OG .fa
filtered_aa = f"{working_dir}/mtHSP70_aa_tbs.fa"
lb = "AA"  
genes_f_aa = filter_genes(genes_filter, OG_aa_fasta, filtered_aa, lb)
filtered_cds = f"{working_dir}/mtHSP70_cds_tbs.fa" 
lb = "CDS" 
genes_f_cds = filter_genes(genes_filter, OG_cds_fasta, filtered_cds, lb)
#%% replace headings with ID
gene_id_mapping = f"{working_dir}/id_mapping.txt"

input_aa= filtered_aa   
id_aa = f"{working_dir}/mtHSP70_aa_id.fa"
mapping = rename_genes_with_ID(genes_filter, input_aa, id_aa, gene_id_mapping)

input_cds= filtered_cds   
id_cds = f"{working_dir}/mtHSP70_cds_id.fa"
mapping = rename_genes_with_ID(genes_filter, input_cds, id_cds, gene_id_mapping)
#alignment = AlignIO.read(id_fasta, "fasta")
#%% make MSA with filtered AA


# Define paths
fasta_file = id_aa
aligned_aa = f"{working_dir}/mtHSP70_aa.fa"

# Run MUSCLE for MSA with MUSCLE 5.0 syntax
try:
    # Use subprocess to execute the MUSCLE command with specific flags for MUSCLE 5.0
    muscle_cmd = [
        "muscle", "-align", fasta_file, "-output", aligned_aa
    ]
    subprocess.run(muscle_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"Error occurred: {e}")
    print(f"Standard Output: {e.stdout}")
    print(f"Standard Error: {e.stderr}")

# Read and display the alignment if successful
try:
    alignment = AlignIO.read(aligned_aa, "fasta")
    print(alignment)
except Exception as e:
    print(f"Error reading alignment: {e}")

#%% insert gap into cds seq 
gap_fa = f"{working_dir}/mtHSP70_cds_id_gap.fa"             
process_alignment_file(aligned_aa , id_cds, gap_fa)

# filter out seq with N
noN_fa = f"{working_dir}/mtHSP70_cds_id_gap_noN.fa"  
remove_sequences_with_n(gap_fa, noN_fa)

#%% Correct to nuc file
alignment = AlignIO.read(noN_fa, "fasta")

 # Définir le type de molécule (ADN dans ce cas)
for record in alignment:
    record.annotations["molecule_type"] = "DNA"

nexus_file = f"{working_dir}/mtHSP70_cds_id_gap.nuc"
AlignIO.write(alignment, nexus_file, "nexus")
print("Alignement sauvegardé en format NEXUS")

input_file = nexus_file#f"/Users/user/Documents/Alignments/{OG}_cds_aligned.nuc"
#output_file = f"/Users/user/Documents/Alignments/{OG}_cds_aligned_ns_corrected.nuc"
output_file = f"{working_dir}/mtHSP70_cds_id_corrected.nuc"
correct_phylip_format(input_file, output_file)


#%%  Commande pour IQ-TREE
    
phylip_file = output_file
iqtree_command = f"iqtree2 --redo -s {phylip_file} -m MFP -B 1000 -alrt 1000 -st CODON -nt AUTO " #"--redo  --model-joint UNREST
subprocess.run(iqtree_command, shell=True)
#%% change ID back to gene names
    
tree_file = f"{working_dir}/mtHSP70_cds_id_corrected.nuc.treefile"

# # Lire le fichier Newick

# with open(tree_file, "r") as file:
#     newick_content = file.read()

# replace leaf name ID with gene names
# Lire l'arbre IQtree et remplacer les IDs par les noms de gènes
iqtree_output_file = tree_file #
#gene_id_mapping = "{working_dir}/id_mapping.txt"
output_tree_file =  f"{working_dir}/mtHSP70_cds_gn_corrected.nuc.treefile"

# Create a dictionary to store the gene ID mappings
gene_id_map = {}

# Read the gene ID mapping file and populate the dictionary
with open(gene_id_mapping, 'r') as mapfile:
    for line in mapfile:
        # Assuming the mapping file has two columns: ID and gene name, separated by a tab or space
        id, gene = line.strip().split()  # Adjust the split character if needed
        gene_id_map[id] = gene

# Now read the IQtree output file and replace the IDs with gene names
with open(iqtree_output_file, 'r') as infile, open(output_tree_file, 'w') as outfile:
    for line in infile:
        for id, gene in gene_id_map.items():
            line = re.sub(r'\b' + re.escape(id) + r'\b', gene, line)  # Replace IDs with gene names
        outfile.write(line)



print("Les IDs dans le fichier de l'arbre IQtree ont été remplacés par les noms de gènes.")


#%% PAML


 
# # Définir les chemins des fichiers
# alignment_file = output_file#f"/Users/user/Documents/Alignments/{OG}_cds_aligned_corrected.nuc"
# #alignment_file = f"/Users/user/Documents/Alignments/{OG}_cds_aligned.nuc"
# #tree_file = f"/Users/user/Documents/Alignments/{OG}_cds_aligned_ns.phy.treefile"
# tree_file = f"{working_dir}/mtHSP70_cds_id_corrected.nuc.treefile"

# # Définir le répertoire de travail

# #model = '2'
# #NSsites = '0 1 2 7 8'
# # Chemin du fichier codeml.ctl
# ctl_file_path = os.path.join(working_dir, "codeml.ctl")
# result_path = os.path.join(working_dir, "results.txt")
# # Créer un fichier de configuration codeml.ctl
# with open(ctl_file_path, "w") as ctl_file:
#     ctl_file.write(f"""
#       seqfile = {alignment_file}            * Path to the alignment file
#      treefile = {tree_file}           * Path to the tree file
#       outfile = {result_path}            * Path to the output file
   
#         noisy = 3              * How much rubbish on the screen
#       verbose = 1              * More or less detailed report

#       seqtype = 1              * Data type
#         ndata = 1           * Number of data sets or loci
#         icode = 0              * Genetic code 
#     cleandata = 0              * Remove sites with ambiguity data?
# 		
#         model = 2      * Models for ω varying across lineages
# 	  NSsites = 2         * Models for ω varying across sites
#     CodonFreq = 7        * Codon frequencies
# 	  estFreq = 0        * Use observed freqs or estimate freqs by ML
#         clock = 0          * Clock model
#     fix_omega = 0         * Estimate or fix omega
#         omega = 0.5        * Initial or fixed omega
#  RateAncestor = 1   * (0,1,2): rates (alpha>0) or ancestral states
#     """)
    
# paml_command = f"codeml {ctl_file_path}"

# subprocess.run(paml_command, shell=True)

