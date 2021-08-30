from Bio import SeqIO
import pandas as pd
from genetic_code import DNA_Codons



def translate_DNA(seq, init_pos = 0):
    protein_sequence = []
    seq = seq[init_pos:]
    for i in range(0, len(seq) - 2, 3):
        protein_sequence.append(DNA_Codons[seq[i:i + 3]])

    return ''.join(protein_sequence)

model_file = 'C:/Users/sigve/Documents/Genome_Data/exon_model_data/exons_chrom_1.fa'
exons = read_exons_to_df(model_file)

t_exons = exons[pd.DataFrame(exons['transcript_ids'].tolist()).isin([t_id]).values]
            print(len(t_exons))
            for j, row2 in t_exons.iterrows():
                t_seq[len(t_ids)-1] += row2['sequence']