import pandas as pd
import numpy as np
from genetic_code import DNA_Codons
import time
import itertools
import ast
from read_exons_to_df import read_exons_to_df
"""
This script takes in the SNPs in coding region results from 'filter.py' and checks if there is a change in AA sequence
"""

def get_exons_by_gene_id(gene_id):
    """Return all the exons for a given gene"""
    return exons[(exons['gene_id'] == gene_id)]


def translate_DNA(seq):
    """Translates DNA, not in use currently"""
    return ''.join([DNA_Codons[seq[i:i + 3]] for i in range(0, len(seq) - 2, 3)])


def compare_and_score_AA(AA):
    """Compares two amino acids and give a score"""
    if AA[0] == AA[1]:
        return 0
    return 1


def assemble_transcripts(gene_data):
    """Assembles the different transcripts for a given gene"""
    unique_transcript_ids = set(itertools.chain.from_iterable(gene_data.transcript_ids))

    data_dict = {}
    ind = 0

    if globals()['current_strand'] == 1:
        strand = True
        coding_pos = "coding_start"
    else:
        strand = False
        coding_pos = "coding_end"

    for transcript in unique_transcript_ids:

        exons_in_transcript = gene_data[
            pd.DataFrame(gene_data.transcript_ids.tolist()).isin([transcript]).any(1).values].sort_values(
            by=['exon_start'], ascending=strand)

        sequence = ''.join(exons_in_transcript.sequence.tolist())

        coding_start = int(exons_in_transcript.loc[exons_in_transcript[coding_pos].ne('').idxmax(), [coding_pos]])

        exon_start = list(filter(None, exons_in_transcript.exon_start.tolist()))
        exon_end = list(filter(None, exons_in_transcript.exon_end.tolist()))

        relative_pos_list = [transcript, coding_start, [], []]
        relative_pos = 0

        for start, stop in zip(exon_start, exon_end):
            relative_pos = relative_pos + abs(stop-start) + 1
            relative_pos_list[2].append(relative_pos)
            relative_pos_list[3].append(stop)

        relative_pos_list.append(sequence)
        relative_pos_list[2] = np.array(relative_pos_list[2])
        relative_pos_list[3] = np.array(relative_pos_list[3])

        data_dict[ind] = relative_pos_list
        ind += 1

    return pd.DataFrame.from_dict(data_dict, orient='index', columns=["transcript_id", "abs_coding_start", "rel_pos",
                                                                      "abs_pos", "sequence"])


def get_rel_pos(atd, pos, strand):
    """Finds the relative position inside a sequence based on an absolute position in the chromosome"""

    abs_exon_pos = atd.iloc[3][atd.iloc[3] >= pos].min()

    if strand == 1:
        rel_pos = int(atd.iloc[2][np.where(atd.iloc[3] == abs_exon_pos)[0][0]] - abs(abs_exon_pos - pos))
    else:
        indices = np.where(atd.iloc[3] > abs_exon_pos)[0]
        if indices.size != 0:
            rel_pos = atd.iloc[2][indices[-1]] + abs(abs_exon_pos - pos)
        else:
            rel_pos = abs(abs_exon_pos - pos)

    return rel_pos


def SNP_eval(atd, pos, var, strand):
    """Evaluates if a SNP produces a change in amino acid sequence for a given transcript"""

    coding_start_pos = get_rel_pos(atd, atd.iloc[1], strand)
    SNP_pos = get_rel_pos(atd, pos, strand) - coding_start_pos

    # print("Exon_end_chosen: " + str(atd.iloc[3][atd.iloc[3] >= atd.iloc[1]].min()))
    # print("Coding_start: " + str(coding_start_pos))

    seq = atd.iloc[4][coding_start_pos-int((1+strand)/2):]
    SNP_codon = seq[SNP_pos - SNP_pos % 3:SNP_pos - SNP_pos % 3 + 3]

    var = var.split('/')
    SNP_codon_var = SNP_codon[:SNP_pos % 3] + var[1] + SNP_codon[(SNP_pos % 3) + 1:]

    amino_acids = [DNA_Codons[SNP_codon], DNA_Codons[SNP_codon_var]]
    diff_score = compare_and_score_AA(amino_acids)

    if not diff_score:
        return 0

    var_seq = seq[:SNP_pos] + var[1] + seq[SNP_pos + 1:]
    var_amino_acid_pos = (SNP_pos // 3) + 1
    len_var_seq = len(var_seq)
    translated_var_seq = []

    for i in range(0, len_var_seq - len_var_seq % 3, 3):
        amino_acid = DNA_Codons[var_seq[i:i+3]]
        translated_var_seq.append(amino_acid)
        if amino_acid == '_':
            break

    # Full sequence if needed
    # translated_var_seq = ''.join(translated_var_seq)[:-1]
    return [amino_acids[0] + "/" + amino_acids[1], var_amino_acid_pos, diff_score]


start_time = time.time()

model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_all_filtered.fa'
exons = read_exons_to_df(model_file)

filtered_SNPs_file = 'C:/Users/Sigve/Genome_Data/results/SNPs_coding.tsv'
SNP_data = pd.read_table(filtered_SNPs_file, index_col=0)

assembled_transcript_data = pd.DataFrame()
current_gene = str()
current_strand = 1

result_list = []

for index, data in SNP_data.iterrows():

    if data['gene_name'] is not current_gene:
        current_gene = data['gene_name']
        gene_data = get_exons_by_gene_id(data['gene_id'])
        current_strand = gene_data.iloc[0, 13]
        assembled_transcript_data = assemble_transcripts(gene_data)

    transcript_ids = ast.literal_eval(data['transcript_ids'])
    data_as_list = data.tolist()

    for transcript in transcript_ids:
        result = SNP_eval(assembled_transcript_data[assembled_transcript_data.transcript_id == transcript].iloc[0], data['chrom_pos'], data['variant_alleles'], current_strand)

        if not result:
            continue

        result_list.append(data_as_list[:6] + [transcript] + data_as_list[7:] + result)


columns = SNP_data.columns.values.tolist() + ['amino_acid_change', 'amino_acid_pos', 'score']
columns[6] = 'transcript_id'

pd.DataFrame(result_list, columns=columns).to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNPs_effect.tsv', sep='\t')

end_time = time.time()
print('SNP evaluation execution time: %.6f seconds' % (end_time-start_time))

