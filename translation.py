from Bio import SeqIO
import pandas as pd
import numpy as np
from genetic_code import DNA_Codons
import time
from operator import itemgetter
import itertools
import ast

# This script takes in the SNPs in coding region results from 'filter.py' and checks if there is a change in AA sequence


def read_exons_to_df(file):

    # Function is the same as the one from filter
    # Can probably be trimmed down

    # Reading to dictionary then convert to dataframe is very fast.
    dictionary = {}
    ind = 0

    for seq_record in SeqIO.parse(open(file, mode='r'), 'fasta'):

        seq_record_attributes = seq_record.description.split("|")
        seq_record_attributes[5] = seq_record_attributes[5].split(';')

        # These are always present as numbers.
        for i in [3, 4, 7, 8, 13]:
            seq_record_attributes[i] = int(seq_record_attributes[i])

        # The phases can be empty, which means the same as phase zero.
        for i in [9, 10]:
            if seq_record_attributes[i] != '':
                seq_record_attributes[i] = int(seq_record_attributes[i])
            else:
                seq_record_attributes[i] = 0

        # Although rare, there can be more than one coding start/stop in a given exon.
        # This picks the closest to the edge of the exon.
        for i in [11, 12]:
            if seq_record_attributes[i] != '':
                pos = list(map(int, seq_record_attributes[i].split(';')))
                if len(pos) > 1:
                    differences = [abs(j - seq_record_attributes[i-4]) for j in pos]
                    seq_record_attributes[i] = pos[min(enumerate(differences), key=itemgetter(1))[0]]

                else:
                    seq_record_attributes[i] = int(seq_record_attributes[i])

        # Adds the sequence
        seq_record_attributes.append(str(seq_record.seq))

        dictionary[ind] = seq_record_attributes
        ind += 1

    return create_df(dictionary)


def create_df(dictionary):

    df = pd.DataFrame.from_dict(dictionary, orient='index', columns=["gene_id", "gene_name", "chrom", "gene_start",
                                                                     "gene_end", "transcript_ids", "exon_id",
                                                                     "exon_start", "exon_end", "phase_start",
                                                                     "phase_end", "coding_start", "coding_end",
                                                                     "strand", "sequence"])
    return df


def get_exons_by_gene_id(gene_id):
    return exons[(exons['gene_id'] == gene_id)]


def translate_DNA(seq):
    protein_sequence = []
    for i in range(0, len(seq) - 2, 3):
        protein_sequence.append(DNA_Codons[seq[i:i + 3]])

    return ''.join(protein_sequence)


def compare_and_score_AA(AA):

    if AA[0] == AA[1]:
        return 0
    return 1


def assemble_transcripts(gene_data):

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

model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_1.fa'
exons = read_exons_to_df(model_file)

filtered_SNPs_file = 'C:/Users/Sigve/Genome_Data/results/SNPs_coding.tsv'
SNP_data = pd.read_table(filtered_SNPs_file, index_col=0)

assembled_transcript_data = pd.DataFrame()
current_gene = str()
current_strand = 1

results = []

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

        results.append(data_as_list[:6] + [transcript] + data_as_list[7:] + result)


columns = SNP_data.columns.values.tolist() + ['amino_acid_change', 'amino_acid_pos', 'score']
columns[6] = 'transcript_id'

results_as_df = pd.DataFrame(results, columns=columns)
results_as_df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNPs_effect.tsv', sep='\t')

end_time = time.time()
print('SNP translation execution time: %.6f seconds' % (end_time-start_time))

