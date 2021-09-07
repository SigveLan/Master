from Bio import SeqIO
import pandas as pd
from genetic_code import DNA_Codons
import time
from operator import itemgetter
import itertools

# This script takes in the SNPs in coding region results from 'filter.py' and checks if there is a change in AA sequence


def read_exons_to_df(file):

    # Same as the one from filter, but can probably be trimmed down
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


def transcript_input_to_list(transcript_str):
    return transcript_str[2:-2].split(';')


def get_exon_by_id(exon_id):
    return exons[(exons['exon_id'] == exon_id)]


def get_exons_by_transcript_id(transcript_id):
    return exons[(exons['transcript'])]


def get_exons_by_gene_id(gene_id):
    return exons[(exons['gene_id'] == gene_id)]


def translate_DNA(seq):
    protein_sequence = []
    for i in range(0, len(seq) - 2, 3):
        protein_sequence.append(DNA_Codons[seq[i:i + 3]])

    return ''.join(protein_sequence)


def translate_codon(codons):
    return [DNA_Codons[codons[0]], DNA_Codons[codons[1]]]


def compare_and_score_AA(AA):

    if AA[0] == AA[1]:
        return 0
    return 1


def identify_codon(pos, start_phase, seq, var):
    var = var.split('/')
    # Only checks for SNPs that are in codons not split between exons.
    if (1 < pos) & (pos < len(seq) - 2):

        if start_phase != '':
            offset = abs(int(start_phase) - 3)
            pos = pos - offset
            seq = seq[offset:]

        var_seq = seq[0:pos] + var[1] + seq[pos + 1:]
        return [seq[pos - pos % 3:pos - pos % 3 + 3], var_seq[pos - pos % 3:pos - pos % 3 + 3]]

    return []


def assemble_transcripts(gene_data):

    unique_transcript_ids = set(itertools.chain.from_iterable(gene_data.transcript_ids))

    data_dict = {}
    ind = 0

    for transcript in unique_transcript_ids:
        exons_in_transcript = gene_data[pd.DataFrame(gene_data.transcript_ids.tolist()).isin([transcript]).any(1).values].sort_values(by=['exon_start'])
        sequence = ''.join(exons_in_transcript.sequence.tolist())

        print(exons_in_transcript.to_string())

        exon_start = list(filter(None, exons_in_transcript.exon_start.tolist()))
        exon_end = list(filter(None, exons_in_transcript.exon_end.tolist()))

        if not list(filter(None, exons_in_transcript.coding_start.tolist())):
            # Skips non coding transcripts
            continue

        relative_pos_list = [transcript, [], []]
        relative_pos = 0

        for start, stop in zip(exon_start, exon_end):
            relative_pos = relative_pos + abs(stop-start)
            relative_pos_list[1].append(relative_pos)
            relative_pos_list[2].append(stop)

        relative_pos_list.append(sequence)

        # Test range
        print('---------------------------------------')
        print(sequence[relative_pos_list[1][0] : relative_pos_list[1][0] + 3])
        print('---------------------------------------')

        # Test range end

        data_dict[ind] = relative_pos_list
        ind += 1



    return pd.DataFrame.from_dict(data_dict, orient='index', columns=["transcript_id", "relative_position", "global_position", "sequence"])



model_file = 'C:/Users/Sigve/Genome_Data/exon_model_data/exons_chrom_1.fa'
exons = read_exons_to_df(model_file)

results_file = 'C:/Users/Sigve/Genome_Data/results/SNPs_coding.tsv'
SNP_data = pd.read_table(results_file, index_col=0)

assembled_transcript_data = pd.DataFrame()
current_gene = str()

ancestral_AA = []
var_AA = []
change_score = []

for index, data in SNP_data.iterrows():

    if data['gene_id'] is not current_gene:
        gene_data = get_exons_by_gene_id(data['gene_id'])
        assembled_transcript_data = assemble_transcripts(gene_data)

    print(assembled_transcript_data.to_string())
    exit()

    exon = get_exon_by_id(data['exon_id'])
    codons = identify_codon(data['chrom_pos']-exon.iloc[0, 7], exon.iloc[0, 9], exon.iloc[0]['sequence'], data['variant_alleles'])

    if not codons:
        ancestral_AA.append('')
        var_AA.append('')
        change_score.append('')
        continue

    translate_result = translate_codon(codons)
    ancestral_AA.append(translate_result[0])
    var_AA.append(translate_result[1])
    change_score.append(compare_and_score_AA(translate_result))



SNP_data['ancestral_AA'] = ancestral_AA
SNP_data['var_AA'] = var_AA
SNP_data['diff_score'] = change_score

SNP_data.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNPs_effect.tsv', sep='\t')