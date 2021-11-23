import pandas as pd
import numpy as np
from src.genetic_code import DNA_Codons, DNA_ReverseComplement
import itertools
from tqdm import tqdm


def SNP_effect_eval(exons: pd.DataFrame, SNP_data: pd.DataFrame) -> pd.DataFrame:
    """This function takes in the SNPs in coding region results from 'SNP_filter.py' and checks if there is a change in
       amino acid sequence."""

    def compare_and_score_AA(AA: list) -> int:
        """Compares two amino acids and give a score."""
        if AA[0] == AA[1]:
            return 0
        return 1

    def assemble_transcripts(gene_data: pd.DataFrame, strand: int) -> pd.DataFrame:
        """Assembles the different transcripts for a given gene by concatenating exons."""
        unique_transcript_ids = set(itertools.chain.from_iterable(gene_data.transcript_ids))

        data_dict = {}
        ind = 0

        if strand == 1:
            strand = True
            coding_pos = "coding_start"
        else:
            strand = False
            coding_pos = "coding_end"

        for transcript_id in unique_transcript_ids:

            exons_in_transcript = gene_data[
                pd.DataFrame(gene_data.transcript_ids.tolist()).isin([transcript_id]).any(1).values].sort_values(
                by=['exon_start'], ascending=strand)

            sequence = ''.join(exons_in_transcript.sequence.tolist())

            coding_start = int(exons_in_transcript.loc[exons_in_transcript[coding_pos].ne('').idxmax(), [coding_pos]])

            exon_start = list(filter(None, exons_in_transcript.exon_start.tolist()))
            exon_end = list(filter(None, exons_in_transcript.exon_end.tolist()))

            relative_pos_list = [transcript_id, coding_start, [], []]
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

    def get_rel_pos(atd: pd.DataFrame, pos: int, strand: int) -> int:
        """Finds the relative position inside a sequence based on absolute position in the chromosome."""

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

    def translate_DNA(seq: str) -> str:
        """Translates DNA, not currently in use, but needed for full translation."""

        len_seq = len(seq)
        translated_var_seq = []

        for i in range(0, len_seq - len_seq % 3, 3):
            amino_acid = DNA_Codons[seq[i:i + 3]]
            translated_var_seq.append(amino_acid)
            if amino_acid == '_':
                break

        return ''.join(translated_var_seq)[:-1]

    def SNP_translation(atd: pd.DataFrame, pos: int, var: str, strand: int) -> list:
        """Evaluates if a SNP produces a change in amino acid sequence for a given transcript."""

        coding_start_pos = get_rel_pos(atd, atd.iloc[1], strand)
        SNP_pos = get_rel_pos(atd, pos, strand) - coding_start_pos

        # Prints for checking output
        # print("Exon_end_chosen: " + str(atd.iloc[3][atd.iloc[3] >= atd.iloc[1]].min()))
        # print("Coding_start: " + str(coding_start_pos))

        seq = atd.iloc[4][coding_start_pos-int((1+strand)/2):]
        SNP_codon = seq[SNP_pos - SNP_pos % 3:SNP_pos - SNP_pos % 3 + 3]

        var = var.split('/')
        if strand == -1:
            # The SNPs are always forward strand; needs to be reversed to get accurate results
            var[0] = DNA_ReverseComplement[var[0]]
            var[1] = DNA_ReverseComplement[var[1]]

        SNP_codon_var = SNP_codon[:SNP_pos % 3] + var[1] + SNP_codon[(SNP_pos % 3) + 1:]

        amino_acids = [DNA_Codons[SNP_codon], DNA_Codons[SNP_codon_var]]
        diff_score = compare_and_score_AA(amino_acids)

        if not diff_score:
            return []

        var_amino_acid_pos = (SNP_pos // 3) + 1

        # If full translated seq is needed, uncomment below. Remember to add it to results
        #transltated_var_seq = translate_DNA(seq[:SNP_pos] + var[1] + seq[SNP_pos + 1:])

        return [amino_acids[0] + "/" + amino_acids[1], var_amino_acid_pos, diff_score]

    assembled_transcript_data = pd.DataFrame()
    current_gene = str()
    current_strand = 1
    result_list = []

    for index, data in tqdm(SNP_data.iterrows(), total=SNP_data.shape[0]):

        if data['gene_name'] != current_gene:
            current_gene = data['gene_name']
            gene_df = exons[(exons['gene_name'] == data['gene_name'])]
            current_strand = gene_df.iloc[0, 13]
            assembled_transcript_data = assemble_transcripts(gene_df, current_strand)

        transcript_ids = data['transcript_ids']
        data_as_list = data.tolist()

        for transcript in transcript_ids:
            result = SNP_translation(assembled_transcript_data[assembled_transcript_data['transcript_id'] == transcript]
                                     .iloc[0], data['chrom_pos'], data['variant_alleles'], current_strand)

            if result:
                result_list.append(data_as_list[:6] + [transcript] + [data_as_list[7]] + ['non_synonymous'] + result)
            else:
                result_list.append(data_as_list[:6] + [transcript] + [data_as_list[7]] + ['synonymous'])

    columns = SNP_data.columns.values.tolist() + ['amino_acid_change', 'amino_acid_pos', 'score']
    columns[6] = 'transcript_id'
    columns[8] = 'SNP_type'

    # Returns a list of two dataframes, one with synonymous SNPs and one with non-synonymous SNPs
    return pd.DataFrame(result_list, columns=columns)
