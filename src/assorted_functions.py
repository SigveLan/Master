import pandas as pd
from Bio import SeqIO
from operator import itemgetter


def SNP_sort(SNP_file: str, write_results_to_file=True) -> pd.DataFrame:
    """A function for sorting SNPs by chromosome and position. Also removes duplicates and non-SNP mutations,
    as well as SNPs with more than one variant base."""

    SNP_df = pd.read_table(SNP_file)

    SNP_df.drop_duplicates(subset=['Variant name'], inplace=True)
    SNP_df.sort_values(by=['Chromosome/scaffold name', 'Chromosome/scaffold position start (bp)'], inplace=True)
    SNP_df.set_index(['Variant name'], drop=True, inplace=True)
    SNP_df['Variant alleles'] = SNP_df['Variant alleles'].apply(str)
    SNP_df = SNP_df[SNP_df['Variant alleles'].str.contains("^[ACTG]/[ACTG]$")]

    SNP_df['Chromosome/scaffold position start (bp)'] = SNP_df['Chromosome/scaffold position start (bp)'].apply(int)
    SNP_df['Chromosome/scaffold position end (bp)'] = SNP_df['Chromosome/scaffold position end (bp)'].apply(int)
    SNP_df['Strand'] = SNP_df['Strand'].apply(int)

    if write_results_to_file:
        SNP_df.to_csv('C:/Users/Sigve/Genome_Data/SNP_data/biomart_chrom_all_test.tsv', sep='\t')

    return SNP_df


def read_exons_to_df(exon_file: str) -> pd.DataFrame:
    """Reads in genome data in the form of an exon fasta file.

    Header setup:
    >gene_id|gene_name|chromosome|gene_start(bp)|gene_end(bp)|transcript_id(s)|exon_id|exon_start(bp)|exon_end(bp)|
    phase_start|phase_end|coding_start(bp)|coding_end(bp)|strand

    Header example:
    >ENSG00000008130.15|NADK|1|1751232|1780457|ENST00000341426.9;ENST00000341991.7;ENST00000378625.5|
    ENSE00003487616|1761952|1762035|2|2|1761952|1762035|-1
    """

    # Reading to dictionary then convert to dataframe is very fast
    dictionary = {}
    ind = 0

    for seq_record in SeqIO.parse(open(exon_file, mode='r'), 'fasta'):

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

        # There can be more than one coding start/stop in a given exon. This picks the closest to the edge of the exon.
        for i in [11, 12]:
            if seq_record_attributes[i] != '':
                pos = list(map(int, seq_record_attributes[i].split(';')))
                if len(pos) > 1:
                    differences = [abs(j - seq_record_attributes[i-4]) for j in pos]
                    seq_record_attributes[i] = pos[min(enumerate(differences), key=itemgetter(1))[0]]

                else:
                    seq_record_attributes[i] = int(seq_record_attributes[i])

        # Sequence is included for the translation/SNP eval script
        seq_record_attributes.append(str(seq_record.seq))

        dictionary[ind] = seq_record_attributes
        ind += 1

    return pd.DataFrame.from_dict(dictionary, orient='index', columns=["gene_id", "gene_name", "chrom", "gene_start",
                                                                     "gene_end", "transcript_ids", "exon_id",
                                                                     "exon_start", "exon_end", "phase_start",
                                                                     "phase_end", "coding_start", "coding_end",
                                                                     "strand", "sequence"])


def add_cbm_id(model_file: str, SNP_data: pd.DataFrame, write_results_to_file=True) -> pd.DataFrame:
    """A function to add gene number from the CBM model to the SNP results"""

    model_df = pd.read_table(model_file)
    model_df = model_df[['gene_number', 'ensembl_gene']]

    model_df['gene_number'] = model_df['gene_number'].astype(int)
    model_df.drop_duplicates(subset=['gene_number'], inplace=True)

    SNP_data.set_index(SNP_data['gene_id'].apply(lambda x: str(x).split('.')[0]), inplace=True)

    joined_df = SNP_data.join(model_df.set_index('ensembl_gene'))
    joined_df.dropna(subset=['gene_number'], inplace=True)

    joined_df.reset_index(drop=True, inplace=True)
    joined_df['gene_number'] = joined_df['gene_number'].astype(int)

    if write_results_to_file:
        joined_df.to_csv(path_or_buf='C:/Users/Sigve/Genome_Data/results/SNPs_effect_mod.tsv', sep='\t')

    return joined_df
