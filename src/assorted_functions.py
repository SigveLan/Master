import pandas as pd
from Bio import SeqIO
from operator import itemgetter


def SNP_sort(SNP_file: str) -> pd.DataFrame:
    """A function for sorting SNPs by chromosome and position. Also removes duplicates and non SNP mutations,
    as well as SNPs with more than one variant base."""

    SNP_df = pd.read_table(SNP_file)
    SNP_df = SNP_df[['#chr', 'pos', 'variation', 'snp_id']]
    SNP_df = SNP_df[SNP_df['#chr'] != '#chr']
    SNP_df.reset_index(inplace=True, drop=True)

    SNP_df['variation'] = SNP_df['variation'].apply(lambda x: x.split('>'))
    SNP_df[['ancestral', 'variation']] = pd.DataFrame(SNP_df['variation'].tolist())
    SNP_df['variation'] = SNP_df['variation'].apply(lambda x: x.split(','))
    SNP_df = SNP_df.explode('variation', ignore_index=True)
    SNP_df['snp_id'] = SNP_df['snp_id'].apply(lambda x: 'rs' + str(x))
    SNP_df['variation'] = SNP_df['ancestral'] + '/' + SNP_df['variation']
    SNP_df.drop(columns='ancestral', inplace=True)
    SNP_df = SNP_df[['snp_id', '#chr', 'pos', 'variation']]
    SNP_df.rename(columns={'#chr': 'chr', 'variation': 'var'}, inplace=True)

    SNP_df = SNP_df[SNP_df['var'].str.contains("^[ACTG]/[ACTG]$")]
    SNP_df.sort_values(by=['chr'], inplace=True)
    SNP_df.reset_index(inplace=True, drop=True)

    return SNP_df


def read_exons_to_df(exon_file_path: str) -> pd.DataFrame:
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

    for seq_record in SeqIO.parse(open(exon_file_path, mode='r'), 'fasta'):

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


def add_cbm_id(model_data_path: str, SNP_data: pd.DataFrame) -> pd.DataFrame:
    """A function to add gene number from the CBM (Recon 3D) to the SNP results."""
    model_df = pd.read_table(model_data_path)
    model_df = model_df[['gene_number', 'symbol']]

    model_df['gene_number'] = model_df['gene_number'].astype(int)
    model_df.drop_duplicates(subset=['gene_number'], inplace=True)

    SNP_data.set_index(SNP_data['gene_name'].apply(lambda x: str(x).split('.')[0]), inplace=True)
    joined_df = SNP_data.join(model_df.set_index('symbol'))

    # It filters out genes without a number.
    df_to_tsv(joined_df, 'C:/Users/Sigve/Genome_Data/results/test.tsv')
    joined_df.dropna(subset=['gene_number'], inplace=True)
    joined_df[['amino_acid_pos', 'protein_length', 'gene_number']] = joined_df[['amino_acid_pos', 'protein_length', 'gene_number']].astype(int)

    joined_df.rename(columns={'gene_number': 'model_gene_number'}, inplace=True)
    joined_df.sort_values(by=['chrom', 'chrom_pos'], inplace=True)
    joined_df.reset_index(drop=True, inplace=True)

    return joined_df


def df_to_tsv(df: pd.DataFrame, file_path: str) -> None:
    """Simply writes a dataframe to file as a tsv."""
    df.to_csv(path_or_buf=file_path, sep='\t')


def split_df_results_by_column(filter_results: pd.DataFrame, column: str, values: list) -> list:
    """Takes in the results from SNP_filter and splits it into three dataframes"""
    return [filter_results[filter_results[column] == value].reset_index(drop=True) for value in values]
