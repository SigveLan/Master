import pandas as pd

"""Produces individual combinations based on individual data. Possible to add SNP filter here as well, so only 
certain SNPs are included. Can also choose between homo- and hetero-zygote. Also combines identical combinations to 
single entry. For larger number of SNPs, combinations are virtually guaranteed to be unique."""


def main():

    ind_data = pd.read_table('C:/Users/Sigve/Genome_Data/results/ind_combinations/all_missense_ind_data.tsv', index_col=0)
    output_file = 'C:/Users/Sigve/Genome_Data/results/ind_combinations/start_stop_comb_het.tsv'

    # Can filter here (currently filtered fir SNPs affecting start/stop codons).
    if True:
        ind_data = ind_data[ind_data['amino_acid_change'].str.contains('_') |
                            (ind_data['amino_acid_pos'] == 1)].reset_index(drop=True)

    res_dict = {}
    n = 0
    for column in ind_data:
        n += 1
        if n < 13:

            continue

        # For hetero-zygote: set '> 1' to '> 0'.
        res_dict[column] = [';'.join(ind_data['variant_name'][ind_data[column] > 0].to_list()), ';'.join(set([i.split('.')[0] for i in ind_data['gene_id'][ind_data[column] > 0].to_list()]))]

    results = pd.DataFrame.from_dict(res_dict, orient='index', columns=['snp_list', 'gene_ids'])
    t = results.reset_index().groupby(by='gene_ids')[['index', 'gene_ids']].agg(set).reset_index(drop=True)
    t.rename(columns={'index': 'sample_ids'}, inplace=True)
    t.to_csv(path_or_buf=output_file, sep='\t')


if __name__ == '__main__':
    main()
