def check_if_SNP_in_gene1(SNP_pos):
    global exons
    genes = pd.DataFrame()
    for i, row in exons.iterrows():

        if len(genes) != 0:
            if row['gene_id'] != genes.iloc[0, 1]:
                break

        if (row['gene_start'] <= SNP_pos) & (row['gene_end'] >= SNP_pos):
            genes.append(row)

        elif row['gene_start'] < SNP_pos:
            #continue
            exons.drop(i, inplace=True)

        else:
            return None

    return genes

def check_if_SNP_in_gene2(SNP_pos):
    global exons
    gene_id = ''
    gene = {}
    for i, row in exons.iterrows():

        if (row['gene_start'] <= SNP_pos) & (row['gene_end'] >= SNP_pos):
            if gene_id == '':
                gene_id = row['gene_id']
            if row['gene_id'] == gene_id:
                gene[i] = row

        elif row['gene_start'] < SNP_pos:
            # Dropping exons as the program runs for increased speed.
            exons.drop(i, inplace=True)

        else:

            if not bool(gene):
                return None
            else:
                return create_df(gene)