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

# From translation - old
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


'''
    print(exon_end)
    print(pos)
    print(rel_pos)
    print(coding_start_pos)
    print(SNP_codon)

    print(seq)
    print(atd.iloc[3])
'''
'''
if not list(filter(None, exons_in_transcript.coding_start.tolist())):
    # Skips non coding transcripts
    continue

'''
"""
exon_coding_start = atd.iloc[3][atd.iloc[3] >= atd.iloc[1]].min()

if strand == 1:
    coding_start_pos = int(
        atd.iloc[2][np.where(atd.iloc[3] == exon_coding_start)[0][0]] - abs(exon_coding_start - atd.iloc[1]))
else:
    indices = np.where(atd.iloc[3] > exon_coding_start)[0]
    if indices.size != 0:
        coding_start_pos = indices[-1] + abs(exon_coding_start - atd.iloc[1])
    else:
        coding_start_pos = abs(exon_coding_start - atd.iloc[1])

exon_end = atd.iloc[3][atd.iloc[3] >= pos].min()

    if strand == 1:
        SNP_pos = int(atd.iloc[2][np.where(atd.iloc[3] == exon_end)[0][0]] - abs(exon_end - pos))
    else:
        indices = np.where(atd.iloc[3] > exon_end)[0]
        if indices.size != 0:
            SNP_pos = indices[-1] + abs(exon_end - pos)
        else:
            SNP_pos = abs(exon_end - pos)
"""