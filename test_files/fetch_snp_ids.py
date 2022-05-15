import pandas as pd
from tqdm import tqdm
from Bio import Entrez

snp_data = pd.read_table('C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/result_chrom_22.tsv', index_col=0)

first_value = snp_data['pos'].iloc[0]
last_value = snp_data['pos'].iloc[0]

result = []

for index, data in tqdm(snp_data.iterrows(), total=snp_data.shape[0]):

    if data['pos'] - last_value < 1000:
        last_value = data['pos']
    else:
        result.append(''.join(['22[Chromosome] AND ', str(first_value-2)+'[CHRPOS] : ', str(last_value+2)+'[CHRPOS]']))
        first_value = data['pos']
        last_value = data['pos']
    if len(result) == 2:
        break

print(*result, sep=',')

#with open('C:/Users/Sigve/Genome_Data/SNP_data/1000_genomes/result_chrom_22_regions.txt', 'w') as text_file:
    #text_file.write(','.join(result))

Entrez.email = ''
handle = Entrez.esearch(db='snp', term=result[0])
records = Entrez.read(handle)
for record in records:
    print(record.IdList)


handle.close()

