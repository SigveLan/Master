# Master
A repository for all code related to my masters project


# Pipeline Flow

####Prepare genome data
- Download data from Ensembl with the biomart tool, using the template as given:<br /><br />

    - Header setup:<br />
    
    >gene_id|gene_name|chromosome|gene_start(bp)|gene_end(bp)|transcript_id(s)|exon_id|exon_start(bp)|exon_end(bp)|
    phase_start|phase_end|coding_start(bp)|coding_end(bp)|strand
    
    - Header example:<br />
    
    >ENSG00000008130.15|NADK|1|1751232|1780457|ENST00000341426.9;ENST00000341991.7;ENST00000378625.5|
    ENSE00003487616|1761952|1762035|2|2|1761952|1762035|-1
    
    
- Filter the genome data using exon_filter_model using the general metabolic model not tissue specific model. 
  This removes all data for genes not in the model as those are not useful anyways.
  
####Filter SNPs using individual data:


- Using the exon data to filter out SNPs for vcf files: vcf_reader.py
    - Requires Linux

- Combine all SNPs to single file: snps_combine.py
- Running the main SNP filter: main.py
- Extracting individual data for all missense SNPs (or other SNPs, depending): fetch_ind_data.py
    - Requires Linux

- Produce combinations: individual_combinations.py

####Filter SNPs using general SNP data. 

- This data will then just be a list of SNPs

- Then run the main SNP filter: main.py
    - You might need to clean up the SNP data a bit for this to work,
      making sure the columns have correct labels and such

- The results will then need to be filtered for affected genes using test_files.SNP_result_processing.py
    - This gets all genes that have been affected by at least one SNP
- The affected genes are then used to produce combinations using combination_generator.py
    - There are two options here
        - Extensive combinations, which produces all combinations of certain sizes, ex: all size 3 and 4 combinations
          NB: produces large numbers of combinations

        - Random combinations, which will be a set number of combinations within a range of sizes,
          ex: 200 combinations with sizes ranging from 3 to 6.

- Combinations can then be used to run FBA as above, with results being of the same format.

####FBA
- Run FBA with tasks using FBA_with_tasks.py
    - Select which model tissues to use.
    - Select either to filter for (non)essential genes or to only filter for genes in model. 
      Nonessential genes that have been task filtered will need to be
      prepared beforehand. They are also evailable in the supplementary data for models used in the project. Can also use 'cobra.analysis.find_essential_genes' to
      get genes that produces zero as solution when knocked out. This will however not get genes essential to tasks.
      - The filter here is necessary because the different tissue models have different genes in them. Genes not in the
        model will essentiallt be ignored for that input.
  

- Interesting results can be filtered out with filter_results.py
  and SNP ids can be retrieved using 'get_SNPs_from_id.py'
- Filtering out samples with containing specific genes in their combination is also possible with 
  get_sample_from_genes.py



####Additional Task Creation




