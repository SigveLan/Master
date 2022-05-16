# Human Phenotype Prediction Through use of Variant Data in Metabolic Modeling
## A Master's Project by Sigve Strand Landa

This is a repository for all code related to my master's project


# Pipeline Flow

#### Prepare genome data
- Download data from Ensembl with the biomart tool, using the template as given or use data from the supplementary files:<br /><br />

    - Header setup:<br />
    
    >gene_id|gene_name|chromosome|gene_start(bp)|gene_end(bp)|transcript_id(s)|exon_id|exon_start(bp)|exon_end(bp)|
    phase_start|phase_end|coding_start(bp)|coding_end(bp)|strand
    
    - Header example:<br />
    
    >ENSG00000008130.15|NADK|1|1751232|1780457|ENST00000341426.9;ENST00000341991.7;ENST00000378625.5|
    ENSE00003487616|1761952|1762035|2|2|1761952|1762035|-1
    
    
- Filter the genome data using exon_filter_model using the general metabolic model not tissue specific model. 
  This removes all data for genes not in the model as those are not useful anyways.
  
#### Filter SNPs using individual data:




- This data will then just be a list of SNPs

#### Running the SNP filter

- Tho run the main SNP filter: main_filter.py
    - Input is a list of SNPs
    - For this to work the input will need to be in a certain format.
    - Multiprocessing is available for the filer, whic hcan be useful for large inpput files.
    - There are a total of five output files:
        - noncoding SNPs
        - SNPs in transcript - noncoding
        - SNPs in coding part
            - The SNPs in this part is further scrutinized into:
            - synonymous SNPs: SNPs that do not produce amino acid change
            - missense SNPs: SNPs that do produce amino acid change 


##### VCF extraction
- Using the exon data to filter out SNPs for vcf files: vcf_reader.py
    - Requires Linux

- Combine all SNPs to single file: snps_combine.py
    - Only necessary if multiple files was produced from vcf_reader.py
    - 
- Running the main SNP filter: main.py
- Extracting individual data for all missense SNPs (or other SNPs, depending): fetch_ind_data.py
    - Requires Linux

- Produce combinations: individual_combinations.py

#### PheWAS extraction
- A clean list of SNPs was produced throug extraction from the catalogue. Both files are found in the supplementary data
- Then combinations are generated based on phecodes and all given SNPs associated for the phecode with phewas_extraction.py
    - There is the option of using multiple SNP filter outputs, not just missense SNPs
    - Adding a filter for specific SNPs is possible, but will require some small changes.
    - Combination changes can be changed in srs.phewas_mp_functions.
        - Multiprocessing was added for large combinations
        - If large combinations, and there are many unique SNPs, is to be produced, a phecode filter should be added to phewas_extraction.py
    - The output of phewas_extraction.py can the nbe used for knockout FBA.


#### Combination preparation using general SNP data. 

- The results will then need to be filtered for affected genes using test_files.SNP_result_processing.py
    - This gets all genes that have been affected by at least one SNP
- The affected genes are then used to produce combinations using combination_generator.py
    - There are two options here
        - Extensive combinations, which produces all combinations of certain sizes, ex: all size 3 and 4 combinations
          NB: produces large numbers of combinations

        - Random combinations, which will be a set number of combinations within a range of sizes,
          ex: 200 combinations with sizes ranging from 3 to 6.

- Combinations can then be used to run FBA as above, with results being of the same format.

#### FBA
- scripts/FBA_scripts

- Run FBA with tasks using FBA_with_tasks.py
    - Select which model tissues to use.
    - Select either to filter for (non)essential genes or to only filter for genes in model. 
      Nonessential genes that have been task filtered will need to be
      prepared beforehand. They are also available in the supplementary data for models used in the project. Can also use 'cobra.analysis.find_essential_genes' to
      get genes that produces zero as solution when knocked out. This will however not get genes essential to tasks.
      - The filter here is necessary because the different tissue models have different genes in them. Genes not in the
        model will essentially be ignored for that input.
  

- Interesting results can be filtered out with filter_results.py


- Sub combinations can be created and run with subcombinations_gen_and_FBA.ipynb
  - Selected gene combination will have to pe copy/pasted in

- SNP ids can be retrieved using 'get_SNPs_from_id.py'

- Filtering out samples which contain specific genes in their combination is also possible with 
  get_sample_from_genes.py

  
#### Additional Task Creation
- Use the scripts/full_task_lists_create.ipynb under scripts.


#### FVA / Model explore
- Use the scripts/FBA_scripts/model.explore.ipynb


