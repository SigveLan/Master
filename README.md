# Human Phenotype Prediction Through use of Variant Data in Metabolic Modeling
## A Master's Project by Sigve Strand Landa

This is a repository for all code used in my master's project.
Included is a list of data sources as well as simple instructions on how to use the code.

## Data Sources

- Genome Data:
    https://www.ensembl.org/biomart/martview/5d15d7e8e172f9cf98740393c44be6a1

- Models:
    https://github.com/SysBioChalmers/Human-GEM

- 1000 Genomes Project:
    https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/

- PheWAS Data:
    https://phewascatalog.org/phewas
    

## Pipeline

#### Prepare genome data

- Download data from Ensembl with the biomart tool, using the template as given or use data from the supplementary files:<br /><br />

    - Header setup:<br />
    
    >gene_id|gene_name|chromosome|gene_start(bp)|gene_end(bp)|transcript_id(s)|exon_id|exon_start(bp)|exon_end(bp)|
    phase_start|phase_end|coding_start(bp)|coding_end(bp)|strand
    
    - Header example:<br />
    
    >ENSG00000008130.15|NADK|1|1751232|1780457|ENST00000341426.9;ENST00000341991.7;ENST00000378625.5|
    ENSE00003487616|1761952|1762035|2|2|1761952|1762035|-1
    
    
- Filter the genome data using scripts/exon_filter_model.py using the general metabolic model genes. 
  This removes all data for genes not in the model as those are not useful here anyways.
  

#### Running the SNP filter

- To run the main SNP filter: main_filter.py
- Input is a list of SNPs with a certain format
  
- Multiprocessing is available for the filer, which can be useful for larger input files.

- There are a total of five output files:
    - noncoding SNPs
    - SNPs in transcript - noncoding
    - SNPs in coding part
        - The SNPs in this part is further scrutinized into:
        - synonymous SNPs: SNPs that do not produce amino acid change
        - missense SNPs: SNPs that do produce amino acid change 


#### VCF extraction

- Files in: scripts/vcf_scripts/

- Using the exon data to filter out SNPs for vcf files: vcf_reader.py
    - Requires Linux

- Combine all SNP files to single file: snps_combine.py

- Running the main SNP filter as above.
- Extracting individual data for all missense SNPs (or other SNPs, depending): fetch_ind_data.py
    - Requires Linux

- Produce combinations: individual_combinations.py


#### PheWAS extraction

- A clean list of SNPs was produced through extraction from the catalogue. Both files are found in the supplementary data.

- Then combinations are generated based on phecodes and all given SNPs associated for the phecode with phewas_extraction.py
    - There is the option of using multiple SNP filter outputs, not just missense SNPs
    - Adding a filter for specific SNPs is possible, but will require some small changes.
  
    - Combination changes can be changed in srs.phewas_mp_functions.
        - Multiprocessing was added for large combinations.
        - If large combinations, and there are many unique SNPs, is to be produced, a phecode filter should be added to phewas_extraction.py
      
    - The output of phewas_extraction.py can then be used for knockout FBA.



#### Combination preparation using general SNP data. 

- Run the SNP filter with the list of SNPs

- The results will then need to be filtered for affected genes using test_files/SNP_result_processing.py
    - This gets all genes that have been affected by at least one SNP.
- The affected genes are then used to produce combinations using scripts/combination_generator.py
    - There are two options here:
        - Extensive combinations, which produces all combinations of certain sizes, ex: all size 3 and 4 combinations
          NB: produces large numbers of combinations

        - Random combinations, which will be a set number of combinations within a range of sizes,
          ex: 200 combinations with sizes ranging from 3 to 6.

- Combinations can then be used to run FBA, with results being of the same format.


#### FBA

- Files in scripts/FBA_scripts/

- Load any file with gene combinations

- If PheWAS data is used, select the phecode to run.

- Run FBA with tasks using FBA_with_tasks.py
    - Select which model tissues to use.
    - Select either to filter for (non)essential genes or to only filter for genes in model. 
      Nonessential genes that have been task filtered will need to be
      prepared beforehand. They are also available in the supplementary data for models used in the project. Can also use 'cobra.analysis.find_essential_genes' to
      get genes that produces zero as solution when knocked out. This will however not get genes essential to tasks.
      - The filter here is necessary because the different tissue models have different genes in them. Genes not in the
        model will essentially be ignored for that input.
        
- If tasks are not used, FBA_simple.py can be utilized for FBA instead. 
  

#### Result Processing

- Result processing varies slightly depending on inputs used.

- Interesting results for any output file can be filtered out with filter_results.py

- Sub combinations can be created and run with: scripts/FBA_scripts/subcombinations_gen_and_FBA.ipynb
  - Selected gene combination will have to pe copy/pasted in

- SNP ids can be retrieved using: scripts/get_SNPs_from_id.py

- Filtering out samples which contain specific genes in their combination is also possible with:
  scripts/get_sample_from_genes.py
  
- Visually clean task results outputs are generated with: scripts/results_clean_task_outputs.py


#### Task Functions
- All functions used to prepare the internal task list from file, and models for task usage are located in src.task_functions.py

- When FBA is done with tasks a special FBA function is used. This functions applies each task to the appropriate model then performs FBA. This function is located in src.mp_function, function: knockout_FBA_w_tasks
  

#### Additional Task Lists Creation
- Use the scripts/full_task_lists_create.ipynb recreate the process of generating additional tasks. Already created lists are also in the supplementary data.

- One task, 'Storage of Glucose in Glycogen' should be removed as it has properties not supported.

- The original full task list is read in and excess columns are removed.

- Essential tasks are filtered away as they are already listed in essential tasks.

- Instructions are listed in the file for removal of tasks with missing metabolites.

- Any task that fails for a given tissue is removed.

- All tasks that passed for a tissue is then written to file.


#### FVA and Model Exploration
- Use the scripts/FBA_scripts/model.explore.ipynb

- In this case the FVA test uses the generic Human1 model.

- A long list of constraints will be added to the model.

- Then FVA is run with and without HIBCH knockout.

- To see what the reactions are, visit: https://metabolicatlas.org/explore/Human-GEM/gem-browser
  - The GEM browser can also be used for exploration 


#### Essential Genes analyses

- Number of essential genes in a tissue compared to total number of genes: 
  - scripts/essential_percentages.py

- Number of affected genes in a combination that are essential, nonessential, or not-in-tissue:
  - scripts/SNP_essential_percent.py


