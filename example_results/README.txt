Result explanation

The first part of the process identifies SNPs inside gene regions as they are defined by Ensembl. The SNPs are then divided into two categories, those in exons and those not in exons.
These results are then written to the two files: 'SNPs_coding' and ' SNPs_non_coding'. (might change this)Note that SNPs inside exons may still not be part of the coding region, but will still be in the coding file.

The results from the 'SNPs_coding' file is then used to check if the SNPs lead to a change in amino acid sequence. This produces the file 'SNPs_effect'.