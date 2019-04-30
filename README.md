# VarCover 

VarCover is an open-source software package and web-interface that provides scientists and clinical laboratories with a simple method to assemble a validation panel with an approximately minimal set of samples needed to cover all possible target variants.

VarCover employs the:  
* SetCoverPy package
* sample weights
* pre-selection of singleton-possessing samples 

to efficiently solve the min-set cover problem. 

## Dependencies  
* python ≥3.6  
* numpy  
* pandas  
* SetCoverPy  
* pandasVCF

<br>

## Command Line Interface
VarCover CLI accepts a VCF file containing the variants of interest and available samples.  

| Argument               | Description      | Default      |   
| -----------------------|:----------------:| :----------------:|
| vcf                    | path to vcf      |    -     |
| -w, --weight (standard/logit)           | use sample allele frequency weighting |  standard  |   
| -s, --singletonReduction | use singleton reduction modification      |  singletonReduction off  |
| -n, --niters (integer) | number of PySetCover iterations      |  20  |
| -o, --output_dir (/path/to/dir) | specifies output directory      |  current working directory  |    

### Optional arguments:
1) -w, --weight: whether sample weights derived from the minor allele frequency spectrum
    * Sample weights derived from the minor allele frequency spectrum can increase the number of variants in the solution set, though it may also break the minimum set cover requirement.  Experiments with ~250 variants led to an increase in the solution sample size by 1.

2) -s, --singletonReduction: whether pre-selection of singleton samples should be used.
    * Pre-selection of samples possessing singleton target alleles may reduce computational processing time if the target set size exceeds 100 alleles.

### Output:
Three files will be created based upon the VarCover solution set:
1) varcover_solution_allele_cover_{randomid}.tsv - CHROM POS REF ALT ALLELE_COUNT
2) varcover_solution_df_{randomid}.tsv - CHROM POS REF ALT SAMPLE1 SAMPLE2 ...
3) varcover_solution_samples_{randomid}.tsv - Sampleids (1 per row)


### Example command:
python run_varcover.py 1000gPh3_hg19_acmg57_vars_allchroms.vcf.gz -s -w logit -o /work/varcover/

* 1000gPh3_hg19_acmg57_vars_allchroms.vcf.gz can be found in the /varcover/varcover/data folder
 





