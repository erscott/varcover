# VarCover 

VarCover is an open-source software package and web-interface that provides scientists and clinical laboratories with a simple method to assemble a validation panel with an approximately minimal set of samples needed to cover all possible target variants.

VarCover employs the:  
* SetCoverPy package
* sample weights
* pre-selection of singleton-possessing samples 

to efficiently solve the min-set cover problem. 

#### Sample weights derived from the minor allele frequency spectrum increased the number of alleles in the solution set. 

#### Pre-selection of samples possessing singleton target alleles reduced computational processing time if the target set size exceeded 100 alleles.

<br>

## Command Line Interface
VarCover CLI accepts a VCF file containing the variants of interest and potential samples.  

| Argument               | Description      |   
| -----------------------|:----------------:| 
| vcf                    | path to vcf      | 
| -w, --weight (standard/logit)           | use sample allele frequency weighting      |   
| -s, --singletonReduction | use singleton reduction modification      |    
| -n, --niters (integer) | number of PySetCover iterations      |    
| -o, --output_dir (/path/to/dir) | specifies output directory      |    

Optional arguments:
1) -w, --weight: whether sample weights derived from the minor allele frequency spectrum

2) -s, --singletonReduction: whether pre-selection of singleton samples should be used.

example command:
python run_varcover.py multisample_multivariant.vcf.gz -s -w logit -o /work/varcover/

