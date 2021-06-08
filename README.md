# BnaFLSs 
scripts belonging to BnaFLSs publication

### Please get in touch if you need help running the scripts on your own data sets: [Hanna Schilbert (email)](mailto:hschilbe@cebitec.uni-bielefeld.de?subject=[GitHub]BnaFLSs_scripts_request) ###

## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces.

## coexp_bn.py
This script calculated co-expressed genes for a gene of interest on the basis of a RNA-Seq count table and a functional annotation. 

```
Usage:
  python coexp_bn.py --in <FILE> --exp <FILE> --out <DIR>
  
  Mandatory:
  
  Input gene IDs
  --in  STR gene IDs, one gene ID per line
  
  
  --exp
  
  Output directory
  --out            STR    Output directory
  
  Optional:
  --ann <FILE> 
```

--in <txt-file that contains a set of gene IDs for which co-expressed genes should be identified, where one gene ID is listed in one row.
--exp <count table containing e.g. the TPM values>
--ann <functional annotation where one coloum contains
--out <specify the output directory where the results should be stored>
