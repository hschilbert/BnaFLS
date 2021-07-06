# BnaFLSs 
scripts belonging to "Characterization of the Brassica napus flavonol synthase gene family reveals bifunctional flavonol synthases" https://doi.org/10.1101/2021.06.30.450533 

### Please get in touch if you need help running the scripts on your own data sets: [Hanna Schilbert (email)](mailto:hschilbe@cebitec.uni-bielefeld.de?subject=[GitHub]BnaFLSs_scripts_request) ###

## Usage

### General recommendation

Full paths should be used to specify input and output files and folders. Sequence names should not contain white space characters like spaces and TABs. Underscores can be used to replace spaces.

## coexp_bn.py
This script identifies co-expressed genes. 

```
Usage:
  python coexp_bn.py --in <FILE> --exp <FILE> --out <DIR>
  
  Mandatory:
  
  Inputs 
  --in     STR     gene IDs, one gene ID per line
  --exp    STR     RNA-Seq count table, columns containing RNA-Seq sample IDs, rows gene IDs
  
  Output directory
  --out    STR     Output directory
  
  Optional:
  --ann    STR     functional annotation file, first column gene IDs, second column functional annotation
```

`--in` txt-file that contains a gene ID or a set of gene IDs for which co-expressed genes should be identified, where one gene ID is listed in one row. The gene ID must match with the gene IDs in the files given at `--exp` and `--ann`.

`--exp` count table/gene expression file containing e.g. the TPM values where the columns contain the RNA-Seq sample IDs and the rows contain the gene IDs. Thus more than one RNA-Seq sample can be analysed at the same time. 

`--out` specify the output directory where the results should be stored

`--ann` functional annotation where one the first column contains gene IDs and the second column contains the functional annotation. Contains a header "Gene_ID\tAnnotation"

All files should be provided in tab separated format.

## heatmap.py
This script generates a heatmap for a gene of interest (or a set of genes of interest) on the basis of a RNA-Seq count table and a functional annotation. 

```
Usage:
  python heatmap.py --exp <FILE> --genes <FILE> --samples <FILE> --out <FILE>
  
  Mandatory:
  
  Inputs 
  --exp       STR     RNA-Seq count table, columns containing RNA-Seq sample IDs, rows gene IDs
  --genes     STR     first column gene IDs, second column corresponding functional annotation
  --samples   STR     first column tissue name, second column corresponding RNA-Seq sample ID(s)
  
  Output directory
  --out       STR     Output file .pdf
  
  Optional:
  --zscore    if provided activates zscore normalization
```

`--exp` count table/gene expression file containing e.g. the TPM values where the columns contain the RNA-Seq sample IDs and the rows contain the gene IDs. Thus more than one RNA-Seq sample can be analysed at the same time. 

`--genes` txt-file that contains a gene ID or a set of gene IDs for which the heatmap should be generated. Gene ID(s) is listed in the first column and the corresponding functional annotation in the second row. The gene ID must match with the gene IDs in the files given at `--exp`.

`--samples` txt-file that contains in the first column the tissue name(s) and in the second column corresponding RNA-Seq sample ID(s). If one tissue comprises multiple RNA-Seq data sets, the RNA-Seq sample IDs should be provided comma separated. The RNA-Seq sample ID(s) must match with the IDs in `--exp`.

`--out` specify the output file (.pdf) where the results should be stored. In addition the raw values used for the generation of the heatmap will be given in a txt-file.

`--zscore` if given zscore normalization is activated

All files should be provided in tab separated format if not stated otherwise.

## pairwise_comp.py
This script calculates pairwise amino acid sequence identities [%].

```
Usage:
  python pairwise_comp.py --in <FILE> --out <FILE>
  
  Mandatory:
  
  Inputs 
  --in     STR     FASTA file, containing amino acid sequences
  
  Output directory
  --out    STR     Output file, where amino acid sequence identities [%] are listed
  
```
