# ts65Dn-gene-integration

### Detailed walkthrough of GTEx integration analysis for ts65Dn DEGs to infer endocrine signaling.  Full integration analysis available in R script using datasets below
### All preprocessed files used for analysis are available here: https://drive.google.com/drive/folders/1xeybWL5frJ0SWHBmz8xED-wTX8pk266G?usp=sharing

### Datasets used in this study consist of:
#### 1. R environment containing GTEx data filtered for gene expression across all tissues: GTEx NA included env.RData
#### source: https://doi.org/10.7554/eLife.76887 and https://doi.org/10.7554/eLife.76887

#### 2. List of human genes annotated to encode secreted proteins: human secreted proteins.tab
#### source: UniProt

#### 3. List of mouse and human orthologue genes: Mouse Gene info with Human Orthologues
#### source: MGI

#### 4. Go term annotations for human genes: uniprot-human-genes and goterms mapping.tab
#### source: UniProt

#### 5. List of Differential expression results from pan-tissue RNA-seq of Ts65Dn mice compared to WT: mouse DEGS/
#### source: this study
