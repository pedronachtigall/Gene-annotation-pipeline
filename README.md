# Gene-annotation-pipeline
Bioinformatics pipeline and tutorial for performing gene annotation in genome assemblies using GALBA

:construction:	**Under construction!** :construction:	



## Dependencies
 - [Python](https://www.python.org/) and [biopython](https://biopython.org/)
 - [GALBA](https://github.com/Gaius-Augustus/GALBA)
 - [CodAn](https://github.com/pedronachtigall/CodAn)
 - [StringTie](https://github.com/gpertea/stringtie)
 - [STAR](https://github.com/alexdobin/STAR/releases)
 - [Orthofinder](https://github.com/davidemms/OrthoFinder)

Ensure that all dependencies are installed and working properly.

<!---
## Model species
We will use the Golden lancehead (*Bothrops insularis*) as a model for this tutorial.

This genome is available in the NCBI ([PRJNA679826](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679826/)) and [figshare](https://figshare.com/projects/Bothrops_insularis_genome/237995) repositories.

We used the soft-maked primary assembly to perform the gene annotation, which can be obtained following the "[Repeat Annotation](https://github.com/pedronachtigall/Repeat-annotation-pipeline)" tutorial.

We will use the following RNA-seq data as transcript evidence:
The raw data is listed below:
| Sample ID | Data type | NCBI accession |
| :-------- | :-------: | :------------: | 
| SB1851_VG_rna  | RNA-seq | SRR32358140 |
| SB1851_HG_rna  | RNA-seq | SRR32358139 |
| SB1851_ILG_rna  | RNA-seq | SRR32358138 |
| SB1851_Pancreas_rna  | RNA-seq | SRR32358137 |
| SB1851_Muscle_rna  | RNA-seq | SRR32358136 |
| SB1851_Heart_rna  | RNA-seq | SRR32358151 |
| SB1851_Brain_rna  | RNA-seq | SRR32358150 |
| SB1851_Spleen_rna  | RNA-seq | SRR32358149 |
| SB1851_Kidney_rna  | RNA-seq | SRR32358148 |
| SB1851_Ovary_rna  | RNA-seq | SRR32358147 |

## Retrieve transcripts and proteins from RNA-seq data
```
#trim reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o ${SAMPLE}_tg ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz


#index gthe genome
${SAMPLE}_R1_val_1.fastq.gz ${SAMPLE}_R2_val_2.fastq.gz

#map reads

#retrieve transcripts

```

### Predict proteins using CodAn
```
wget https://github.com/pedronachtigall/CodAn/blob/master/models/VERT_full.zip
unzip VERT_full.zip
codan.py -t transcripts.fa -o output_folder -m VERT_full/
```

 - It is the transcript-derived proteins for the target species, which will be used as transcript evidence.

## Retrieve proteins from other snake species
Survey protein databases to ensure a compelte protein set to be used as protein evidence.

Here we are surveying the ENSEMBL, Uniprot, and NCBI databases.
```
ADD CODE
```

## Merge the transcript-dervied proteins to the other snake species proteins
```
ADD CODE
```

## Perform gene annotation using GALBA
```
ADD CODE
```
