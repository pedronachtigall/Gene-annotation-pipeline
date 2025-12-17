# Gene-annotation-pipeline
Bioinformatics pipeline and tutorial for performing gene annotation in genome assemblies using GALBA

## Dependencies
 - [Python](https://www.python.org/) and [biopython](https://biopython.org/)
 - [GALBA](https://github.com/Gaius-Augustus/GALBA)
 - [CodAn](https://github.com/pedronachtigall/CodAn)
 - [StringTie](https://github.com/gpertea/stringtie)
 - [HISAT2](https://github.com/DaehwanKimLab/hisat2)
 - [GffRead](https://github.com/gpertea/gffread)
 - [Orthofinder](https://github.com/davidemms/OrthoFinder)

Ensure that all dependencies are installed and working properly.

## Model species
We will use the Golden lancehead (*Bothrops insularis*) as a model for this tutorial ([Nachtigall et al., 2025](https://doi.org/10.1093/gbe/evaf243)).

This genome is available in the NCBI ([PRJNA679826](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA679826/)) and [figshare](https://figshare.com/projects/Bothrops_insularis_genome/237995) repositories.

It is recommended that you run GALBA on genomic sequences that have been softmasked for Repeats. You can follow the "[Repeat Annotation](https://github.com/pedronachtigall/Repeat-annotation-pipeline)" tutorial to soft-mask the genome.

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

#index the genome
hisat2-build -p 20 Binsularis_primary_chromosomes.fasta GINDEX

#map reads
hisat2 -p 20 --rg-id ${SAMPLE} --rg SM:${SAMPLE} --summary-file ${SAMPLE}_hisat2_summary.txt -x GINDEX -1 ${SAMPLE}_R1_val_1.fastq.gz -2 ${SAMPLE}_R2_val_2.fastq.gz -S ${SAMPLE}.sam
samtools view -@ 20 -b -S -o ${SAMPLE}.bam ${SAMPLE}.sam
rm ${SAMPLE}.sam
samtools view -@ 20 -b -F 4 ${SAMPLE}.bam > ${SAMPLE}_mapped.bam
rm ${SAMPLE}.bam
samtools sort -@ 20 ${SAMPLE}_mapped.bam -o ${SAMPLE}_mapped.sorted.bam"
rm ${SAMPLE}_mapped.bam
samtools index ${SAMPLE}_mapped.sorted.bam

#merge bam files
samtools merge -o merged.bam *.mapped.sorted.bam
samtools index merged.bam

#retrieve transcripts
stringtie -p 20 -o transcripts.gtf merged.bam
gffread -w transcripts.fasta -g Binsularis_primary_chromosomes.fasta transcripts.gtf
```

### Predict proteins using CodAn
```
wget https://github.com/pedronachtigall/CodAn/blob/master/models/VERT_full.zip
unzip VERT_full.zip
codan.py -t transcripts.fasta -o CodAn_output/ -m VERT_full/
```
 - The file ```CodAn_output/PEP_sequences.fasta``` is the transcript-derived proteins for the target species, which will be used as transcript evidence.

## Retrieve proteins from other species
Survey protein databases to ensure a compelte protein set to be used as protein evidence.

Here, we are surveying the [ENSEMBL](https://www.ensembl.org/index.html), [UniProt](https://www.uniprot.org/), and [NCBI](https://www.ncbi.nlm.nih.gov/) databases. Specifically, we will retrieve sequences from snakes, lizard, chicken, and mouse. But you can retrieve seqeunces from more species (and also modifiy based on your target lineage).
```
#snakes
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/797/435/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/769/535/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_protein.faa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/pseudonaja_textilis/pep/Pseudonaja_textilis.EBS10Xv2-PRI.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/notechis_scutatus/pep/Notechis_scutatus.TS10Xv2-PRI.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/naja_naja/pep/Naja_naja.Nana_v5.pep.all.fa.gz

#lizards
wget https://ftp.ensembl.org/pub/release-115/fasta/anolis_carolinensis/pep/Anolis_carolinensis.AnoCar2.0v2.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/pogona_vitticeps/pep/Pogona_vitticeps.pvi1.1.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/varanus_komodoensis/pep/Varanus_komodoensis.ASM479886v1.pep.all.fa.gz
wget https://ftp.ensembl.org/pub/release-115/fasta/podarcis_muralis/pep/Podarcis_muralis.PodMur_1.0.pep.all.fa.gz

#chicken
wget https://ftp.ensembl.org/pub/release-115/fasta/gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz

#mouse
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz

#merge those proteins into a single file
zcat *.gz > ProteinDB.fasta
```

## Merge the transcript-derived proteins to the proteins available to other species
```
cat CodAn_output/PEP_sequences.fasta ProteinDB.fasta > FINALProteinDB.fasta
```

## Perform gene annotation using GALBA
```
galba.pl --threads 20 --genome=Binsularis_primary_chromosomes.softmasked.fasta --prot_seq=FINALProteinDB.fasta
```
 - It will take a while to finish.
