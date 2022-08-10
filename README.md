# MHCtyping 

### Backgroud:

The Major Histocompatibility Complex (MHC) locus contains many genes associated with antigen presentation, including the ‘classical’ MHCI and MHCII genes, which encode molecules that bind and present peptide fragments to CD8+ and CD4+ T-cells respectively. As key regulators of the antigen-specificity of T-cell immune responses, characterization of the MHCI and MHCII genes is important for many aspects of immunology and associated disciplines such as vaccine development.


In cattle the MHC locus and the repertoire of MHCI/MHCII (BoLA-I/BoLA-II) genes remains only partially characterized. Based on sequence analyses it had been proposed that there were 6 classical BoLA-I loci in cattle, however the number of loci expressed varies between different haplotypes, with different permutations of the six loci represented in different haplotypes. 
Cattle BoLA-DRA is monogenic and although there are 3 BoLA-DRB loci, only BoLA-DRB3 is considered to be functional. As in other species, from a functional perspective BoLA-DR is essentially mono-morphic, whereas BoLA-DRB3 shows high levels of polymorphism with 384 alleles recorded in the [IPD database](https://www.ebi.ac.uk/ipd/mhc/group/BoLA). The genetic organization of BoLA-DQ is more complex, a maximum of 2 loci are expressed in any single haplotype. As with MHCI there is also variability in the number of loci expressed in different haplotypes - with approximately half of the haplotypes expressing single DQA/DQB genes, and the other half expressing 2. As a consequence of intra-and inter-haplotyping pairing of DQA and DQB molecules the presence of duplicate loci enables the generation of multiple DQ molecules.


Due to the polygenic and polymorphic nature of bovine MHCI, DRB3, DQA and DQB, generation of high resolution sequence data using Sanger sequencing has required costly and laborious sub-cloning procedures which has limited large-scale studies of the allelic repertoire and also studies that attempt to sequence complete MHCI/MHCII haplotypes. 


This bioinformatic pipeline is developed for high-throughput NGS bovine MHC genotyping to facilitate a rapid and cost-effective way to characterise the MHCI repertoires of different cattle populations. The sets of novel ‘pan-MHC’ primers that would allow amplification and sequencing of bovine MHCI, DRB3, DQA and DQB alleles using the Illumina MiSeq platform were used. These snakemake pipeline and interative perl scripts is developed to comprehensively analyse the MHC amplicon sequencing data or any similar custom amplicon sequencing data.


### Workflow
---
<img width="490" alt="Screenshot 2022-08-10 at 14 02 11" src="https://user-images.githubusercontent.com/8590103/183908092-7f16b969-ca8a-473e-be77-e8087aaf688a.png">

As shows in a workflow figure, it is split into two seperate analysis steps:

## 1. Snakemake pipeline (steps in green backgroud)
This pipeline analyse the raw sequencing data for each sample and produces the final results for each primer pairs used in the samples.

   - It starts with quality trimming of raw data. 
   - The high quality trimmed reads are then overlapped to generate extedened amplicon sequences.
   - Each extended sequence is then searched for PCR primers used. It searched for the forward and reverse primers, remove the excat primer sequences from reads to remove PCR sites. 
   - 

### Data received from MiSeq (FASTQ):

We receive data for each individual index and the folder names are given sample/library ID in samplesheet spreadsheet during sample submission. Each sample have two sequencing read files – forward and reverse. 

All sequencing facility has their own format of data files and folder structure. Data we receive from MiSeq through Edinburgh Genomics have two read files (fastq.gz) – forward and reverse and text file with the information of both reads (.count file)

  1.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.count*
  2.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.gz*
  3.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_2.fastq.count*
  4.	*180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_2.fastq.gz*

If you want to check what these files have, simply use command less <file> command.
For example, `less 180831_M05898_0019_000000000-BYR6F_1_11441CT0099L01_1.fastq.count`

If you are using same index for multiple samples (like we did with Pig, Cattle and Sheep recently), we must make sure we assign samples with the correct files in the config file (see below).

Raw reads files must end with *_1.fastq.gz* and *_2.fastq.gz* as my script recognises forward and reverse read by this extension. 

### Conda environment:

We need following packages installed in environment: Perl, Blast, Flash, Sickle. 
There are two ways to setup conda environment. 1) using yml file 2) manully installing all packages. 
1. Using yml file: I have attached yml file. Following command will help to install from this file: `conda env create -f environment.yml`
2. Manually creating env and installing packages using following commands. “mhc” is the name of the environment, you can choice any name.
```
conda create --name mhc
proceed ([y]/n)?
y
source activate mhc
conda install -c conda-forge perl
conda install -c bioconda blast
conda install -c bioconda flash
conda install -c bioconda sickle-trim
conda install -c bioconda snakemake
```
Everytime when you want to run scripts or pipeline, use following command to activate environment first: 
`source activate mhc`

  
### Folders and files: 
  
Create a main folder with project name and in this project folder should have following folders:
- **fastq**: It should have all sample folders with paired end data
- **fasta**: 
    - It should have all primer sequences in fasta format for all primer pairs used and indexed database.
    - To index database, use this command: `makeblastdb -in <database.fa> -dbtype nucl`
    - Primer sequences in fasta format. Make sure forward primers have “for” and reverse primers have “rev” in header. It doesn’t matter what else they are called in header but for and rev should be there. 
    - Make sure you remove all ambiguous IUPAC nucleotide and put all possible primer sequences in fasta file. Check out the primer sequences in fasta folder.
- **results**: This directory will have all output files generated by pipeline
- **summary**: This folder will be used to create summary files for each sample
- **scripts**: All the perl scripts that will be used by Snakemake and other steps will be here
- **Config.yaml**: It should be in your main project folder. Check out the attached example config file. You can modify it according to data you are analysing.
    - Make sure you don’t use sample names starts with numeric digits. Sample names must start with alphabets. If sample names start with digits, the simple trick is to put an alphabet in front of all samples. That will work.
- **Snakefile**: It should be in your main project folder.
- **samplesheet.txt**: This file is the list of samples.

### Running pipeline:
Before running pipeline, make sure you have following steps done.
  * In main folder, all subfolders are present.
  * In main folder, Snakefile and config file are present.
  * Config file must have correct information in correct format. Use my test config file to cross check the format. 
  * All reads should be in separate sample folders in fastq folder
  * Database and primer sequences should be in fasta format in fasta folder
  * Primer sequences should only have ATGC sequences and forward and reverse should be names as for and rev
  * Sequence database must be indexed. 
  * All the perl scripts must be in scripts folder

Once above steps are done, follow the commands below:

  1. To check if snakemake pipeline is ready to run and how many jobs will be run, use this command:
         `snakemake -p -n`
  2. To run pipeline
         `snakemake -p -j 3`
      
      Depending on the how many cores/processors you can use on system, change -j. 
      Alternatavely, you can submit individual jobs for each sample on cluster (Check out the script *submit_snakemake_eddie.pl*
  
  It will take a while to finish running pipeline depending on number of clusters you want to check for artefacts in config file and number of samples you are analysing. 
  
  Once the pipeline is finished running, there will be many files presence in result folder for each sample. These outputs are are created at each steps of the pipeline: 
  
## 2. Multiple scripts (Stpes shown in yellow backgroud)

  ----
## Citations
Vasoya D, Law A, Motta P, Yu M, Muwonge A, Cook E, Li X, Bryson K, MacCallam A, Sitt T, Toye P, Bronsvoort B, Watson M, Morrison WI, Connelley T. Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations using next-generation sequencing. Immunogenetics. 2016 Nov;68(10):765-781. [doi: 10.1007/s00251-016-0945-7](https://link.springer.com/article/10.1007/s00251-016-0945-7). Epub 2016 Aug 11. PMID: [27516207](https://link.springer.com/article/10.1007/s00251-016-0945-7)

Vasoya D, Oliveira PS, Muriel LA, Tzelos T, Vrettou C, Morrison WI, de Miranda Santos IKF, Connelley T. High throughput analysis of MHC-I and MHC-DR diversity of Brazilian cattle populations. HLA. 2021 Aug;98(2):93-113. [doi: 10.1111/tan.14339](https://onlinelibrary.wiley.com/doi/10.1111/tan.14339). Epub 2021 Jun 17. PMID: [34102036](https://onlinelibrary.wiley.com/doi/10.1111/tan.14339).

