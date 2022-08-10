#MHCtyping 

###Backgroud:
The Major Histocompatibility Complex (MHC) locus contains many genes associated with antigen presentation, including the ‘classical’ MHCI and MHCII genes, which encode molecules that bind and present peptide fragments to CD8+ and CD4+ T-cells respectively. As key regulators of the antigen-specificity of T-cell immune responses, characterization of the MHCI and MHCII genes is important for many aspects of immunology and associated disciplines such as vaccine development.

In cattle the MHC locus and the repertoire of MHCI/MHCII (BoLA-I/BoLA-II) genes remains only partially characterized. Based on sequence analyses it had been proposed that there were 6 classical BoLA-I loci in cattle, however the number of loci expressed varies between different haplotypes, with different permutations of the six loci represented in different haplotypes. 
Cattle BoLA-DRA is monogenic and although there are 3 BoLA-DRB loci, only BoLA-DRB3 is considered to be functional. As in other species, from a functional perspective BoLA-DR is essentially mono-morphic, whereas BoLA-DRB3 shows high levels of polymorphism with 384 alleles recorded in the IPD database (https://www.ebi.ac.uk/ipd/mhc/group/BoLA). The genetic organization of BoLA-DQ is more complex, a maximum of 2 loci are expressed in any single haplotype. As with MHCI there is also variability in the number of loci expressed in different haplotypes - with approximately half of the haplotypes expressing single DQA/DQB genes, and the other half expressing 2. As a consequence of intra-and inter-haplotyping pairing of DQA and DQB molecules the presence of duplicate loci enables the generation of multiple DQ molecules.

Due to the polygenic and polymorphic nature of bovine MHCI, DRB3, DQA and DQB, generation of high resolution sequence data using Sanger sequencing has required costly and laborious sub-cloning procedures which has limited large-scale studies of the allelic repertoire and also studies that attempt to sequence complete MHCI/MHCII haplotypes. 

This bioinformatic pipeline is developed for high-throughput NGS bovine MHC genotyping to facilitate a rapid and cost-effective way to characterise the MHCI repertoires of different cattle populations. The sets of novel ‘pan-MHC’ primers that would allow amplification and sequencing of bovine MHCI, DRB3, DQA and DQB alleles using the Illumina MiSeq platform were used. These snakemake pipeline and interative perl scripts is developed to comprehensively analyse the MHC amplicon sequencing data or any similar custom amplicon sequencing data.


###Workflow
<img width="490" alt="Screenshot 2022-08-10 at 14 02 11" src="https://user-images.githubusercontent.com/8590103/183908092-7f16b969-ca8a-473e-be77-e8087aaf688a.png">




Deepali Vasoya, Andy Law, Paolo Motta, Mingyan Yu, Adrian Muwonge, Elizabeth Cook, Xiaoying Li, Karen Bryson, Amanda MacCallam, Tatjana Sitt, PhilipToye, Barend Bronsvoort, Mick Watson, W. Ivan Morrison and Timothy Connelley. "Rapid identification of bovine MHCI haplotypes in genetically divergent cattle populations Using Next-Generation Sequencing." PMID: [27516207] (https://link.springer.com/article/10.1007/s00251-016-0945-7)
