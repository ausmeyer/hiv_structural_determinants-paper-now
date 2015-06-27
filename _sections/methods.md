---
title: "Methods"
order: 3
---

### Obtaining and preparing HIV gene sequences

Gene sequences came from the HIV database at Los Alamos National Laboratory ([HIV Database](http://www.hiv.lanl.gov/content/index)). We used only pre-made gene alignments. The alignments were assumed correct and therefore not changed for all of the genes in this study. To identify the genes in the whole genome sequences, we used the annotation landmarks available in the sequence database ([HIV Landmarks](http://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html)). With [this script](data/gp120/gp120_sequences/translate_dna_sequences.py), sequences were filtered to include only those with entirely canonical bases, and to change stop codons in the original alignment to gaps. We used Biopython for all sequence manipulation and cleaning ([Cock et al. 2009](https://dx.doi.org/10.1093/bioinformatics/btp163)).

### Calculating the site-wise evolutionary rate ratio

We built phylogenetic trees with FastTree 2.1.7 with the following input line ([Price et al. 2009](https://dx.doi.org/10.1371/journal.pone.0009490)).

<pre> <span style="font-family:Courier">fasttree -nt -gtr -nosupport alignment.fasta > alignment.tree</span> </pre>

In addition, we made minor changes to the FastTree source code to prevent rounding errors in the tree branch lengths (the branch length issue was detailed  [here](http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)).

Given a codon alignment and evolutionary tree, we used the phylogenetic software HyPhy to calculate the evolutionary rate ratio $$ dN/dS $$ ([Kosakovsky-Pond et al. 2005](https://dx.doi.org/10.1093/bioinformatics/bti079)). We used the built-in one-rate fixed effects likelihood (FEL) method with the default model ([Muse 1994](http://www.ncbi.nlm.nih.gov/pubmed/7968485)). The FEL method calculates an independent $$ dN/dS $$ value at each column in the alignment. To calculate a one-rate $$ dN/dS $$ value, the model computes a single $$ dS/dt $$ value for the entire alignment and an individual $$ dN/dt $$ at each site. It then normalizes each $$ dN/dt $$ value by dividing by the average $$ dS/dt $$ to obtain $$ dN/dS $$ ([Yang 2006](http://download.bioon.com.cn/view/upload/month_0808/20080811_979a0656719f9157b466IruZP7mlp2U9.attach.pdf)). 

### Obtaining and mapping protein structures

Protein structures were obtained from the RCSB Protein Databank (PDB) ([Berman et al. 2000](https://dx.doi.org/10.1093/nar/28.1.235), [Berman 2008](https://dx.doi.org/10.1107/S0108767307035623)). For each of the six proteins, we obtained both the monomeric structure and, if quaternary structure was functionally important, the full biological assembly. This study included all three of the enzymes encoded by the HIV *pol* gene: reverse transcriptase (PDBID: 1HYS, [Sarafianos et al. 2001](https://dx.doi.org/10.1093/emboj/20.6.1449)), integrase (PDBID: 1EX4 and 1K6Y, [Chen et al. 2000](https://dx.doi.org/10.1073/pnas.150220297) and [Wang et al. 2001](https://dx.doi.org/10.1093/emboj/20.24.7333)), and protease (PDBID: 1HPV, [Kim et al. 2000](https://dx.doi.org/10.1021/ja00108a056)). In addition we included the three HIV structural proteins that did not bind nucleic acid. Those included the capsid (*p24*, PDBID: 3H47), the matrix (*p17*, PDBID: 1HIW, [Hill et al.](http://www.ncbi.nlm.nih.gov/pubmed/8610175)), and the major component of the HIV receptor binding protein (gp120, PDBID: 4TVP, [Pancera et al. 2014](https://dx.doi.org/10.1038/nature13808)). We broadly followed the HIV structure suggestions curated by the PDB ([HIV Structures](http://www.rcsb.org/pdb/education_discussion/educational_resources/struct_bio_hiv_lores.pdf)).

Integrase constituted a special case as it lacked a full length structure in a single experiment and it was the protein and showed no structural predictors. As a result, we tried a number of different structures to find structural predictors. We first tried the 1EX4 structure alone; this structure include the C-terminal domain, but lacks the N-terminal region. In addition, we tried using the 1K6Y structure alone which includes the N-terminal region, but lacks a portion of the C-terminal domain. Then, in an attempt to have a full length model, we aligned the center structural domain of 1EX4 with that of 1K6Y. Fortunately, there is a relatively large gap from site 46 to site 56 in the N-terminal domain of 1K6Y and the central domain of the two structures align with $$ RMSD = 0.686 $$. Therefore, the core of the two proteins aligned almost perfectly making it possible to simply translate and add the coordinates of the N-terminal domain atoms from 1K6Y to the more complete 1EX4 structure. This procedure produced a nearly full length model of integrase that was used to construct the further models.

To map protein structures onto the existing alignment, we used a developmental version of the sequence alignment software MAFFT ([Katoh and Frith 2012](https://dx.doi.org/10.1093/bioinformatics/bts578), [Katoh and Standley 2013](https://dx.doi.org/10.1093/molbev/mst010), [Katoh and Standley 2014](https://dx.doi.org/10.1007/978-1-62703-646-7_8)). This version included the ability to add a sequence to an existing alignment while removing any sites where there was an insertion in the protein structure sequence. Therefore, we used the following input line.

<pre> <span style="font-family:Courier">mafft --mapout --addfragments aas.fasta alignment.fasta > added_alignment.fasta</span> </pre>

### Calculating structural predictors

We used two independent structural predictors from the protein structures. The first was relative solvent accessibility (RSA) which has been detailed extensively previously ([Scherrer et al. 2012](https://dx.doi.org/10.1186/1471-2148-12-179), [Meyer and Wilke 2012](http://dx.doi.org/10.1093/molbev/mss217), [Meyer et al. 2013](http://dx.doi.org/10.1098/rstb.2012.0334), [Meyer and Wilke 2015](http://dx.doi.org/10.1371/journal.ppat.1004940)). We used the program DSSP to compute the absolute solvent accessibilities of each residue in the protein ([Kabsch and Sander 1983](http://www.ncbi.nlm.nih.gov/pubmed/6667333)). Then, we used maximum absolute accessibility values for each residue to normalize solvent exposure to a value between 0 and 1 ([Tien et al. 2013](https://dx.doi.org/10.1371/journal.pone.0080635)). For all models in this study, we used the RSA of a chain in the functional multimeric state of the biological assembly. The second metric was distance to a reference point ([Meyer and Wilke 2015](http://dx.doi.org/10.1371/journal.ppat.1004940)). To compute the set of distances, we used each C-alpha in the protein as a reference point, and calculated the distance to every other C-alpha. Thus, the distances from a single C-alpha to every other C-alpha constituted a single distance set. We repeated this calculation using every amino acid to generate a symmetric matrix of distances where the diagonal is always zero as it is simply the distance from an amino acid to itself. For gp120, we included glycosylation sites as a categorical predictor. To find sites covered by sugar moieties, we used the glycosylations found in the 4TVP structure to calculate all amino acid sites that were within four angstroms of any sugar atom. Then, we included the set of these sites are a predictor in the linear models below.

### Constructing linear models and cross validation

It has been shown previously that RSA is a relatively strong structural predictor for site-wise $$ dN/dS $$ both in viruses and large enzyme datasets. Thus, we started with a linear model that predicted site-wise $$ dN/dS $$ with site-wise RSA. Then, to find the best reference point among the entire set of reference points, we constructed combined linear models with both RSA and the distance set from a single reference point. To guard against overfitting, we trained the model on a randomly chosen 75% of the data and reserved 25% for validation of the best reference point found in the training set. To be clear, all of the possible reference points were available for training, but for each training set we only used $$ dN/dS $$, RSA, and distances for 75% of the sites. We repeated the training and validation 100 times for each protein. We found that the training and validation $$ R^{2} $$ and p-value we very similar. Therefore, for each protein, we used the site that was found most frequently during training as the best site for the final prediction. All of our models and figures were generated in the statistical language R utilizing the ggplot2 package ([Ihaka and Gentleman 1996](https://dx.doi.org/10.1080/10618600.1996.10474713), [Wickham 2009](http://had.co.nz/ggplot2/book)).  