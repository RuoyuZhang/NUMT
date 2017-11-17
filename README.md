# NUMT
A probabilistic approach for mtDNA/NUMTs reads separation in NGS data

The main code is mtReadLikelihood.py. It will take an alignment file (sam format) and calculate the log likelihood of each alignment in the file.

insertMutation.py can help generate mutated genome, which is required to simulate sequencing data.

mt.variants.frequency.txt and numt.variants.frequency.txt are mutation frequency file for mitochondrial genome and NUMT respectively.

chrM.fa and NUMTs.fa are mitochondrial and NUMT reference sequences (from hg19).
