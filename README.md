The R markdown script of Solak et al. 2024 research paper. All the custom functions were written by Jakub Kreisinger (see: https://www.researchgate.net/profile/Jakub-Kreisinger). 


The script contains the steps below:

1 - Demultiplexing and Trimming with Skewer [use Bash Terminal]

2 - DADA2 Quality Filtering and Creating OTU Table

2A - DADA2 Quality Filtering and Creating OTU Table using a linux cluster (TRUBA), use the script below [use Bash Terminal,not R]

3 - Dada2 get REF.fasta [use Bash Terminal, not R]

4 - Filter the chimeras with Usearch [use Bash Terminal, not R]

5 - DADA2 Assign Taxonomy [might take long time in PC, better to  use on linux server]

6 - Create phyloseq object

7 - Number of sequences per sample

8 - Compare Duplicates

9 - Merge Duplicates

9A - re-assign taxonomy [not mandatory read the text below]

9B - Number of sequences per sample for the final database

9C - DATEST RELATIVE ABUNDANCE

10 - Barplots

11 - ALPHA and BETA DIVERSITY

12 - MIXED MODELLING

13 - OTU differential abundance testing with DESeq2

14 - Effect of the sampling year, body mass, and sex

15 - Microsatellite Dataset, The correlation between host genetics and microbiome


If you have any problems, contact me on mr.solak@hotmail.com.
