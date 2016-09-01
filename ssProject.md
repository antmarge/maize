# Project To Do List

### Signatures of Selection Project (Rafael V)

- [ ] Create signature of selection LD/Haplotype [figure](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006178)
- [x] Plot qvalues
- [x] Plot LD blocks from Haploview
- [ ] Sliding window (7 SNPs) for focal SNP 
- [x] Plot difference in lines (Parental and DH lines) for allele frequencies
- [ ] Calculate theta threshold as done in Muir 2016 Plos Genetics

Update from meeting with Dr. Muir

- [ ] Sliding window (7 SNPs) for focal SNP -- Fst
- [ ] Pval for fst then qval the -log(qval)
- [ ] Running average for allele freq
- [x] Redo Haploview haplotype blocks with HOM (remove HET in Tassel)
    

### Inflorescence Tassel GWAS Project

- [x] Get genotypes from Panzea and split by chr
- [IP] Run CHR hapmap files through hapmap2Matrix.R and snpQC() to get format for NAM
- [ ] Population structure analysis (with Travis' code)
- [ ] Run GWAS in NAM with imputed genotype matrix and population structure and BLUPs
