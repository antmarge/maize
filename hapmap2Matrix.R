#Margaret Antonio 16.08.30
#Adapted from code by Travis Beckett

#install.packages("NAM")
library(NAM)

# Loading the datasets
# Panzea imputed genotypes H5 file was imported into Tassel5.
# Filtering: 
    #(a) Data->Homozygotes 
    #(b) Filter->Sites: Min count is 80% of genotypes, MAF>.05, Remove minor SNP states
# Data-> Separate (split by chr). Remove chr0 (non-matching)
# Export data as Hapmap WITHOUT ANNOTATIONS, NOT DIPLOID

hapmap2Matrix<-function(chr){
      
      #FUTURE EDITS: separate snpQC and hapmap2
      #OPTION1: Use a loop to got through all files in directory 
      #for (i in 1:chr){
      #Hardcoded to go through an entire directory where all filenames are identical
      #except for chromosome number
      
      #OPTION2: Want to run chromosomes on separate cores so function just takes chr#
      #Note that all files in the directory may only differ by chr#, or i
      i=chr
      
      input=paste0("panzea_imp_chr/panzea_imp_2107genos_c",i,".hmp.txt")
      
      print(paste0("Reading in Hapmap file: ,",input))
      hap <- read.delim(input)

      # Store first 4 columns: chromosome, position, and allele variants at each marker
      marker_info <- hap[,1:4] 
      # Create a vector of alleles 
      alleles=as.character(hap[,2]) 
      # Remove any forward slashes
      alleles<-gsub('/','',alleles)
      
      #Make markers rownames
      rownames(hap)=hap[,1] 
      hap<-hap[,-1]
      #View(hap)
      

      #Transpose the genotype data so m x n matrix is geno x marker
      #Do not include columns 1 to 10 because they don't have genotypes
      print ("Transposing matrix so rows=genotypes, cols=markers")
      geno<-t(hap[,-c(1:10)]) 
      colnames(geno)<-rownames(hap)
      
      #View(genos)
      
      # Fixing missing values
      geno[geno=="N"]<-NA
      
      # Creates a matrix 'geno' from elements of Gen.
      m=nrow(geno) 
      n=ncol(geno)
      genbin=matrix(NA,m,n)
      
      # Loop that recodes genotypes as numbers 
      # 0 is missing, 1 is HET, 2 is HOM
      print ("Recoding genotype data to 0,1,2")
      for(i in 1:n){ 
        A1=strsplit(alleles[i],'')[[1]][1]
        A2=strsplit(alleles[i],'')[[1]][2]
        BB=paste(A1)
        bb=paste(A2)
        M=as.character(geno[,i])
        genbin[M==BB,i]=2
        genbin[M==bb,i]=0
      }
      
      colnames(genbin)=paste(hap[,3],hap[,4],hap[,2],sep='.')
      #colnames(genbin)<-colnames(geno)
      rownames(genbin)<-rownames(geno)
      matrixfile=paste0("genoMatrix/chr",i,"_genoMatrix.txt")
      print(paste0("Writing matrix (not imputed by NAM) to: ",matrixfile))
      write.table(genbin,matrixfile,row.names=TRUE,quote=FALSE,append=FALSE)
      
      # See the package 'NAM' in R for details about the functions
      igenbin <- snpQC(genbin,psy=1,MAF=0.05,remove=TRUE,impute=TRUE) 
      
      # Write NAM imputed matrix to a new file
      newfile=paste0("namImp/chr",i,"_namImp.txt")
      print(paste0("Outputting NAM imputed matrix to file: ",newfile))
      write.table(igenbin,newfile,row.names = TRUE,quote=FALSE,append=FALSE)
  }





