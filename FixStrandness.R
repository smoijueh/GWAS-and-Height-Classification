## How to use this R script


## FUNCTION DEFINITIONS ###
# First complie all the functions below. Then move to the main section
# to implement your data.

## FUNCTION get major and minor allels for opensp user, negative strand of DNA
################# Correct Standedness ###################

get_osnp_MAF<- function(myMatrix){
  MAFs<-vector()
  for (k in 1:(ncol(myMatrix)/2)){ 
    loci<- table(SVM[209:521,c(2*k -1, 2*k)])
    total_alleles<- sum(loci)
    loci<- loci/total_alleles
    
    if (loci[1] > loci[2]){
      major.allele<- names(loci[1])
      minor.allele<- names(loci[2])
    } else{
      major.allele<- names(loci[2])
      minor.allele<- names(loci[1])
    }
    
    MAFs<-cbind(MAFs,cbind(major.allele,minor.allele))
  }
  return(MAFs)
}

#####################################################################


## FUNCTION
## copy matrix and modify incorrect strandness in HPGP Dataset #####

correct_strandness<-function(badMatrix,osnp_matrix){
  myMatrix<-badMatrix
  ## Make haash to find complement of DNA base
  complement_hash = new.env(hash=TRUE, parent=emptyenv(), size=100L)
  assign('A', 'T', complement_hash); assign('T', 'A', complement_hash);
  assign('G', 'C', complement_hash); assign('C', 'G', complement_hash);
  #assign('-', 'NA', complement_hash);
  
  for (snp in 1:(ncol(myMatrix)/2)){ 
    #loci<- table(SVM[1:208,c(2*snp -1, 2*snp)])
    allele1.index = 2*snp-1
    allele2.index = 2*snp
    osnp_loci<-c(osnp_matrix[,allele1.index],osnp_matrix[,allele2.index])
    major.allele<-osnp_loci[[1]]
    minor.allele<-osnp_loci[[2]] 
    
    ##Loop through this single snp(2 columns) and fix bad mappings
    for (user_row in 1:208){
      allele1<-myMatrix[user_row,allele1.index]
      allele2<-myMatrix[user_row,allele2.index]
      
      ## Map missing values to heterozygous
      if(allele1 == '-' || allele2 == '-'){
        myMatrix[user_row,allele1.index]<-major.allele
        myMatrix[user_row,allele2.index]<-minor.allele
        allele1<-major.allele
        allele2<-minor.allele
      }
      ## Allele1
      
      if(allele1 != major.allele && allele1 != minor.allele ){ 
        myMatrix[user_row,allele1.index]<-get(allele1, complement_hash)
      }
      ## Allele2
      
      if(allele2 != major.allele && allele2 != minor.allele ){ 
        myMatrix[user_row,allele2.index]<-get(allele2, complement_hash)
      }
      
    }
  }
  return(myMatrix)
}



#############################################################################


### Helper Function for Finding the Wildtype/Mutant Allele for Each SNP 
compute_MAF<- function(k,alleleMatrix=SVM){
  loci<- table(alleleMatrix[,c(2*k -1, 2*k)])
  total_alleles<- sum(loci)
  loci<- loci/total_alleles
  
  if (loci[1] > loci[2]){
    major.allele<- names(loci[1])
    minor.allele<- names(loci[2])
  } else{
    major.allele<- names(loci[2])
    minor.allele<- names(loci[1])
  }
  return(cbind(major.allele,minor.allele))
}

##########################################################################
###############################   MAIN    ################################

##1 Read in geno.data matrix that combines HPGP and OSNP datasets
## HPGP row 1-208, OSNP rows 209-159
SVM <- as.matrix(read.csv("path/to/genotypes_HPGP_and_OSNP.txt", header=FALSE))
SVM2 <- SVM[,-1]
rownames(SVM2) <- SVM[,1]
SVM<-SVM2; rm(SVM2)

## 2. Run Function to pull out wt/mut alleles in OSNP
## OSNP is our correct reference of alleles mapped to the neg strand
MAF_vector<-get_osnp_MAF(SVM)


## 3. Run Function. Pass in confused matrix(HPGP) and vector (OSNP) 
##    and funct returns matrix mapped all mapped to the same strand.
SVM_fixed<-correct_strandness(SVM,MAF_vector)

## Run GWAS BELOW