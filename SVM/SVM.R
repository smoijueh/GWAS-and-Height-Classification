# Samuel Moijueh
# Demarcus Briers
# November 26, 2014

setwd("~/Desktop/SVM/")

library(hash)

############################# Import Demarcus' 23andMe SNP genotype data ##############################
demarcus <- read.csv("demarcus.snp", header=FALSE, stringsAsFactors=FALSE)
demarcus<-replace(demarcus,demarcus==TRUE, 'T')

################## Import OpenSNP meta data including Height Class for Every User ######################
load("height_class.Rda")
colnames(master_opensnp_meta) <- c("userID", "height", "height_class")


##################### Import Significant SNP genotype data for Every Individual ########################
SVM <- as.matrix(read.csv("sigSNPsgenodata.txt", header=FALSE))
SVM2 <- SVM[,-1]
rownames(SVM2) <- SVM[,1]
SVM<-SVM2; rm(SVM2)

num_of_SNPs<-ncol(SVM)

############### Subset the OpenSNP meta data to obtain SNP genotypes of Reference Population ###########

ref_pop <- subset(master_opensnp_meta, master_opensnp_meta[,3] == 2)

# make a matrix of the SNP genotype data for each person in the reference population
ref_pop <- as.matrix(read.csv("~/Desktop/SVM/ref_pop.txt", header=FALSE))
ref_pop2<- ref_pop[,-1]
rownames(ref_pop2) <- ref_pop[,1]
ref_pop<-ref_pop2; rm(ref_pop2)

################# Helper Function for Finding the Wildtype/Mutant Allele for Each SNP ##################
compute_MAF<- function(k){
  loci<- table(SVM[,c(2*k -1, 2*k)])
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


###################### Function for Creating the WT/MUT matrix for Each SNP ###########################
create_WT_MUT_matrix<- function(SNP_matrix){
  count<-0
  score_matrix<-matrix(, nrow = nrow(SNP_matrix), ncol = ncol(SNP_matrix))
  h = new.env(hash=TRUE, parent=emptyenv(), size=100L)
  
  for (j in 1:ncol(SNP_matrix)){  
    
    # if j is odd then update count
    if (j %% 2 != 0){
      count<-count + 1
    }
    
    for (i in 1:nrow(SNP_matrix)){
      data<-compute_MAF(count)
      wild_type<- data[1]
      mutant <- data[2] 
      h = new.env(hash=TRUE, parent=emptyenv(), size=100L)
      assign(wild_type, 'W', h)
      assign(mutant, 'M', h)
      score_matrix[i,j] <- get(SNP_matrix[i,j], h)
    }
    print(paste("j = ",j," count = ",count))
  }
  return(score_matrix)
}

########### Helper Function for Calculating Euclidian_distance used in the SVM Scoring Matrix ############

## Score matrix to reperesent differences in alleles between a user and a user the reference population.
M1 = t(matrix(
  c(1,2.5,10,2.5,1,7.5,10,7.5,1),  
  nrow=3, ncol=3)
)

## Hash used to map letter combo to Score Matrix
score_hash = new.env(hash=TRUE, parent=emptyenv(), size=100L)
assign('WW', 1, score_hash); assign('WM', 2, score_hash);
assign('MW', 2, score_hash); assign('MM', 3, score_hash);


# Bootstrapping-median used to obtain score
euclidian_mapping <- function(genotype1,genotype2,snp_number,numSamples = 150){
  ##set hashkey for current user
  genotype_user<-paste0(genotype1,genotype2, sep = "")
  
  # Bootstrapping: select a random index, or random row of user in the reference population.
  rand_refuser_index <- sample(1:nrow(ref_pop_WT_MUT_matrix),numSamples,replace = T)
  
  score<-vector()
  
  for (ref_index in rand_refuser_index ){
  ##set hashkey for random reference user. given snp index get col index using 2n-1,2n
  snp_col1<-(2 *snp_number) - 1
  snp_col2<-snp_col1 + 1
  
  ref_allele_a<-ref_pop_WT_MUT_matrix[ref_index,snp_col1]
  ref_allele_b<-ref_pop_WT_MUT_matrix[ref_index,snp_col2]
  genotype_ref<-paste(ref_allele_a,ref_allele_b,sep = "")
  
  score<-rbind(score,M1[get(genotype_ref,score_hash),get(genotype_user,score_hash)])
  }
  return(median(score))
}


######################### Function for Converting SNP alleles to score matrix ###########################
compute_wt_mutant_scores<-function(user_matrix){
SVM_scores<-matrix(, nrow = nrow(user_matrix), ncol = (ncol(user_matrix)/2))

for (snp_idx in 1:(ncol(user_matrix)/2)){
  for (row in 1:nrow(user_matrix)){
    allele1<-user_matrix[row,(2*snp_idx-1)] #must go from snp idex to allele matrix
    allele2<-user_matrix[row,(2*snp_idx)]
    SVM_scores[row,snp_idx] = euclidian_mapping(allele1,allele2,snp_idx)
    
  }
  print(paste("snp_idx = ",snp_idx))
}
return(SVM_scores)
}

################################################# MAIN #############################################

d<- create_WT_MUT_matrix(demarcus)

wildtype_mutant_matrix<- create_WT_MUT_matrix(SVM)

ref_pop_WT_MUT_matrix<- create_WT_MUT_matrix(ref_pop)

SVM_scores <- compute_wt_mutant_scores(wildtype_mutant_matrix)

write.table(SVM_scores, file="SVM_scores.output", row.names=FALSE, col.names=FALSE)
write.table(master_opensnp_meta[,3], file="height_classes.output", row.names=FALSE, col.names = FALSE)

SVM_scores<-compute_wt_mutant_scores(d)

write.table(SVM_scores, file="demarcus_SVM_scores.output", row.names=FALSE, col.names =FALSE)

