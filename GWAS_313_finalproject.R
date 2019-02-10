# Samuel Moijueh
# BS 845Applied Statistical Modeling and Programming in R
# Fall 2014 Final Project
# December 8, 2014

library(genetics)
library(parallel)
library(qqman)

setwd("~/Desktop/GWAS/")

################# Function Definitions ###################

#' compute the Major and Minor Allele
compute_MAF<- function(j){
  loci<- NULL
  loci<- table(as.vector(geno.data[,c(2*j -1, 2*j)]))
  total_alleles<- sum(loci)
  loci<- loci/total_alleles
  
  # major and minor allele frequencies for one locus
  if (loci[[1]] == 0.5 & length(loci) == 2){
    major_allele<-names(loci[1]); minor_allele<-names(loci[2])
  } else {
    major_allele<-names(which.max(loci)); minor_allele<-names(which.min(loci))
  }
  return(cbind(major_allele,minor_allele))
}

#' Vectorized: perform a Hardy Weinburg Equlibrium Chi Square Test for a given SNP and return its p-value
get.HWEpval <- function(count) {
  g2 <- genotype(paste0(geno.data[,2*count-1], "/", geno.data[,2*count]))
  if (all(homozygote(g2))){
    return(1);
  } else
    HWE.chisq(g2)[3]   # calculate the HWE pvalue for one SNP
}

#' SNP to geno.data index
getGenomeIndex<-function(x){
  c(x * 2 - 1, x * 2)
}

############# MAIN: Start Here ############

############ Read in/Load the GWAS data files ###################
# loading the data files into R takes 40 seconds
load("~/Desktop/GWAS/import_files.RData")

# geno.data is the SNP data. each row is an individual, every pair of columns is a SNP
# pheno.data is the height of each individual in centimeters
# chromosome_bin is a vector of all the SNPs indexed by chromosome. It is used later in generating the 
# the Manhattan Plot

######## Examine the Phenotype Distribution and Test for Normality #########
pd <- pheno.data$V1
mean(pd)
var(pd)

hist(pd, xlab = "Height in Centimeters of N = 313 individuals",
main = "Distribution of Heights", prob = TRUE, col = "pink")

# normal density curve
lines(density(pd, adjust = 2), col="blue", lwd=2)

# kernel density curve
lines(density(pd), col="darkgreen", lty="dotted", lwd=2)

legend("topleft", inset=.05, title="Density Curves",
       c("Normal Density Curve", "Kernel Density Curve"), fill=c("blue", "darkgreen"), horiz=FALSE)

par(mfrow=c(1,2))
distn <- rnorm(n = 1000,mean = mean(pd),sd = sd(pd))

# density curve from the observed data
plot(density((pd)), main = "Theoretical Density (red) \n vs Empirical Kernel Density")
lines(density(distn),col="red")  # density curve from the simulated data
legend("bottom", inset=.05, title="Density Curves",
       c("Theoretical Density", "Empirical Kernel Density"), fill=c("red", "black"), horiz=FALSE)

# test normality with qq plot
qqnorm(pd)
qqline(pd)  # the phenotype is normally distributed
par(mfrow=c(1,1))

t.test(distn,pd)
# p-value is approx. 0.7   # data most likely came from the same distribution


########### Genome Filtering Protocol ##########

# need to check HWE (Hardy W-berg Equilibrium: performed on each SNP) to check the quality. SNP specific quality

# 1. Perform HWE chiq test on each SNP to ensure Genome Quality
remove<- c()
Num_of_SNPs<-ncol(geno.data) /2

for (count in 1:Num_of_SNPs){   
  HWE_pvalues[count]<- get.HWEpval(count)   # stores the indices of SNPs that fail HWE  
}

#HWE_pvalues<- unlist(mclapply(1:Num_of_SNPs, function(x) get.HWEpval(x), mc.preschedule = TRUE, mc.cores = 4))

save.image("HWE_pvalues_backup.RData")  # save as a checkpoint

# 2. store the indices of SNPs that fail HWE chi square test
remove<- sort(getGenomeIndex(which(HWE_pvalues<10^-4)))  # remove is 102788

# 3. Remove SNPs that has a Minor Allele Frequency less than 6%
for (j in seq(1, Num_of_SNPs, 2)){
  loci<- table(geno.data[,c(2*j -1, 2*j)])
  total_alleles<- sum(loci)
  loci<- loci/total_alleles
  
  if (min(loci) < 0.06){
    remove<- c(remove, c(j,j+1))   # remove is 18834; these are SNPs to be removed 
  }
}

remove<-sort(unique(remove))    # remove is 20964

indices <- sort(remove[c(F,T)])
indices <- indices/2

save.image("remove_vector.RData")  # another checkpoint in the code

# 4. filter bad SNPs from the geno.data and the Chromosome Bin
if (!is.null(remove)){
  # filter the genotypes
  geno.data<- geno.data[, -remove] # 10482 bad SNPs removed
}

chromosome_bin<- chromosome_bin[-indices,]

# 5. Person-wide missing rate: Remove SNPs with missing data in more than 5% of individuals
badrows<- which( rowSums(geno.data=="--")/nrow(geno.data) >= 0.05 )   # faster 

# 6 . SNP-wide missing rate: Remove SNPs that missing data in more than 5% of the of the SNP
badcols <- which( colSums(geno.data=="--")/ncol(geno.data) >= 0.05 )   # faster

if (length(badcols) != 0){
geno.data<- geno.data[, -badcols]
}

if (length(badrows) != 0){
geno.data<- geno.data[-badrows,]
}

############ Anova genotype coding, Creating EIGENSTRAT file ##############

# Anova Dummy Coding 
# Xa is the additive score; Xd is the dominant score

# matrix for the coded genotypes
Num_of_SNPs<-ncol(geno.data) /2
Xa <- matrix(NA, nrow = nrow(geno.data), ncol = Num_of_SNPs)

EIGENSTRAT_genofile<-vector("list",Num_of_SNPs)

# transforms the genetic loci into numeric values and create EIGENSTRAT genofile.
for (j in 1: Num_of_SNPs){
  
  allele<-compute_MAF(j)
  
  X <- 1 * (geno.data[,c(2*j -1, 2*j)] == allele[2])  # minor allele
  reference_allele_count <- rowSums(geno.data[,c(2*j-1,2*j)]==allele[1], na.rm=TRUE)
  EIGENSTRAT_genofile[[j]]<-reference_allele_count
  
  Xa[,j] <- X[,1] + X[,2] - 1
}

Xd <- 1 - 2*abs(Xa)

EIGENSTRAT_genofile<-do.call(rbind, EIGENSTRAT_genofile)
rownames(EIGENSTRAT_genofile)<- NULL
write.table(EIGENSTRAT_genofile, file="eigenstrat_genofile.txt", row.names=FALSE, col.names=FALSE, sep ="")

############## Midpoint of Analysis ####################
save.image("analysis_midpoint.RData") # loading takes 1 minute

PCA <- as.matrix(read.table("First_TENPC.txt", sep=" "))

### Perform multiple LMs to obtain p-values for each SNP in the GWAS ####
pval <- vector(length = Num_of_SNPs)
m2<-lm(pheno.data$V1 ~ PCA[,5])
for (i in 1 : Num_of_SNPs){
  X.mx <- cbind(Xa[,i],Xd[,i])

  m1<-lm(pheno.data$V1 ~ X.mx + PCA[,5])
  
  pval[i]<- anova(m2,m1)$Pr[2]
}

# test<- lm(pheno.data$V1 ~ X.mx + PCA[,1] + PCA[,2] + PCA[,3] + PCA[,4] + PCA[,5] + PCA[,6] + PCA[,7] + PCA[,8] + PCA[,9] + PCA[,10])

####### Calculate Genomic inflation factor, Lambda #####
chisq <- qchisq(1-pval,1)
lambda = median(na.omit(chisq))/qchisq(0.5,1)

###### Manhanttan and QQ Plots  ##########

# this loads the data (variables and data frames) from the preceding code
load("OH_obtained_pvalues.RData")

pval_manhattan<- cbind(pval,chromosome_bin)
manhat<-data.frame(SNP = pval_manhattan[,2], CHR = as.integer(pval_manhattan[,3]),  
                   BP = as.integer(pval_manhattan[,4]), P = as.numeric(pval_manhattan[,1]), stringsAsFactors=FALSE)

qq(manhat$P, main = "Q-Q plot of GWAS p-values")
#text(locator(1), paste("Lamdba = ",lambda))

manhattan(na.omit(manhat), col = c("black", "dimgrey", "chocolate3"), main = "Manhattan Plot of Human Height")

top_most_SNPs<-manhat[with(manhat, order(manhat$P)),]  # sort in ascending order of p-values
top_500_SNPs<-head(top_most_SNPs, 500)
the_SNPs<-top_500_SNPs$SNP

write.csv(manhat, file = "GWAS_results.csv")

#### Identify Significant SNPs in the GWAS ########

# # indices of p-values which are above the suggestive line
# sig_marker=which(-log10(manhat_log$P) > -log10(1e-05))  
# 
# # markers<- levels(droplevels(manhat[sig_marker,1]))
# markers<- c()
# for (i in 1:length(sig_marker)){
#   markers[i]<-manhat[sig_marker[i],]$SNP
# }
# 
# 
# # highlight the significant snps on the Manhattan Plot
# manhattan(manhat_log, col = c("black", "dimgrey", "chocolate3"), main = "Human Height", highlight = markers)