# GWAS and Height Classification
###### Boston University &bull; ENG BE 562: Computational Biology &bull; Samuel Moijueh, Demarcus Briers &bull; 12/2014

---

![height](https://user-images.githubusercontent.com/42754056/53932718-97d70500-4060-11e9-962a-b2a67b59cc27.jpg)

**Read the final paper** on [Google Drive](https://drive.google.com/open?id=1uC6GpkJsjCKZOAEV_LSqSgIYeNIrbNMz)

&nbsp;
## Motivation

Human height is a complex polygenic trait that has been difficult to predict due to the influence of potentially thousands of interacting genetic loci. Genetics accounts for approximately 80% of a person height (the other 20% is nutrition). Despite the overwhelming influence of hereditary factors, there are currently no genetically based models used to predict a person’s height.

Here we present a predictive model (using Genome Wide Association Study and Support Vector Machine) to unravel this genetic enigma. The results of this study are paramount to the practice of pediatric endocrinology and embryonic screening. This novel approach serves as a proof of concept for classification of other polygenetic traits and diseases.

&nbsp;
## What is a Genome Wide Association Study (GWAS)?
GWAS is a biostatistic approach for interrogating SNPs that commonly arise in populations, and determining whether they are associated with a disease or phenotypic trait.

In this project, GWAS was a linear regression of the principal components. (**See Method Overview for more details**)

Geneticists and Epidemiologists have used GWAS to study genetic variations in diseases such as asthma, cancer, diabetes, heart disease, among others.

**The GWAS Workflow is Outlined in the Figure Below**

![gwas_workflow](https://user-images.githubusercontent.com/42754056/53932503-7de8f280-405f-11e9-92b1-0d321dd7802e.png)

&nbsp;
## Method Overview
**Read the [final paper](https://drive.google.com/open?id=1uC6GpkJsjCKZOAEV_LSqSgIYeNIrbNMz) for more details**. Although not explicitly stated here, a multiple hypothesis test correction was applied in key steps below.


### Genomic Association

1. Downloaded SNP Array Data of 313 individuals from [OpenSNP](https://opensnp.org/).
2. Perl was used to parse, clean, and wrangle the raw SNP array data into a genotype matrix
3. Genotype and phenotype matrices were imported into R
4. Checked [Linear Regression Model Assumptions](http://r-statistics.co/Assumptions-of-Linear-Regression.html)
5. Performed computationally intensive filtering protocol ([Hardy Weinberg Equilibrium Chi Square Test](http://www.husdyr.kvl.dk/htm/kc/popgen/genetics/2/2.htm)) to ensure SNP genomic qualiy
6. Alleles were encoded into numeric values based on the major and minor allele frequency
7. Created additive genotype (Xa) and dominant genotype (Xd) matrices as required by [EIGENSTRAT](https://github.com/DReichLab/EIG/tree/master/EIGENSTRAT)
7. EIGENSTRAT was used to perform a Principal Component Analysis which corrects for [population-stratification](https://www.nature.com/articles/ng1847)
8. The first 10 principal components were used to fit a linear regression model to the reported heights of the individuals
9. An [ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance) likelihood test was applied to the fitted linear regression model to obtain p-values for each SNP
10. The p-values were used to generate a QQ-plot and Manhattan Plot. Evaluate Results
11. Obtained a list of the top 500 SNP in order of statistical significance (p-value)
11. The alleles associated with the significant SNPs were mapped to the Euclidian Space and assigned a score (**see Paper for an explanation of this step**). These values were scaled and normalized

### Height Class Stratification

1. [k-means clustering](https://en.wikipedia.org/wiki/K-means_clustering) (k=2 and k=3) was used to assign individuals into 2 (binary) and 3 (multi) height classes.
    - The goal here was to stratify the height classes into relatively equal clusters based on variance

### Support Vector Machine

1. Split data into training and validation set
2. Trained a Support Vector Machine for Classification
    - Used the list of p-prioritized SNPs and allele encoded values to make predictions
3. 10-fold cross validation to identity optimal parameters
4. Evaluate model performance on Test Set

## Results

**The results below reflect the best performing model after a 10-fold cross validation**

* The binary classification (tall vs short) model achieved a predictive accuracy of 86%
* The multi-class classification achieved a predictive accuracy of 72%

## Limitations

* only 313 individuals, low [statistical power](https://en.wikipedia.org/wiki/Power_(statistics))
* unable to account for common covariates; individuals were missing AGE and SEX
    - Height is a sexually dimorphic trait; on average men are taller than women in all human populations
    - knowing this information would contribution to the predictive accuracy

## Future Direction

* obtain SNP array data from **more** individuals
* examine data for [Linkage Disequilibrium (LD)](https://en.wikipedia.org/wiki/Linkage_disequilibrium)
   - Regional Gene Plot to Visualize LD
   - [Haplotype Analysis](https://en.wikipedia.org/wiki/Haplotype)
   - Examine for [Epistasis](https://en.wikipedia.org/wiki/Epistasis)

  ![gwas2](https://user-images.githubusercontent.com/42754056/53936829-d674bb80-4070-11e9-98ec-28260ac5c012.jpg)
