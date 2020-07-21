# ssGBLUP with R
## General info 
This is a source code for my Bachelor's thesis. The thesis shows testing of the single-step model. 
This method is characterized by simplicity
and equally good accuracy, as well as for example, the two-step method. With one equation,
we can easily obtain breeding values. This model was tested with real data, which included the
genetic marker SNP(Single Nucleotide Polymorphism) and two workability traits describing
milking speed and temperament. This data contains chromosome 29 of Holesteins, which has
significant genetic potential. That means it can suit well to the genomic selection. The
calculation was performed on the workability traits, which play a very important role in the
selection. It’s because they can greatly increase the profitability of the herd. Due to the
one-step methodology and back-solve method, which from the estimated breeding values can
estimate effect of every marker, I was able to identify the location of QTL connected with
functional characteristics. That means I could identify locations, which have a significant
impact on tested quantitative features. I have analyzed two workability traits: the milking
speed and bull temperament. To summerise: In my thesis, I managed to identify QTL locations
on chromosome 29 for bull’s temperament only.
## Data description
The dataset contains genomic and phenotypic data of 7,646 Holstein's Friesian cattle. 
### Phenotopic data
I've analyzed two workability traits:
* Milking speed
* Milking temperament 
The "plots" folder contains histogram of data's distribution, and boxplot showing basic statistic.
### Genotopic data
The genomic data concern the 29 chromosome. This chromosome is important, because the QTL sites turned out to be the most important on this chromosome for given workability traits. Genomic markers were obtained using the Illumina BovineSNP50 BeadChip Version 2 microarray. The maximum number of SNP per chromosome for one individual is 951. Missing data is 0.007%.
## Technologies
* R
* C++

