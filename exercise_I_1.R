#Cancer Phylogenies I
#https://moodle.ucl.ac.uk/pluginfile.php/3938207/mod_folder/content/0/CancerPhylogeniesI_practical_solutions_exercise1.html?forcedownload=1

snvs <- read.delim(file.choose()) #open "tumour1.vcf"
head(snvs)
colnames(snvs)

#How many mutations are there in this tumour?
nrow(snvs) #15518

#How many mutations per chromosome can we observe?
table(snvs$CHROM)

#Note that there are some strangely looking chromosomes here (e.g. chrUn_gl000211). 
#These are pseudo-chromosomes that make up all the bits of the reference genome assembly that did not fit within the classical chromosomes. 
#See the explanation from the UCSC website: “ChrUn contains clone contigs that cannot be confidently placed on a specific chromosome. 
#For the chrN_random and chrUn_random files, we essentially just concatenate together all the contigs into short pseudo-chromosomes. 
#The coordinates of these are fairly arbitrary, although the relative positions of the coordinates are good within a contig.”

#Ignore the pseudo-chromosomes for now and remove them from the analysis
grepl("gl000",snvs$CHROM) #從snvs$CHROM篩選有含“gl000”的string
which(!grepl("gl000",snvs$CHROM))
snvs.keep <- snvs[which(!grepl("gl000",snvs$CHROM)),]
table(snvs.keep$CHROM) #only normal chromosomes are left
unique(snvs.keep$CHROM)

nrow(snvs.keep) # mutation down from 15518 to 15510 after removing out mutations(num: 8) at pseudo-chromosomes
table(snvs.keep$FILTER) #Every mutations are labeled as 'pass'

#How many base substitutions of each type are there in this tumour?
snvs.keep$BaseSubstitution <- apply(snvs.keep[,c("REF","ALT")],1,       #"snvs.keep$BaseSubstitution" set new column
                                    function(x) paste0(x[1],">",x[2]))  # apply(data, 1(work on row or 2 for working on column)
colnames(snvs.keep) #new column of "BaseSubstitution" is made
table(snvs.keep$BaseSubstitution) #table of numbers in every different substitution

#Plot a graph of base substitution counts per chromosome:
snvs.keep[,c("CHROM","BaseSubstitution")]
counts <- table(snvs.keep[,c("CHROM","BaseSubstitution")])
counts #二維table!! (counts per chromosome)

library(reshape)
df.counts <- melt(counts) #把25 x 12的表格變成1 x 300
df.counts

# Stacked bar plot of total counts:
library(ggplot2)
ggplot(df.counts, aes(x=CHROM,y=value,fill=BaseSubstitution)) +
  geom_bar(position="stack", stat="identity")+
  xlab("Chromosome")+
  ylab("Number of substitutions")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(limits = c(paste0("chr", c(1:22,"X","Y","M")))) #這行可以指定x軸label的順序

# Stacked bar plot of proportions:
ggplot(df.counts, aes(x=CHROM,y=value,fill=BaseSubstitution)) +
  geom_bar(position="fill", stat="identity")+
  xlab("Chromosome")+
  ylab("Number of substitutions")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete(limits = c(paste0("chr", c(1:22,"X","Y","M"))))
#Most chromosomes have a similar distribution of counts, with the exception of chr21 (fewer G>C and T>A, more A>G) and chrM (only T>C, but there is only one variant there).
#Forum in moodle同學有分享如何讓x axis by order of chromosones

#extract the "allele frequency" information. Here is how to do it for a single instance
info <- snvs.keep[1,]$INFO 
info
#Split the string by "VariantAlleleFrequency=" using the strsplit function
splitinfo <- strsplit(info,"VariantAlleleFrequency=")
splitinfo[[1]][2] #The beginning of the second part contains the VAF
vaf <- strsplit(splitinfo[[1]][2],";")
vaf[[1]][1]
# Save all the extracted VAFs into a new column
snvs.keep$VAF <- sapply(snvs.keep$INFO,       #建立新column
                        function(x) strsplit(strsplit(x,"VariantAlleleFrequency=")[[1]][2],";")[[1]][1])
snvs.keep$VAF

snvs.keep$VAF <- as.numeric(snvs.keep$VAF) # Make sure this is seen as numeric

#Plot the histogram
hist(snvs.keep$VAF, breaks=50) #r內建plot

#ggplot
ggplot(snvs.keep, aes(VAF)) + 
  geom_histogram(aes(y=..density..),bins=100) +  #為何y=density???
  geom_density(col = "red") 

#Fitting a Gaussian mixture model
library(mixtools)
my_mix <- normalmixEM(snvs.keep$VAF, k = 4) #telling it to find two gaussians in the observations, k可改2, 3, 4...
my_mix
plot(my_mix, which=2) #which 都是設定2
#The graph shows the proposed mixture of two distributions that these allele frequencies may have been drawn from (which would roughly equate to the subclonal populations the respective mutations arose in).

#Fitting a Dirichlet process
library(dirichletprocess)
# Create a Dirichlet mixture of Gaussians:
dp <- DirichletProcessGaussian(snvs.keep$VAF)
# Fit the Dirichlet process object (10 iterations, ideally 1000s):
dp <- Fit(dp, its = 20) #1000 iteration要跑超久
# Plot the estimated density of the data:
plot(dp)

# Show outcome (number of clusters found):
dp
dp$numberClusters
table(dp$clusterLabels)
