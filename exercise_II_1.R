##Exercise 1
#install.packages("ape")
# First, load the library for phylogenetic analysis:
library(ape)

# Build Newick string:
text.string<-
  "((EGFR,ERBB2),KRAS);"
# Read tree from string:
vert.tree1<-read.tree(text=text.string)
# Tree output
vert.tree1
str(vert.tree1)
#Visualize the first tree
plot(vert.tree1,no.margin=TRUE,edge.width=1)

#We have a cancer phylogeny where a subclone driven by KRAS emerged earlier, and the rest of the clonal population then split into another 2 subclones, driven by either EGFR or ERBB2.


# Build Newick string:
text.string<-
  "(((((((EGFR, ERBB2),KRAS),(TP53,(CDKN2A,KIF13))),(PIK3CA,SMAD2)),APC),BRAF),HRAS);"
vert.tree2<-read.tree(text=text.string)
str(vert.tree2)
plot(vert.tree2,no.margin=FALSE,edge.width=2)

#Highlight the TP53-containing subtree in red
clcolr <- rep("darkgrey", nrow(vert.tree2$edge))
clcolr[c(10:14)] <- "red"
plot(vert.tree2, edge.color=clcolr, no.margin=FALSE,edge.width=2)

#Extract the subtree as indicated (重打一次tree)
# Build Newick string (the third tree):
text.string<-
  "(TP53,(KIF13,CDKN2A));"
# Read tree from string:
vert.tree3<-read.tree(text=text.string)
plot(vert.tree3,no.margin=TRUE,edge.width=1)

#We can, for example, check to see if our trees are equal. Trees with rotated nodes are equal. Re-rooted trees are not; however they are the same if unrooted. 
## Check if all three trees are equal:
all.equal(vert.tree1,vert.tree2,vert.tree2)
## Check if the first and third tree are equal
all.equal(vert.tree1,vert.tree3)

text.string<-"((KIF13,CDKN2A),TP53);"
vert.tree4<-read.tree(text=text.string)
plot(vert.tree4,no.margin=TRUE,edge.width=1)

## Are the 3rd and 4th tree equal?
all.equal(vert.tree3,vert.tree4) #True, 圖找不一樣沒關係，genes之間的相關性是一樣的就好


##Exercise 2
#install IRanges from Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("IRanges")
#BiocManager::install("limma")
##install devtools if you don't have it already
#install.packages("devtools")
#library(devtools)
#install_github("genome/bmm")
#install_github("genome/sciClone")

# Reading the variant allele frequencies for 2 tumours
library(sciClone)
load("vafs1.RData")
load("vafs2.RData")
load("cn1.RData") #copy number file
load("cn2.RData")

# Reading the variant allele frequencies for 2 tumours: 
# Explore the data tables:
head(vafs1) #dim: 5000 x 5
head(vafs2) #dim: 3981 x 5
head(cn1) #dim: 32 x 4
head(cn2) #dim: 88x 4

#In the VAF files, the chromosome location and exact position of each mutation are listed in each line. 
#RC_ref represents the read counts for the reference allele, RC_alt the read counts for the alternate allele. 
#VAF stands for the variant allele frequency. The copy number files contain the chromosome location, start and end of the segment with altered copy number. 
#The CN column indicates the total estimated copy number (not an integer because of the calculation from a noisy signal).

#1D clustering of mutations by copy number states on just one sample
sc = sciClone(vafs=vafs1,
              copyNumberCalls=cn1,
              sampleNames="sample1")
#Create output
writeClusterTable(sc, "clusters1.txt")
sc.plot1d(sc)

#1D clustering of mutations by copy number states for the other sample:
sc2 = sciClone(vafs=vafs2,
               copyNumberCalls=cn2,
               sampleNames="sample2")
#Create output
writeClusterTable(sc2, "clusters2.txt")
sc.plot1d(sc2)

#You will notice a difference in cluster distribution between samples S1 and S2, both in terms of changes in allele frequency as well as mutation distribution based on copy number state. 
#This suggests changes in clonal composition between the two samples (e.g. before/after therapy).

#2D clustering using two samples
sc = sciClone(vafs=list(vafs1,vafs2),
              copyNumberCalls=list(cn1,cn2),
              sampleNames=c("S1","S2"))
writeClusterTable(sc, "clusters2.txt")
sc.plot1d(sc) #1-dimensional
sc.plot2d(sc) #2-dimensional

#Inspect the output. 
#Notice that the clonal reconstruction using both samples together enables a straightforward identification of subclones 1 and 2 that had a low frequency in sample S1 and expanded to ~80% in sample S2. 
#The 2D comparison plot highlights this much more clearly than the VAF vs copy number plot.(左上角的圓形和三角形)


##Exercise 3: Constructing and visualizing clonal cancer evolution
library(devtools)
#install_github("chrisamiller/fishplot")
#install_github('hdng/clonevol')
#install.packages('gridBase')
#install.packages('gridExtra')
#install.packages('ggplot2')
#install.packages('igraph')
#install.packages('packcircles')
#install_github('hdng/trees')
#install_github("chrisamiller/fishplot")

# Load library:
library("clonevol")
# Load data frame with information on mutation clusters for the AML patient:
load("variants.aml.RData")
# Explore the data structure:
head(variants.aml)
#Information is provided about the clusters identified, the genes associated with each cluster (driver and passanger events, see is.driver), along with variant allele frequencies (vaf) and cancer cell fractions (ccf) of the primary (P) and relapse (R) samples.

#Next, infer the phylogenetic tree models that comply with these data
# Check the parameters for the function used to construct clonal models:
?infer.clonal.models #載不下來？？？！！！

# prepare sample grouping information:
sample.names <- c("P","R")
sample.groups <- sample.names
names(sample.groups) <- c("P.vaf","R.vaf")

# Defining some colours for the subclones
clone.colors <- c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')

# Infer clonality for the given data (check documentation for parameter definitions):
y = infer.clonal.models(variants = variants.aml,
                        cluster.col.name = 'cluster',
                        vaf.col.names = c("P.vaf","R.vaf"),
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        clone.colors = clone.colors,
                        min.cluster.vaf = 0.01,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)

#Build a consensus tree based on driver events:
consTree <- transfer.events.to.consensus.trees(y,
                                               variants.aml[variants.aml$is.driver,],
                                               cluster.col.name = 'cluster',
                                               event.col.name = 'gene')
# Create a list of fish objects 
library(fishplot)
f = generateFishplotInputs(results=consTree)
fishes = createFishPlotObjects(f)

# Plot with fishplot
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
           vlines=seq(1, length(sample.names)), vlab=sample.names, pad.left=0.5)
}
#ClonEvol can plot both node-based tree (each clone is a node), or branch-based tree (each branch represents the evolution of a clone from its parental clone, and each node represents a point where the clone is established/founded). 
#Before we can draw the latter tree, we need to prepare it.
branchTree <- convert.consensus.tree.clone.to.branch(consTree, branch.scale = 'log2')
plot.all.trees.clone.as.branch(branchTree, branch.width = 0.5,
                               node.size = 0.5, node.label.size = 0.3)






