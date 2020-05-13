#2012 Nov 1

rm(list=ls())
#library(Biobase)
require(GEOquery)

#http://www.ncbi.nlm.nih.gov/geo/browse/
#GSE20444  GSE27222  GSE4295   http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2529
gset <- getGEO("GSE3821", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset) #This is the expression matrix

dictionary = gset@featureData@data[, c('ID', 'ORF')]  #This is a lookup table for ID and ORF 

#########
# A simple approach to create an expression matrix with ORFs as row names
# This approach takes only one probe for each ORFs, which is often true for cDNA arrays
ORFs = unique(as.character(dictionary$ORF))
ex2 = ex[match(ORFs, dictionary$ORF), ]   
rownames(ex2) = ORFs
head(ex2) #Now, expression matrix is named by ORFs

##########
#Another approach is to calculate the average sigals for all the probes in the same ORFs

multipleORFs = NA;
ex3 = ex2 #This is just a template
# orf = 'YLR331C'
for (orf in ORFs) {
  myrows = as.character( dictionary$ID[dictionary$ORF==orf] )
  if (length(myrows) > 1) {
    print (orf)
    multipleORFs = c(multipleORFs, orf)
    ex[myrows, ] = apply(ex[myrows,], 2, mean) 
  }else {
    ex3[orf, ] = ex[myrows[1], ]
  }
}
multipleORFs = multipleORFs[-1]

######
#normalizaion is needed. 


#########
# now, have a look at the signals
hist(ex3[,1], br=100)
ex4 = log2(ex3)
hist(ex4[,3])
ex4[ex4<0] = NaN

#############
#calculate coefficient of variation
myVar = apply( ex4, 1, FUN=var)
myStddev = sqrt(myVar)
myMean = apply( ex4, 1, FUN=mean)
myCV = myStddev / myMean
myarray= data.frame(cbind( myStddev, myMean, myCV))
myarray$ORF = ORFs
myarray = myarray[, c(4, 1:3)]

write.csv(myarray, "GSE3821_log2CV.csv", row.names=F)
test = read.csv("GSE3821_log2CV.csv", colClasses = c('character', NA, NA, NA))
str(test)
hist(test$myCV, br=100)
hist(test$myStddev, br=100)
hist(test$myMean, br=100)
