 library(limma);
 library(GEOquery);

 gse = getGEO("GSE3821");
 gpl = getGEO("GPL90")
 Meta( gse );
 names(GSMList(gse));

 probesets = Table( GPLList(gse)[[1]])$ID

 data.matrix = do.call("cbind", lapply(GSMList(gse),function(x) {
	tab = Table(x)
	mymatch = match( probesets, tab$ID_REF)
	return( tab$VALUE[mymatch] )
}))  #no log2 transformation

 require(Biobase);
 rownames( data.matrix) = probesets;
 colnames( data.matrix) = names( GSMList(gse) );

 rawtotals = apply(data.matrix, 2, sum); 
 CONST = 1.5E6
 scales = CONST / rawtotals
 names(scales) = colnames( data.matrix );

 data.matrix2 = data.matrix;
 for( j in 1:length(scales) ){
  data.matrix2[,j] = data.matrix[,j] * scales[j]; 
 }

ORFs = Table( GPLList(gse)[[1]])$ORF

 #### pearson correlation 21:31->32 ?? done within 1 minute! 
 pearsoncor.m = cor( t(data.matrix2) );
 rownames( pearsoncor.m ) = ORFs;
 colnames( pearsoncor.m ) = ORFs;   

 #clean up the momory
 rm(gse);
 rm( data.matrix );

#Now, map pearson data to ppi data
 dip.tb = read.table( "large.dip.pairs.tab", sep="\t");
 names( dip.tb ) = c("orf1", "orf2", "n");

 dip.tb$pearson = NA;

 for( i in 1:length(dip.tb[,1]) ) {
  dip.tb$pearson[i] = pearsoncor.m[ match( dip.tb[i, 1], ORFs ),  match( dip.tb[i, 2], ORFs ) ];
 }

### write.table( dip.tb[,c(1,2,4)], "dip.pearson.glucose.pulse.csv", quote=F,sep="\t", row.name=F);

 #### parse mitochondrial gene subsets
 #localization.tb = read.table( "localization/OSheaLocalization.WeissmanAbundance.tab", header=T,sep="\t", fill=T);
 #summary ( localization.tb[,c("mitochondrion")] )
 #mito.orfs = localization.tb$orfid[ localization.tb$mitochondrion ==T ]
 #tmp = mito.orfs[ ! is.na(mito.orfs) ]
 #mito.orfs = as.character( tmp );

 file = "/home/hqin//projects/BEI08.hainan/key.data/orfs.mito.plus.glucose.path.tab";
 mito.tb = read.table(file, header=F )
 mito.orfs = as.character( mito.tb[,1] )

 # ORFs[ match( mito.orfs, ORFs ) ];
 matched.pos = match( mito.orfs, ORFs )
 matched.pos2 = matched.pos[ ! is.na(matched.pos) ]
 mito.orfs2 = mito.orfs[ is.na(matched.pos)==F  ];

 mito.matrix = data.matrix2[ matched.pos2,  ]
 rownames( mito.matrix ) = mito.orfs2;  #I need to double check this. 

 length(mito.orfs); # [1] 500
 length(mito.orfs2); # [1] 497

 ### pearson correlation
 mito.pearson = cor( t(mito.matrix) ); 
 hist(mito.pearson); #strong correlation due to temporal correlation

 ### now partial correlation for mito genes
 require(ggm);
 require(GeneNet);
 mito.shrinkage = ggm.estimate.pcor( t(mito.matrix) ); 
  #Estimating optimal shrinkage intensity lambda (correlation matrix): 0.3955
 hist(mito.shrinkage); #most pearson correlation disappear

 ### map mitopcor to ppi
 match( mito.orfs2,  as.character( dip.tb[, 1] )  ) #test 

 dip.tb$mitopcor = NA;
 #  i = 498
 for( i in 1:length(dip.tb[,1]) ) {
  dip.tb$mitopcor[i] = mito.shrinkage[ match( as.character(dip.tb[i,1]), mito.orfs2 ),  
+ match( as.character(dip.tb[i, 2]), mito.orfs2 )  ];
 }

 summary( dip.tb );
 hist( dip.tb$mitopcor );
 hist( dip.tb$pearson );

 write.table( dip.tb[,c(1,2,5)], "dip.mitopcor.glucose.pulse.csv", quote=F,sep="\t", row.name=F);
 
 save.image("RData.112407")
