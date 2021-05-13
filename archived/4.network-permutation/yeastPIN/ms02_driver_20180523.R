# 20180522 revise to deal with highly connected networks
# 20180523 add gaussian white noises

rm(list=ls())

# R -f ms02_driver_20180523.R --args bg.phys.csv bg.3.csv 100 0.025 9 > log.3.txt
#R -f file --args input_network_csv_filename output_csv_filename in_ncycles noise_weight debug

options(echo=TRUE) # if you want see commands in output file 
args <- commandArgs(trailingOnly = TRUE)
#print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandsArgs(trailingOnly=FALSE))
infile = args[1]
#infile = "bg.phys.csv"; 
outfile = args[2]
in_ncycles = as.integer(args[3])
in_noise_weight = as.numeric(args[4])
debug = as.integer(args[5])

#set.seed(20180523)
print(paste("infile=[",infile, "]\n"));

f_match_degree = function (inTb, noise.weight = 0.025 ) {
  inTb$DegreePercentile1 = degree$DegreePercentile[match(inTb$id1, degree$ORF)];
  inTb$DegreePercentile2 = degree$DegreePercentile[match(inTb$id2, degree$ORF)];
  inTb$DegreePercentileProduct = inTb$DegreePercentile1 * inTb$DegreePercentile2; 
  if (noise.weight > 0 ) {
    inTb$DegreePercentileProduct =
           inTb$DegreePercentileProduct + rnorm(n=length(inTb[,1]), sd=noise.weight); 
  }
  inTb = inTb[order(inTb$DegreePercentileProduct, decreasing = FALSE), ];
  return (inTb); 
}

f_redo_rest.tb = function ( inTb ) {
  redo.tb = inTb[ inTb$selfpairs==1 | inTb$tag_counts>1, ]
  rest.tb = inTb[ inTb$selfpairs==0 & inTb$tag_counts==1, ]
  return(list("redo.tb"=redo.tb, "rest.tb" = rest.tb))
}

# recursive_permutation_no_selfpairing_v0.01 rpns
#This function only perform permutation, but does NOT check for correctness of permutation. 
#Correctness of permutation should be check outside of this function!
recursive_permutation_no_selfpairing_v0.01 = function( inpairs,  ncycles=5, indebug=0, preserve_rate = 0.1, noise.weight=0.025 ) { 
    if(indebug>0) {
      cat(paste( "\n**", '(rpns) Start ncycles=', ncycles, ", preserve_rate=", preserve_rate, 
                 ", length(inpairs[,1])=",length(inpairs[,1]), "**\n\n" )); 
    }
    longids = c(as.character(inpairs[,1]), as.character(inpairs[,2]) )
    longids = sample(longids)
    len = length(inpairs[,1])
    
    newpairs0 = data.frame( cbind( longids[1:len], longids[(len+1): (2*len)]) )
    newpairs = t(apply(newpairs0, 1, sort)); #oder id1 and id2
    newpairs = data.frame(newpairs); 
    if(indebug > 2) {#check ids ordering results
      #print(paste("===before sort ===="),NULL); 
      #print(newpairs0); print(paste("===after sort ===="),NULL); 
      #print(newpairs);
      print(paste("===Before:Aftersort  ===="),NULL); 
      print(cbind( newpairs0, newpairs));
    }   
    names(newpairs) = c('id1', 'id2')
    newpairs$id1 = as.character( newpairs$id1)
    newpairs$id2 = as.character( newpairs$id2)    
    
    newpairs$tag =  paste(newpairs[,1], newpairs[,2], sep="_")
    counts = table( newpairs$tag )
    newpairs$tag_counts = counts[newpairs$tag]
    # if(indebug>8) {    counts;    }
    newpairs$selfpairs = ifelse( newpairs$id1 == newpairs$id2, 1, 0 )
    
    newpairs = f_match_degree(newpairs, noise.weight = noise.weight ); 
    
    redo.tb = newpairs[ newpairs$selfpairs==1 | newpairs$tag_counts>1, ]
    rest.tb = newpairs[ newpairs$selfpairs==0 & newpairs$tag_counts==1, ]
    
    if(indebug>8) {
      print(paste('(rpns)ncycles=', ncycles, "===redopairs===="),NULL); 
      print (redo.tb);
      print(paste('(rpns)ncycles=', ncycles, "===restpairs===="),NULL); 
      print(rest.tb);
      print(paste("================="),NULL)
    }
    
    if( length(redo.tb[,1])>=1 ) {
       if ( ( ncycles == 0) | (length(rest.tb[,1]) < 1) ) { 
        print(paste("ncycles reached zero, ncycles, OR, not enough data in rest.tb for randomization", ncycles) );
        if (indebug > 8 ) {
          write.csv(redo.tb, "tmp/redo_tb.csv");
          write.csv(rest.tb, "tmp/rest_tb.csv");
        }
        return( newpairs ); # no more randomization inside of this function
      } else  {
        degreeProduct.cutoff = quantile(rest.tb$DegreePercentileProduct, prob = (1 - preserve_rate)); 
        unchangedpairs = rest.tb[ rest.tb$DegreePercentileProduct >degreeProduct.cutoff, ]
        selectedpairs = rbind(redo.tb,rest.tb[ rest.tb$DegreePercentileProduct <=degreeProduct.cutoff , ] )   
        
        if (indebug > 0) {
          print(paste('(rpns) ncycles=', ncycles, ", degreeProduct.cutoff=",degreeProduct.cutoff, 
                      ", length(redo.tb[,1])=",length(redo.tb[,1]),
                      ", length(rest.tb[,1])=",length(rest.tb[,1]) ),NULL); 
        }
        ncycles = ncycles - 1; 
        return( rbind(unchangedpairs, 
                      recursive_permutation_no_selfpairing_v0.01(selectedpairs, ncycles,indebug = indebug, preserve_rate = preserve_rate )))#20180522, recursive trap bug
      }
    } else if (length(redo.tb[,1])==0) { #20180522
      if(indebug>0) {
        print(paste('(rpns) SUCCESS, ncycles=', ncycles, 
                    ", length(redo.tb[,1])=",length(redo.tb[,1]),
                    ", length(rest.tb[,1])=",length(rest.tb[,1]) ),NULL); 
      }      
      return (newpairs )
    } 

}#end of function


success_flag = 0; 
global_cycles = 5;

net = read.csv( infile, colClasses = c("character", "character") )
head(net)

longids = c( net[,1], net[,2] );
degree = sort( table( longids ), decreasing = TRUE ); 
degree = data.frame(degree)
names(degree) = c("ORF", "degree"); 
degree$ORF = as.character( degree$ORF); 
# Pick highly connected hub nodes, and prioritize their permutations
# Notes of caution: This means degree-degree profile are not random in permuted networks
degree$DegreePercentile = 0;
for ( i in 1:length(degree[,1])) {
  degree$DegreePercentile[i] = 1 - length(which(degree$degree> degree$degree[i]))/ length(degree[,1]);
}

#while (( success_flag == 0 & global_cycles > 0 ) ){
#  global_cycles = global_cycles - 1; 

  perNet0 = recursive_permutation_no_selfpairing_v0.01(net, ncycles = in_ncycles, indebug = debug, noise.weight = in_noise_weight, preserve_rate = 0.1); 
  #perNet0$DegreePercentile1 = degree$DegreePercentile[match(perNet0$id1, degree$ORF)];
  #perNet0$DegreePercentile2 = degree$DegreePercentile[match(perNet0$id2, degree$ORF)];
  #perNet0$DegreePercentileProduct = perNet0$DegreePercentile1 * perNet0$DegreePercentile2; 
  # summary( lm( perNet0$tag_counts ~ perNet0$DegreePercentileProduct) );
  perNet0 = f_match_degree( perNet0 );
  
  x = f_redo_rest.tb( perNet0 );
  redo.PerNet0 = x$redo.tb; 
  rest.PerNet0 = x$rest.tb;
  
  if ( length(redo.PerNet0[,1])==0 ) {
    write.csv(perNet0, outfile, quote=F, row.names=F)
  } else {
    system(paste("touch ", outfile, "_err")); 
  }
  
  
  #redo.PerNet0b = perNet0[ perNet0$selfpairs==1 | perNet0$tag_counts>1, ]
  #rest.PerNet0b = perNet0[ perNet0$selfpairs==0 & perNet0$tag_counts==1, ]

  #if (debug > 0 ) {
  # quantile(rest.PerNet0$DegreePercentileProduct, prob=c(0.1, 0.5, 0.80, 0.9)); 
  #}
  
  #keep half of the unqiue pairs with higher degrees unchanged
  #good.perNet =    rest.PerNet0[ rest.PerNet0$DegreePercentileProduct > median(rest.PerNet0$DegreePercentileProduct), ];
 
  # add half of the unique pairs back to the mixing pool
  #perNet0.redo = 
  #  rest.PerNet0[ rest.PerNet0$DegreePercentileProduct <= median(rest.PerNet0$DegreePercentileProduct), ];

  #redo.PerNet0b = rbind( redo.PerNet0, perNet0.redo  ); #redo this part
  #perNet1 = recursive_permutation_no_selfpairing_v0.01(redo.PerNet0b, ncycles = 10, indebug = 2, preserve_rate = 0.5); 
  #perNet1 = f_match_degree( perNet1 )
  #y = f_redo_rest.tb( perNet1 );
  #redo.PerNet1 = y$redo.tb; 
  #rest.PerNet1 = y$rest.tb;
  
  # check permutation to set success_flag
#}






