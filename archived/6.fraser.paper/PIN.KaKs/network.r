#2014 April 9, Hong Qin hqin@spelman.edu

#todo current simulation does not consider essen-essen interaction

#2015 Fall

summarize_mean_from_files = function(infiles, inputdir){
  debug = 0; 
  #calcualte means of networking aging data in a vector of files
  myfiles = infiles; inputdir = inputdir;
  if( debug) { print (myfiles)}; 
  outtb = data.frame(myfiles)
  outtb$mean = NA;
  for( i in 1:length(myfiles)) {
    currentFile = paste( inputdir, '/', myfiles[i], sep='');
    tb = read.csv(currentFile)
    outtb$mean[i] = mean(tb[,1])
  }
  outtb; 
}

multiple_network_failure = function(INpopSize, INlambda1, INlambda2=INlambda1/10, INthreshold=5, probability, INpairs, INessenLookupTb ) {
  popSize = INpopSize; lambda1=INlambda1; lambda2=INlambda2;
  threshold=INthreshold; p= probability; pairs=INpairs; essenLookupTb=INessenLookupTb; 
  
  popAges = numeric(popSize)
  j=1; count = 0; 
  while ((j <= popSize) && ( count < popSize*30)) {
    count = count + 1;      
    print(paste("count=",count))
    currentNetworkAge = single_network_failure_v2(lambda1, lambda2, degreeThreshold, p, pairs, essenLookupTb)
    if (currentNetworkAge > 0) {
      popAges[j] = currentNetworkAge      
      j = j+1
    } 
  }# end of j while-loop, population loop
  
 return(popAges)
}#END


############# 
# To use lambda1 for all edges, choose threshold = 0
single_network_failure_v2 = function(lambda1, lambda2=lambda1/10, threshold=5, p, pairs, essenLookupTb ) {
  # single network failure simulation, 20151013Tue
  # lambda1: First exponential constant failure rate for edges with degree >= threshold
  # lambda2: Second exponential constant failure rate for edges with degree < threshold
  # threshold: degree threshold for lambda1 and lambda2
  # pairs: network in pairwide format, using numeric NOs 20151013
  # essenLookupTb: lookup table for essential and nonessential genes, numeric values 
  ## for debug:   lambda1 = 1/50; lambda2= lambda1/10; threshold=4; p=0.8
  
  inpairs = pairs[,1:2] #bookkeeping  
  names(inpairs) = c('No1','No2')
  
  #get connectivities per node
  degreeTb = data.frame( table(c(inpairs$No1, inpairs$No2)))
  names(degreeTb) = c("No", "degree")
  degreeTb$moduleAge = NA;
  
  for( i in 1:length(degreeTb[,1])){
    if ( essenLookupTb[ degreeTb$No[i] ] != 0) { #essential node
      lambda = ifelse( degreeTb$degree[i] >= threshold, lambda1, lambda2)
      age = rexp( degreeTb$degree[i], rate=lambda ) #exponential age
      if(degreeTb$degree[i] >= threshold){
        active = runif(degreeTb$degree[i])  #uniform interaction stochasticity
        active = ifelse( active<=p, 1, NA  ) #pick active interactions
        if( sum(active, na.rm=T) > 0 ){ #there should be at least 1 active intxn
          age = age * active # only active interactions for modular age estimation
          degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
        } else {# when no active intxn is available 
          degreeTb$moduleAge[i] = 0; #this module is born dead.
        }
      } else { # for degree < threshold, no stochasticity is applied. 
        degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
      }
    } else {# non-essential node
      degreeTb$moduleAge[i] = NA 
    }
  }
  
  summary(degreeTb)
  currentNetworkAge = min(degreeTb$moduleAge, na.rm=T)
}


#20140408 old ms02_singlerun() did not check id1-id2 versus id2-id1. 
# So, I wrote v2 and wrapp the old function to v2 function call. 
#permute.pairs.wo.selfpairs = function( inpairs,  ncycles=10, debug=1 ) {
ms02_singlerun = function( inpairs,  ncycles=10, indebug=0 ) { # Renamed, 2014 Feb 12
  return( ms02_singlerun_v2( inpairs,  ncycles=ncycles, indebug=indebug ))
}

ms02_singlerun_v2 = function( inpairs,  ncycles=10, indebug=0 ) { 
  if (ncycles >= 1 ) {
    if(indebug>0) {
      print(paste('ncycles=', ncycles))
    }
    longids = c(as.character(inpairs[,1]), as.character(inpairs[,2]) )
    longids = sample(longids)
    len = length(inpairs[,1])
    newpairs2 = data.frame( cbind( longids[1:len], longids[(len+1): (2*len)]) )
    newpairs2 = t(apply(newpairs2, 1, sort))
    newpairs2 = data.frame(newpairs2)
    names(newpairs2) = c('id1', 'id2')
    newpairs2$id1 = as.character( newpairs2$id1)
    newpairs2$id2 = as.character( newpairs2$id2)    
    
    newpairs2$tag =  paste(newpairs2[,1], newpairs2[,2], sep="_")
    counts = table( newpairs2$tag )
    newpairs2$tag_counts = counts[newpairs2$tag]
    
    newpairs2$selfpairs = ifelse( newpairs2$id1 == newpairs2$id2, 1, 0 )
    
    redo.tb = newpairs2[ newpairs2$selfpairs==1 | newpairs2$tag_counts>1, ]
    rest.tb = newpairs2[ newpairs2$selfpairs==0 & newpairs2$tag_counts==1, ]
    if(indebug>0) {
      print(paste("===redopairs===="),NULL);      print(redo.tb);
      #print(paste("===restpairs===="),NULL);      print(rest.tb);
      print(paste("================="),NULL)
    }
    if( length(redo.tb[,1])>=1 ) {
      if ( ncycles == 0) { 
        #return (c(NA,NA, NA) );
        print(paste("ncycles reached zero, ncycles"),ncycles)
        print(paste("Abort!"),NULL)
        stop; 
      } else {
        ncycles = ncycles - 1
        splitPos = round( length(redo.tb[,1]) * sqrt(ncycles) ) + 5
        splitPos = min( splitPos, (length(rest.tb[,1])-1 ) )
        selectedpairs = rbind(redo.tb,  rest.tb[1: splitPos, ] )   #20140408, potential bug. always take initial section
        unchangedpairs = rest.tb[ (splitPos + 1): length(rest.tb[,1]), ] #20140408, potential bug. 
        return( rbind(unchangedpairs, ms02_singlerun_v2(selectedpairs, ncycles)))  #2014 Feb 12
      }
    } else {  
      return (newpairs2 )
    }
  } else {
    return( c(NA,NA,NA )) 
  }
}#end of ms02 v2

#old version
single_network_failure = function(lambda, p, pairs, runningORFs) {
  # single network failure simulation
  # lambda: exponential constant failure rate for edges
  # pairs: network in pairwide format
  # runningORFs: GooddEssentialORFsPPI  
  
  inpairs = pairs[,1:2] #bookkeeping  
  names(inpairs) = c('id1','id2')
  
  #stochasticity into pairs   
  inpairs$active = runif(length(inpairs[,1]))  #uniform
  # tmp = pairs$active > 1-p
  # table(tmp) / length(tmp)  ; #double-check, very good. 
  
  inpairs$age = rexp( length(inpairs[,1]), rate=lambda )  #exponential ages for pairs
  inpairs$age = ifelse(inpairs$active > (1-p), inpairs$age, NA ) #if not active, intxn is excluded. 
  #pairs$age = ifelse(pairs$active > (1-p), pairs$age, 0 )  # in what situations, can non-ative intxn be treat as 0-age?
  
  ModuleTb = data.frame(runningORFs) #buffer for module ages    
  #loop every essential genes to identify the module age
  for (i in 1:length(runningORFs)) {
    myORF = runningORFs[i]
    pos1 = grep(myORF, inpairs$id1)
    pos2 = grep(myORF, inpairs$id2)  #id1,2 to ORF1,2 is a really bad choice. 
    if( length( c(pos1,pos2))>=1 ) {
      ModuleTb$age.m[i] = max( inpairs$age[c(pos1,pos2)], na.rm=T )   #maximal intxn age -> module age
    } else {
      ModuleTb$age.m[i] = NA; 
    }
  }
  #head(ModuleTb); 
  summary(ModuleTb)
  ModuleTb$age.m[ ModuleTb$age.m== -Inf] = 0; #dead births occur when links are not active
  currentNetworkAge = min(ModuleTb$age.m)
}


