####################################################################################
#
# OrderedProbitModel.R
# runs the ordinal probit model on response data, and also plots some graphs

#based on 
# John K. Kruschke, January-August 2018. 
# Sept 18, 2018: Revised by John K. Kruschke; graphical output revised.
# Sept 28, 2018: Revised by John K. Kruschke; graphical output revised more, and
#                argument caseIDcolName added.
#
### Description:
#
# This R script defines a function, ordinalAnalysis(), that analyzes
# groups of ordinal data with an ordered probit model.
#
# The function is used for an example in the article, "Analyzing ordinal data
# with metric models: What could possibly go wrong?" by Torrin M. Liddell and
# John K. Kruschke. Files can be found at https://osf.io/53ce9/
#
### Usage:
#
# source("OrderedProbitModel.R") # source this file to load the function
#
# ordinalResults = ordinalAnalysis( dataFileName , yColNames ,
#                                            caseIDcolName=NULL ,
#                                            compareCases=NULL , 
#                                            graphFileType="pdf" ,
#                                            hierarchSD=FALSE )
#
# For a detailed example of using this function, see the script
# OrderedProbitModel-Example.R
#
### Arguments:
#
# dataFileName: The name of the data file. Assumed to be comma separated values.
#               See below for more details.
#
# yColNames: A vector of column names of the frequencies of ordinal levels.
#
# caseIDcolName: Column name of row identifiers or labels. Optional, defaults to
#                NULL.
#
# compareCases: A list of case pairs to compare. Optional, defaults to NULL.
#
# graphFileType: The graphics format of saved graphs. Defaults to "pdf".
#
# hierarchSD: Whether or not to put hierarchical structure on the group standard
#             deviations. Logical TRUE/FALSE, defaults to FALSE.
#
### Output (Value):
#
# A list of two components:
#
# OrdSummaryMat: A matrix, with rows being the parameters in the ordinal model
# and columns being summary statistics of the posterior distribution.
#
# OrdMcmcMat: A matrix, with rows being steps in the MCMC chain and columns
# being the chain number and the parameters in the ordinal model.
#
# N.B.: The function saves the output to the working directory and also produces
# numerous graphs that are saved to the working directory.
#
### Important Dependencies:
#
# This script assumes a particular data structure: The dependent (predicted)
# variable is an ordinal value from a single item. The independent (predictor)
# variable is group. For example, the dependent variable could be a 1 to 5 star
# rating with the independent variable being different movies. (The model makes
# no assumptions about repeated measures, that is, the model does not represent
# any identities of responders that might be the same across groups.) For each
# group, the data are specified as a count of each level of rating scale. The
# data file must be strucured with named columns, one column for the group
# identifier, and K columns for the counts of each level of the 1 through K
# ordinal levels. See OrderedProbitModel-Example.R for more info.
#
# The file DBDA2E-utilities.R is needed for various functions such as
# gammaShRaFromModeSD(), diagMCMC(), saveGraph(), plotPost(), HDIofMCMC(), etc.
# It also installs and loads the R package runjags, so you must be connected to
# the internet if you do not already have runjags installed. The file
# DBDA2E-utilities.R should accompany this file at the article's OSF repository,
# or you can find it in the DBDA2Eprograms.zip file from Step 5 of
# https://sites.google.com/site/doingbayesiandataanalysis/software-installation
# Get the DBDA2Eprograms.zip file, extract it, and copy the file
# DBDA2E-utilities.R into the same folder (working directory) as this file.
#
### Example:
#
# For a detailed example of using this function, see the script
# AnalyzeValidationAndBlackBoxDataWithOrderedProbitModel.R
#
####################################################################################

source("DBDA2E-utilities.R")

ordinalAnalysis = function( datFrm , yColNames , 
                            caseIDcolName=NULL ,
                            compareCases=NULL , 
                            graphFileType=c("pdf","png","eps","jpg")[2] ,
                            hierarchSD=FALSE,
                            doBlackBox, doingT, tDF) {
  fileNameRoot = paste0("OrderedProbitModel-")
  #-----------------------------------------------------------------------------------
  # Set up data for JAGS
  #datFrm = read.csv( dataFileName , header=TRUE )
  # yColNames must be vector of names IN ORDER from 1 to nYlevels
  y = as.matrix(datFrm[,yColNames]) # y is matrix
  isThree = as.matrix(datFrm[,"isThree"])
  conclusionType = as.matrix(datFrm[,"conclusionType"])
  
  #  doingT = TRUE
  # tDF = 5
  
  z = rowSums(y)
  x = 1:nrow(datFrm)
  Ncases = nrow(datFrm)
  #doingT = doingT
  #nu = tDF
  stimID = as.matrix(datFrm[,"stimid"])
  stimID = as.vector(stimID)
  #numStims = nrow(unique(stimID))
  numStims = length(unique(stimID))
  # Threshold 1 and nYlevels-1 are fixed; other thresholds are estimated.
  # This allows all parameters to be interpretable on the response scale.
  nYlevels = length(yColNames) 
  
  threethresh = rep(NA,nYlevels-1)
  threethresh[1] = 1.5
  threethresh[2] = 1.5 #not fit for this model
  threethresh[3] =  nYlevels-1 + 0.5 #not fit for this model
  threethresh[nYlevels-1] = nYlevels-1 + 0.5
  
  threshETrad = rep(NA,nYlevels-1)
  # threshETrad[1] = 1 + 0.5
  #  threshETrad[nYlevels-1] = nYlevels-1 + 0.5
  
  threshESOS = rep(NA,nYlevels-1)
  # threshESOS[1] = 1 + 0.5
  #threshESOS[nYlevels-1] = nYlevels-1 + 0.5
  
  # Specify the data in a list, for later shipment to JAGS:
  gammaShRa = unlist( gammaShRaFromModeSD( mode=3.0 , sd=3.0 ) )
  # For the ordered probit model:
  if( hierarchSD ) {
    OrdDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      threethresh = threethresh ,
      threshETrad = threshETrad ,
      threshESOS = threshESOS ,
      conclusionType = as.double(conclusionType),
      isThree = as.double(isThree),
      x = x ,
      z = z ,
      stimID = stimID, 
      numStims = numStims,
      Ncases = Ncases,
      gammaShRa = gammaShRa
    ) 
  } else {
    OrdDataList = list(
      y = y ,
      nYlevels = nYlevels ,
      threethresh = threethresh ,
      threshETrad = threshETrad ,
      threshESOS = threshESOS ,
      conclusionType = as.double(conclusionType),
      isThree = as.double(isThree),
      x = x ,
      z = z ,
      stimID = stimID,
      numStims = numStims,
      Ncases = Ncases
      
    )
  }
  
  #-----------------------------------------------------------------------------------
  # THE *ORDERED PROBIT* MODEL:
  #pT with df = 5, df=1 cauchy       zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] , 1/zsigma^2 , nu )
  #    nu ~ dexp(1/30.0)
  
  #throw out weak and see what happens
  #change prior for mu
  
  modelString = paste0("
  model {
    for ( i in 1:Ncases ) {
      y[i, ] ~ dmulti( pr[i,1:nYlevels] , z[i] )
                       ", 
                       ifelse(doingT, 
                              "
                              pr[i,1] <- (isThree[i]) * (pt( threethresh[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) ) + 
        (1-isThree[i]) * (conclusionType[i]) * (pt( threshESOS[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) ) + 
        (1-isThree[i]) * (1-conclusionType[i]) * (pt( threshETrad[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) )",
                              
                              "pr[i,1] <- (isThree[i]) * (pnorm( threethresh[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ) ) + 
        (1-isThree[i]) * (conclusionType[i]) * (pnorm( threshESOS[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ) ) + 
        (1-isThree[i]) * (1-conclusionType[i]) * (pnorm( threshETrad[1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ))"
                       ),
                       "
     
      #pr[i,1] <- 1 
 
      for ( k in 2:(nYlevels-1) ) {", 
                       ifelse(doingT, 
                              "
                              pr[i,k] <- (isThree[i]) * (max( 0 ,  pt( threethresh[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 )
                             - pt( threethresh[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) ) ) +
                      (1-isThree[i]) * (conclusionType[i]) * (max( 0 ,  pt( threshESOS[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 )
                          - pt( threshESOS[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) ) ) +  
                       (1-isThree[i]) * (1-conclusionType[i]) * (max( 0 ,  pt( threshETrad[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 )
                             - pt( threshETrad[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ) ) )
        ",
                              "
         pr[i,k] <- (isThree[i]) * (max( 0 ,  pnorm( threethresh[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 )
                             - pnorm( threethresh[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ) ) ) +
                      (1-isThree[i]) * (conclusionType[i]) * (max( 0 ,  pnorm( threshESOS[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 )
                          - pnorm( threshESOS[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ) ) ) +  
                       (1-isThree[i]) * (1-conclusionType[i]) * (max( 0 ,  pnorm( threshETrad[ k ] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 )
                             - pnorm( threshETrad[k-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ) ) )"
                       ),
                       "
              #      pr[i,k] <- 1        
      }
      ", 
                       ifelse(doingT, 
                              "
      pr[i,nYlevels] <- (isThree[i]) * (
       1 - pt( threethresh[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 )) + 
       (1-isThree[i])*(conclusionType[i]) *  (
       1 - pt( threshESOS[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 )) + 
       (1-isThree[i])*(1-conclusionType[i]) *  (
       1 - pt( threshETrad[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2, 5 ))
                              ",
                              "
                                pr[i,nYlevels] <- (isThree[i]) * (
       1 - pnorm( threethresh[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 )) + 
       (1-isThree[i])*(conclusionType[i]) *  (
       1 - pnorm( threshESOS[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 )) + 
       (1-isThree[i])*(1-conclusionType[i]) *  (
       1 - pnorm( threshETrad[nYlevels-1] , mu[stimID[i]] , 1/sigma[stimID[i]]^2 ))"
                       ),
                       "
      #pr[i,nYlevels] <- 1
    }
    for ( j in 1:numStims ) { 
      mu[j] ~ dnorm( (1+nYlevels)/2 , 1/(nYlevels)^2 )#play around and see if it is sensitive 1/(10*nYlevels)^2)- spread outs out black box data, reduces LRs
      sigma[j] ~ dgamma( sigmaSh , sigmaRa )
    }
    sigmaSh <- 1 + sigmaMode * sigmaRa
    sigmaRa <- ( ( sigmaMode + sqrt( sigmaMode^2 + 4*sigmaSD^2 ) ) 
                  / ( 2*sigmaSD^2 ) ) ",
                       ifelse( hierarchSD ,
                               "sigmaMode ~ dgamma( gammaShRa[1] , gammaShRa[2] ) 
       sigmaSD ~ dgamma( gammaShRa[1] , gammaShRa[2] ) " ,
                               "sigmaMode <- 3.0
       sigmaSD <- 3.0" ) , " # open quote for next line
    #for ( k in 2:(nYlevels-2) ) {  # 1 and nYlevels-1 are fixed, not stochastic
    #  threethresh[k] ~ dnorm( k+0.5 , 1/2^2 )
    #} 
    for ( k in 1:(nYlevels-1) ) {  #expanded scales estimated all four thresholds
      threshETrad[k] ~ dnorm( k+0.5 , 1/2^2 )
    } 
    for ( k in 1:(nYlevels-1) ) {  #expanded scales estimated all four thresholds
      threshESOS[k] ~ dnorm( k+0.5 , 1/2^2 )
    }
  }") # close quote for modelString paste
  # Write out modelString to a text file
  writeLines( modelString , con=paste0(fileNameRoot,"-OrdModel.txt") )
  #-----------------------------------------------------------------------------------
  # Run the ordered-probit model:
  if ( hierarchSD ) {
    if (doBlackBox)
    {
      parameters = c( "mu" , "sigma" , "threethresh" ,  "sigmaMode" , "sigmaSD" )
    }
    else
    {
      parameters = c( "mu" , "sigma" , "threethresh" , "threshETrad" , "threshESOS" , "sigmaMode" , "sigmaSD" )
      
    }
  } else {
    if (doBlackBox)
    {
      parameters = c( "mu" , "sigma" ,  "threethresh" ) 
      
    }else
    {
      parameters = c( "mu" , "sigma" ,  "threethresh" , "threshETrad" , "threshESOS" ) 
    }
  }
  adaptSteps = 500
  burnInSteps = 1000
  #  numSavedSteps = 24000
  numSavedSteps = 48000
  thinSteps = 5
  nChains = 4
  runJagsOut <- run.jags( method="parallel" ,
                          model=paste0(fileNameRoot,"-OrdModel.txt") , 
                          monitor=parameters , 
                          data=OrdDataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  OrdCodaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(fileNameRoot) ) {
    if (doBlackBox)
    {
      if (doingT)
      {
        save( OrdCodaSamples , file=paste0(fileNameRoot,"-OrdModel-McmcBlackBoxTDist.Rdata") )}
      else{
        save( OrdCodaSamples , file=paste0(fileNameRoot,"-OrdModel-McmcBlackBox.Rdata") )
      }
    }else
    {
      if (doingT)
      {
        save( OrdCodaSamples , file=paste0(fileNameRoot,"-OrdModel-McmcTDist.Rdata") )}
      else{
        save( OrdCodaSamples , file=paste0(fileNameRoot,"-OrdModel-Mcmc.Rdata") )}
    }
  }
  
  
  #-----------------------------------------------------------------------------------
  # For output, consider case names:
  if ( !is.null(caseIDcolName) ) { 
    if (doBlackBox)
    {
      caseNames = as.character( datFrm[,"Pair_ID"]  )
      
    }
    else
    {
      #caseNames = as.character( datFrm[,caseIDcolName]  )
      caseNames = paste( "Pair" , 1:Ncases )
    }
    
  } else {
    caseNames = paste( "Case" , 1:Ncases )
  }
  #-----------------------------------------------------------------------------------
  # Display diagnostics of chains:
  if (doBlackBox)
  {
    PathToSubfolder = "BlackBoxMCMCDiagnosticPlots/"
  }else
  {
    PathToSubfolder = "ValidationModelMCMCDiagnosticPlots/"
  }
  OrdParameterNames = varnames(OrdCodaSamples) 
  if (doingT)  {suffix = "TDist"}else {suffix = ""}
  if ( FALSE ) {
    for ( parName in OrdParameterNames ) {
      diagMCMC( codaObject=OrdCodaSamples , parName=parName ,
                saveName=paste0(PathToSubfolder, fileNameRoot,"-OrdModel", suffix) , #
                saveType=graphFileType )
      graphics.off()#close graphs to avoid overflow
    }
  }
  
  #-----------------------------------------------------------------------------------
  # Display and output posterior information:
  OrdMcmcMat = as.matrix(OrdCodaSamples,chains=TRUE)
  OrdChainLength = NROW( OrdMcmcMat )
  
  # Create numerical summary object for returning from this function:
  postStats = c("Mean","Median","Mode","HDImass","HDIlow","HDIhigh","ESS") # +"psrf" below
  OrdSummaryMat = matrix( 0 , nrow=length(OrdParameterNames) , ncol=length(postStats)+1 )
  colnames(OrdSummaryMat) = c(postStats,"psrf")
  rownames(OrdSummaryMat) = OrdParameterNames
  for ( parName in OrdParameterNames ) {
    OrdSummaryMat[parName,postStats] = summarizePost(OrdMcmcMat[,parName],credMass=0.95)[postStats]
    OrdSummaryMat[parName,"psrf"] = coda::gelman.diag( OrdCodaSamples[,parName] )$psrf[1,1]
  }
  #show(OrdSummaryMat)
  
  
  
  return(list( OrdSummaryMat=OrdSummaryMat, 
               OrdMcmcMat=OrdMcmcMat, 
               fileNameRoot=fileNameRoot, 
               hierarchSD = hierarchSD,
               graphFileType = graphFileType,
               OrdDataList = OrdDataList,
               doingT = doingT,
               caseNames) )
  
}

plotThresholds = function (threshCols, OrdMcmcMat, fileNameRoot, condName, hierarchSD, graphFileType, yColNames){
  threshMean = rowMeans( OrdMcmcMat[,threshCols] )
  xLim = range(OrdMcmcMat[,threshCols])
  xLim = c( xLim[1]-0.5 , xLim[2]+0.5 )
  xTickPos = apply( OrdMcmcMat[,threshCols] , 2 , median )
  nPtToPlot = 1000
  plotIdx = floor(seq(1,nrow(OrdMcmcMat),length=nPtToPlot))
  openGraph(width=7.0,height=5.0)
  layout( matrix( c( rep(1,length(xTickPos)-2) , (1:(length(xTickPos))+1) ) ,
                  nrow=2 , byrow=TRUE ) )
  par( mar=c(3.5,3.5,3,1) , mgp=c(2.25,0.7,0) )
  # Plot thresh mean x thresh:
  plot( OrdMcmcMat[plotIdx, threshCols[1]] , threshMean[plotIdx] ,
        xlab="Latent Scale" , ylab="Mean Threshold" , cex.lab=1.5 ,
        xlim=xLim , xaxt="n" , #xaxp=xTickPos ,
        main=bquote("Thresholds ("*theta[k]*")") , cex.main=1.75 ,
        col="skyblue" )
  abline(v=mean(OrdMcmcMat[plotIdx,threshCols[1]]),lty="dashed",col="skyblue")
  for ( i in 2:length(threshCols) ) {
    points( OrdMcmcMat[plotIdx,threshCols[i]] , threshMean[plotIdx] , col="skyblue" )
    abline(v=mean(OrdMcmcMat[plotIdx,threshCols[i]]),lty="dashed",col="skyblue")
  }
  axis( 1 , at=xTickPos , labels=round(xTickPos,2) , cex.axis=1.5 )
  for ( k in 1:(length(threshCols)+1) ) {
    text( x=mean(c(1.0,xTickPos,length(xTickPos)+1)[c(k,k+1)]) , 
          y=mean(threshMean) , labels=paste0("'",k,"'") , cex=1.5 )
  }
  # Plot marginals of individual thresh's:
  for ( tIdx in (1:(length(xTickPos))) ) {
    plotPost( OrdMcmcMat[,threshCols[tIdx]] , main=bquote(theta[.(tIdx)]) ,
              xlab="Latent Scale" , cex.main=1.75 )
  }
  if ( !is.null(fileNameRoot) ) {
  }
}

plotDataWithMeans = function (z, y, stimID, 
                              numStims, whichScale, thisConclusionType, 
                              OrdMcmcMat, datFrm,
                              fingerprints_summary_allUnfiltered,
                              yColNames, caseNames,fileNameRoot, graphFileType,
                              fileDesc, NForEachPairNumber,
                              NForEachPairNumberUnfiltered,
                              titleForGraph,
                              doBlackBox, doingT, tDF)
{  
  #-----------------------------------------------
  titleCaseNameMaxLength = 10 # maximum characters to display from label
  
  actuallyPlotMeans = TRUE
  
  
  # Plot data with ORDERED PROBIT posterior predictions and latent mu,sigma annotation:
  nCell = numStims
  #nCol=ceiling(sqrt(nCell))
  #nRow=ceiling(nCell/nCol)
  nCol = 10
  nRow = 6
  #nCol = min(nCol,10)
  #nRow = min(nRow,6)
  inchPerPanel = 1.7
  
  numGraphsNeeded = floor((numStims-1)/((nCol*nRow)))+1
  graphsPerPanel = (nCol*nRow)
  
  sortMuIdx = 1:length(grep("mu",colnames(OrdMcmcMat)))
  if(doBlackBox)#randomize the order to get a selection on the first page
  {
    sortMuIdx = sample(sortMuIdx)
  }
  for (thisGraph in c(1:numGraphsNeeded))
  {
    
    openGraph(width=nCol*inchPerPanel,height=nRow*inchPerPanel)
    layout(matrix(1:(nCol*nRow),nrow=nRow,ncol=nCol,byrow=TRUE))
    par( mar=c(3,3,4,1) , mgp=c(1.8,0.7,0) )
    
    
    #sortMuIdx = sort( apply( OrdMcmcMat[,grep("mu",colnames(OrdMcmcMat))] , 2 , median ) ,
    #                  index.return=TRUE )$ix
    
    #just get sequence for current panel
    currentMusToGraph = sortMuIdx[seq((thisGraph-1)*graphsPerPanel+1,min((thisGraph)*graphsPerPanel, numStims))]
    
    for ( caseIdx in currentMusToGraph ) {
      if (whichScale == 3)
      {
        dataForThisCase = datFrm %>% 
          filter(stimid == caseIdx) %>% 
          filter(condition == whichScale) 
        
        dataForThisCaseUnfiltered = fingerprints_summary_allUnfiltered %>% 
          filter(stimid == caseIdx) %>% 
          filter(condition == whichScale) 
      }
      if (whichScale == 5)
      {
        dataForThisCase = datFrm %>% 
          filter(stimid == caseIdx) %>% 
          filter(condition == whichScale) %>% 
          filter(conclusionType == thisConclusionType)
        
        dataForThisCaseUnfiltered = fingerprints_summary_allUnfiltered %>% 
          filter(stimid == caseIdx) %>% 
          filter(condition == whichScale) %>% 
          filter(conclusionType == thisConclusionType)
      }
      # if (whichScale == 3)
      # {thisY = data.frame(colSums(dataForThisCase,2))}
      
      thisY = c(sum(dataForThisCase$`1`), sum(dataForThisCase$`2`), sum(dataForThisCase$`3`), sum(dataForThisCase$`4`), sum(dataForThisCase$`5`) )
      thisYUnfilered = c(sum(dataForThisCaseUnfiltered$`1`), sum(dataForThisCaseUnfiltered$`2`), sum(dataForThisCaseUnfiltered$`3`), sum(dataForThisCaseUnfiltered$`4`), sum(dataForThisCaseUnfiltered$`5`) )
      
      n = sum(thisY)
      nUnfiltered = sum(thisYUnfilered)
      NForEachPairNumber[caseIdx] = NForEachPairNumber[caseIdx] + n
      NForEachPairNumberUnfiltered[caseIdx] = NForEachPairNumberUnfiltered[caseIdx] + n
      
      if (actuallyPlotMeans)
      {
        
        plot( 1:length(yColNames) , thisY , 
              xlab="Rating" , ylab="Frequency" , cex.lab=0.75 ,
              xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(thisY)) ,
              type="h" , lend=1 , lwd=10 , col="pink" )
        medMu = median(OrdMcmcMat[,paste0("mu[",caseIdx,"]")])
        medSigma = median(OrdMcmcMat[,paste0("sigma[",caseIdx,"]")])
        
        if (doBlackBox)
        {      
          titleCaseName = paste0( caseIdx , ":" , 
                                  caseNames[caseIdx] )
        }else{titleCaseName = paste0( caseIdx  )}
        
        title( main=bquote( atop( .(titleCaseName)
                                  *", N="*.(n)*""
                                  
                                  #  *", NTotal="*.(nUnfiltered) *"  "
                                  ,
                                  mu*"="*.(round(medMu,2))
                                  *", "*sigma*"="*.(round(medSigma,2))
        ) ) , cex.main=1.0 )
        # superimpose posterior predictions:
        predProb = matrix(0,nrow=nrow(OrdMcmcMat),ncol=length(yColNames))
        for ( stepIdx in 1:nrow(OrdMcmcMat) ) {
          if (whichScale == 3)
          {threshVec = OrdMcmcMat[stepIdx,grep("threethresh",colnames(OrdMcmcMat))]}
          if (whichScale == 5)
          {
            if (thisConclusionType == 0)
            {
              threshVec = OrdMcmcMat[stepIdx,grep("threshETrad",colnames(OrdMcmcMat))]
            }
            if (thisConclusionType == 1)
            {
              threshVec = OrdMcmcMat[stepIdx,grep("threshESOS",colnames(OrdMcmcMat))]
            }
          }
          thisMu = OrdMcmcMat[stepIdx,paste0("mu[",caseIdx,"]")]
          thisSigma = OrdMcmcMat[stepIdx,paste0("sigma[",caseIdx,"]")]
          predProb[stepIdx,] = ( pnorm( (c(threshVec,Inf)-thisMu)/thisSigma )
                                 - pnorm( (c(-Inf,threshVec)-thisMu)/thisSigma ) )
        }
        points( 1:length(yColNames) , 
                apply(predProb,2,median) * n  , 
                cex=1 , lwd=2 , col="skyblue" )
        for ( respIdx in 1:length(yColNames) ) {
          lines( rep(respIdx,2) , HDIofMCMC( predProb[,respIdx] ) * n , 
                 col="skyblue" , lwd=3 )
        }
      }
    }
    if (doingT==TRUE)  {suffix = "TDist"}else {suffix = ""}
    
    if ( !is.null(fileNameRoot) ) {
      saveGraph( file=paste0(fileNameRoot,fileDesc,"-OrdModel-PostPred", titleForGraph,suffix, "-", thisGraph ), type=graphFileType)
    }
  }#looping across graphs
  return(list(NForEachPairNumber=NForEachPairNumber, NForEachPairNumberUnfiltered=NForEachPairNumberUnfiltered))
}

plotParametersAfterMCMC = function(OrdModelResultsList, datFrm, fingerprints_summary_allUnfiltered, yColNames)
{
  
  OrdMcmcMat = OrdModelResultsList$OrdMcmcMat
  OrdSummaryMat = OrdModelResultsList$OrdSummaryMat
  fileNameRoot = OrdModelResultsList$fileNameRoot
  hierarchSD = OrdModelResultsList$hierarchSD
  doingT = OrdModelResultsList$doingT
  graphFileType = OrdModelResultsList$graphFileType
  OrdDataList = OrdModelResultsList$OrdDataList
  numStims = OrdDataList$numStims
  y = OrdDataList$y
  stimID = OrdDataList$stimID
  # caseNames = OrdModelResultsList$caseNames
  caseNames = OrdModelResultsList$fingerprints_summary_all$Pair_ID
  
  z = OrdDataList$z
  
  #-----------------------------------------------
  # Plot thresholds of ordered-probit model:
  threethreshCols = grep("threethresh",colnames(OrdMcmcMat),value=TRUE)
  
  plotThresholds(threethreshCols, OrdMcmcMat, fileNameRoot, "Three Scale", hierarchSD, graphFileType, yColNames)
  
  if (doBlackBox == FALSE)
  {
    threshETradCols = grep("threshETrad",colnames(OrdMcmcMat),value=TRUE)
    
    plotThresholds(threshETradCols, OrdMcmcMat, fileNameRoot, "Expanded Traditional Scale", hierarchSD, graphFileType, yColNames)
    
    threshESOSCols = grep("threshESOS",colnames(OrdMcmcMat),value=TRUE)
    
    plotThresholds(threshESOSCols, OrdMcmcMat, fileNameRoot, "Expanded SOS Scale", hierarchSD, graphFileType, yColNames)
  }
  #keep track of N for each pair number
  NForEachPairNumber = rep(0, numStims)
  NForEachPairNumberUnfiltered = rep(0, numStims)
  
  whichScale = 3
  thisConclusionType = 0
  if (doBlackBox == FALSE)
  {
    titleForGraph = "Trad"
  }else{
    titleForGraph = "BlackBoxTrad"
  }
  NForEachPairNumberList =  plotDataWithMeans(z, y, stimID, numStims, whichScale, 
                                              thisConclusionType, OrdMcmcMat, datFrm,
                                              fingerprints_summary_allUnfiltered,
                                              yColNames, caseNames, fileNameRoot, graphFileType, "ThreeScale",
                                              NForEachPairNumber, NForEachPairNumberUnfiltered,
                                              titleForGraph, doBlackBox, doingT)
  NForEachPairNumber = NForEachPairNumberList$NForEachPairNumber
  NForEachPairNumberUnfiltered = NForEachPairNumberList$NForEachPairNumberUnfiltered
  
  if (doBlackBox == FALSE)
  {
    whichScale = 5
    thisConclusionType = 0
    titleForGraph = "ExpandTrad"
    
    NForEachPairNumberList =  plotDataWithMeans(z, y, stimID, numStims, whichScale, 
                                                thisConclusionType, OrdMcmcMat, datFrm,
                                                fingerprints_summary_allUnfiltered,
                                                yColNames, caseNames, fileNameRoot, graphFileType, "FiveScaleETrad",
                                                NForEachPairNumber,NForEachPairNumberUnfiltered,
                                                titleForGraph, doBlackBox, doingT)
    
    NForEachPairNumber = NForEachPairNumberList$NForEachPairNumber
    NForEachPairNumberUnfiltered = NForEachPairNumberList$NForEachPairNumberUnfiltered
    
    whichScale = 5
    thisConclusionType = 1
    titleForGraph = "SOS"
    
    NForEachPairNumberList= plotDataWithMeans(z, y, stimID, numStims, whichScale, 
                                              thisConclusionType, OrdMcmcMat, datFrm,
                                              fingerprints_summary_allUnfiltered,
                                              yColNames, caseNames, fileNameRoot, graphFileType, "FiveScaleESOS",
                                              NForEachPairNumber,NForEachPairNumberUnfiltered,
                                              titleForGraph, doBlackBox, doingT)
    
    NForEachPairNumber = NForEachPairNumberList$NForEachPairNumber
    NForEachPairNumberUnfiltered = NForEachPairNumberList$NForEachPairNumberUnfiltered
  }
  #-----------------------------------------------
  if ( hierarchSD ) {
    # Plot sigmaMode, sigmaSD of ordered-probit model:
    openGraph(width=7,height=3.5)
    layout(matrix(1:2,nrow=1,byrow=TRUE))
    par( mar=c(4,2,4,1) , mgp=c(2.25,0.7,0) )
    plotPost( OrdMcmcMat[,"sigmaMode"] , xlab=bquote(omega[sigma[i]]) ,
              main="Ordered-Probit Modal Sigma" )
    plotPost( OrdMcmcMat[,"sigmaSD"] , xlab=bquote(sigma[sigma[i]]) ,
              main="Ordered-Probit SD of Sigma[i]" )
    if (doingT==TRUE)  {suffix = "TDist"}else {suffix = ""}
    
    if ( !is.null(fileNameRoot) ) {
      saveGraph( file=paste0(fileNameRoot,"-OrdModel-sigmaMode", titleForGraph, suffix), type=graphFileType)
    }
    
  }
  
  return( list(NForEachPairNumber = NForEachPairNumber, NForEachPairNumberUnfiltered) )
} # end of function definition
