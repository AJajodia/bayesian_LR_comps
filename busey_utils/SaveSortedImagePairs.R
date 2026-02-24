
#this creates summary images for each image from the Busey et al (2022) dataset
#also creates a summary dataframe for the table for the paper
#occasionally I get an error from 
#  Stats\temp.png': Permission denied @ error/blob.c/OpenBlob/2924 
# which I think is due to onedrive, but YMMV.

#note that this only applies to the Busey et al (2022) dataset. We don't have access to the black box images.

#this also creates the data table for Table 1. Table 2 is created only for the black box data by running ComputeLikelihoodRatios.R

rm(list=ls())

library(magick)
library(tidyverse)
source("DBDA2E-utilities.R")

ActuallySaveImages = TRUE #

load(file= "LikelihoodData.rdata")

PathToImages = "stimsCorrected/"
PathToSortedImages = "ImagesCombinedAndSorted/"

load("MCMCParameterSummary.rdata")

tradyColNames = c("Ex","N/A","Inc","N/A","ID")
EtradyColNames = c("Ex","SDS","Inc","SCS","ID")
ESOSyColNames = c("ESSDS","SDS","Inc","SCS","ESSCS")


muMedians = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"Median"]
muHDILow = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"HDIlow"]
muHDIHigh = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"HDIhigh"]

sortMuIdx = sort( muMedians,index.return=TRUE )$ix
rankMuIdx = rank( muMedians)
stimid = c(1:60)
mated = if_else(stimid %% 2 == 0, "Mated", "Nonmated")

plotOneGraphWithPosteriors = function (whichScale, thisConclusionType, datFrm, 
                                       fingerprints_summary_allUnfiltered,
                                       caseIdx, OrdMcmcMat)
{
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
  thisYUnfiltered = c(sum(dataForThisCaseUnfiltered$`1`), sum(dataForThisCaseUnfiltered$`2`), sum(dataForThisCaseUnfiltered$`3`), sum(dataForThisCaseUnfiltered$`4`), sum(dataForThisCaseUnfiltered$`5`) )
  
  n = sum(thisY)
  nUnfiltered = sum(thisYUnfiltered)
  
  
  if (whichScale == 3)
  {
    ConditionTitle = "Traditional Scale"
    yColNames = tradyColNames
  }
  else
  {
    if (thisConclusionType == 0)
    {
      ConditionTitle = "Expanded Traditional"
      yColNames = EtradyColNames
    }
    else
    {
      ConditionTitle = "Expanded Strength of Support"
      yColNames = ESOSyColNames
    }
  }
  
  medMu = median(OrdMcmcMat[,paste0("mu[",caseIdx,"]")])
  medSigma = median(OrdMcmcMat[,paste0("sigma[",caseIdx,"]")])
  EstimatedLikelihood = approx(LikelihoodData$AllLatentX, LikelihoodData$likelihood, medMu)
  
  if(ActuallySaveImages)
  {
    openGraph(width=7.0,height=5.0)
    
    plot( 1:length(yColNames) , thisY , 
          xlab="" , ylab="Frequency" , cex.lab=0.75 ,
          xlim=c(0.5,length(yColNames)+0.5) , ylim=c(0,1.2*max(thisY)) ,
          type="h" , lend=1 , lwd=30 , col="pink", xaxt = "n", cex.lab =1.2,cex.axis = 1.2,  )
    axis(1, at=1:5, labels = yColNames, cex.axis = 1.5)
    
    titleCaseName = paste0( ConditionTitle , ":" , 
                            caseIdx )
    
    
    title( main=bquote( atop( .(titleCaseName)
                              *", N="*.(n)*""
                              
                              *", NTotal="*.(nUnfiltered) *"  "
                              ,
                              mu*"="*.(round(medMu,2))
                              *", "*sigma*"="*.(round(medSigma,2))
                              *", "*LR*"="*.(round(EstimatedLikelihood$y,1))
    ) ) , cex.main=1.5,cex.lab = 10 )
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
            cex=1.5 , lwd=2 , col="skyblue" )
    for ( respIdx in 1:length(yColNames) ) {
      lines( rep(respIdx,2) , HDIofMCMC( predProb[,respIdx] ) * n , 
             col="skyblue" , lwd=2 )
    }
    saveGraph( file="temp", type="png")
  }
  
  returnData = list(mu = medMu, sigma = medSigma, likelihood = EstimatedLikelihood,
                    n = n,nTotal = nUnfiltered,thisY = thisY,thisYUnfiltered = thisYUnfiltered,
                    yColNames = yColNames)
  
  return(returnData)
}

Table1Colnames =  c('pairID', 'mu', 'sigma', 'LR', 'TradExcl', 'TradInc', 'TradID', 'TradNV',
                    'ETradExcl','ETradSDS', 'ETradInc', 'ETradSCS', 'ETradID', 'ETradNV',
                    'SOSESSDS', 'SOSSDS','SOSInc', 'SOSSCS', 'SOSESSCS', 'SOSNV')
Table1DF <- data.frame(matrix(ncol = length(Table1Colnames), nrow = 0))
colnames(Table1DF) = Table1Colnames

allStims = unique(OrdModelResults$fingerprints_summary_all$stimid)
for (pairNumber in allStims)
{
  graphics.off()
  
  leftImage = image_read(paste0(PathToImages, "Pair", sprintf("%02d", pairNumber),"Left.png"))
  rightImage = image_read(paste0(PathToImages, "Pair", sprintf("%02d", pairNumber),"Right.png"))
  
  fullImage = c(leftImage, rightImage)
  
  fullImage  = image_append(fullImage)
  thisRank = rankMuIdx[pairNumber];
  medianLocation = OrdModelResults$OrdSummaryMat[paste0("mu[",pairNumber, "]"), "Median"]
  medianSigma = OrdModelResults$OrdSummaryMat[paste0("sigma[",pairNumber, "]"), "Median"]
  
  EstimatedLikelihood = approx(LikelihoodData$AllLatentX, LikelihoodData$likelihood, medianLocation)
  
  
  fullImage = image_annotate(fullImage, paste0("Pair ", sprintf("%02d", pairNumber),  "  Rank ", 
                                               thisRank, "  Mu = ", sprintf("%3.2f", medianLocation),
                                               " [",sprintf("%3.2f", muHDILow[pairNumber])," to ",sprintf("%3.2f", muHDIHigh[pairNumber]),"]", "  Sigma = ", sprintf("%3.2f",medianSigma), "  ",
                                               "LR=",round(EstimatedLikelihood$y,1), " ",mated[pairNumber]), color = "red", size = 30)
  
  whichScale = 3
  thisConclusionType = 0
  ThreeReturnData = plotOneGraphWithPosteriors (whichScale, thisConclusionType, 
                                                OrdModelResults$fingerprints_summary_all, 
                                                OrdModelResults$fingerprints_summary_allUnfiltered,
                                                pairNumber, OrdModelResults$OrdMcmcMat)
  currentGraphThree = image_read("temp.png")
  
  whichScale = 5
  thisConclusionType = 0
  ETradReturnData = plotOneGraphWithPosteriors (whichScale, thisConclusionType, 
                                                OrdModelResults$fingerprints_summary_all,
                                                OrdModelResults$fingerprints_summary_allUnfiltered,
                                                pairNumber, OrdModelResults$OrdMcmcMat)
  currentGraphFiveTrad = image_read("temp.png")
  
  whichScale = 5
  thisConclusionType = 1
  SOSReturnData = plotOneGraphWithPosteriors (whichScale, thisConclusionType, 
                                              OrdModelResults$fingerprints_summary_all, 
                                              OrdModelResults$fingerprints_summary_allUnfiltered,
                                              pairNumber, OrdModelResults$OrdMcmcMat)
  currentGraphFiveSOS = image_read("temp.png")
  
  allgraphs = c(currentGraphThree, currentGraphFiveTrad, currentGraphFiveSOS)
  allgraphs = image_append(allgraphs, stack = TRUE)
  fullImage = c(fullImage, allgraphs)
  fullImage  = image_append(fullImage)
  
  
  if (ActuallySaveImages)
  {
    image_write(fullImage, 
                paste0(PathToSortedImages,"Rank", sprintf("%02d", thisRank),
                       "Pair", sprintf("%02d", pairNumber), mated[pairNumber],".png"))
  }
  
  #now update dataframe
  Table1DF[pairNumber, 'pairID'] = pairNumber
  Table1DF[pairNumber, 'mu'] = ThreeReturnData$mu
  Table1DF[pairNumber, 'sigma'] = ThreeReturnData$sigma
  Table1DF[pairNumber, 'LR'] = ThreeReturnData$likelihood$y
  
  Table1DF[pairNumber, 'TradExcl'] = ThreeReturnData$thisY[1]
  Table1DF[pairNumber, 'TradInc'] = ThreeReturnData$thisY[3]
  Table1DF[pairNumber, 'TradID'] = ThreeReturnData$thisY[5]
  Table1DF[pairNumber, 'TradNV'] = ThreeReturnData$nTotal-ThreeReturnData$n
  
  Table1DF[pairNumber, 'ETradExcl'] = ETradReturnData$thisY[1]
  Table1DF[pairNumber, 'ETradSDS'] = ETradReturnData$thisY[2]
  Table1DF[pairNumber, 'ETradInc'] = ETradReturnData$thisY[3]
  Table1DF[pairNumber, 'ETradSCS'] = ETradReturnData$thisY[4]
  Table1DF[pairNumber, 'ETradID'] = ETradReturnData$thisY[5]
  Table1DF[pairNumber, 'ETradNV'] = ETradReturnData$nTotal-ETradReturnData$n
  
  
  Table1DF[pairNumber, 'SOSESSDS'] = SOSReturnData$thisY[1]
  Table1DF[pairNumber, 'SOSSDS'] = SOSReturnData$thisY[2]
  Table1DF[pairNumber, 'SOSInc'] = SOSReturnData$thisY[3]
  Table1DF[pairNumber, 'SOSSCS'] = SOSReturnData$thisY[4]
  Table1DF[pairNumber, 'SOSESSCS'] = SOSReturnData$thisY[5]
  Table1DF[pairNumber, 'SOSNV'] = SOSReturnData$nTotal-SOSReturnData$n
  
  
  
}

#now save the dataframe
write.table(Table1DF, "ValidationDataForTable1.csv", sep = ",", col.names = TRUE, row.names = FALSE)

