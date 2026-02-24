#ComputeLikelihoodRatios.R
#created by Tom Busey, June 2022
#create likelihood ratios for either the 


rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off()
library(tidyverse)
library(reshape2)
library(Hmisc)
library(metRology)



source("DBDA2E-utilities.R")

doBlackBox = TRUE #change to fit either the Busey et al data or the FBI/Black Box data

doingT = FALSE #if true, use the t distribution

if (doBlackBox==TRUE)
{
  if (doingT==FALSE)
  {load( file= "MCMCParameterSummaryBlackBox.rdata")}
  else{load( file= "MCMCParameterSummaryBlackBoxTDist.rdata")}
}else
{
  if (doingT==FALSE)
  {load(file= "MCMCParameterSummary.rdata")}
  else{load(file= "MCMCParameterSummaryTDist.rdata")}
}


muMedians = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"Median"]
muHDILow = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"HDIlow"]
muHDIHigh = OrdModelResults$OrdSummaryMat[grep("mu",rownames(OrdModelResults$OrdSummaryMat)),"HDIhigh"]

doBlackBox = OrdModelResults$doBlackBox #probably redundant

sortMuIdx = sort( muMedians,index.return=TRUE )$ix
rankMuIdx = rank( muMedians)
if (doBlackBox == FALSE)
{
  stimid = c(1:60)
  mated = if_else(stimid %% 2 == 0, "Mated", "Nonmated")
}else
{
  stimid = OrdModelResults$fingerprints_summary_all$stimid
  mated = if_else(OrdModelResults$fingerprints_summary_all$mated, "Mated", "Nonmated")
}

plottingTopXValue = 15
plottingBottomXValue = -8


topXValue = 24
bottomXValue = -14
numberOfBins =2000
AllLatentX = seq(bottomXValue , topXValue, length = numberOfBins)

binWidth = (topXValue-bottomXValue)/(numberOfBins-1)

MatedSum = rep(0, length(AllLatentX))
NonmatedSum = rep(0,length(AllLatentX))
graphFileType=c("pdf","png","eps","jpg")[1] 

openGraph(width=10.0,height=8.0)


#strategy- plot all 60 Normal distributions
maxY = 0

#thePlot <- ggplot()
lineAlpha = .5
lineWidth = 1

nonMatedColor = rgb(253/255,127/255,125/255,lineAlpha)
matedColor =rgb(125/255, 132/255, 253/255,lineAlpha)

bigNonMatedColor = rgb(1,0,0,1)
bigMatedColor =rgb(0,0,1,1)

bigLineWidth = 4

if (doBlackBox) #need to set the top of the scale manually. Technically the black box goes above 1.0, but these are relative likelihoods, not probabilities, and a probability density can be greater than 1.0
{maxY = .7}else
{maxY = .7}
if (doBlackBox)
{
  theTitle = "Ordered Probit Model Estimation (Black Box Data)"
}else
{
  theTitle = "Ordered Probit Model Estimation (Busey et al. Data)"
}

for (pairNumber in 1:length(OrdModelResults$NForEachPairNumber))
{
  
  
  medianLocation = OrdModelResults$OrdSummaryMat[paste0("mu[",pairNumber, "]"), "Median"]
  medianSigma = OrdModelResults$OrdSummaryMat[paste0("sigma[",pairNumber, "]"), "Median"]
  
  if (doingT == FALSE)
  {thisSubjectY = dnorm(AllLatentX, mean = medianLocation, sd = medianSigma, log = FALSE)}
  else
    {thisSubjectY = dt.scaled(AllLatentX, 5, mean = medianLocation, sd = medianSigma, ncp = 0, log = FALSE)}
  #dt.scaled(x, df, mean = 0, sd = 1, ncp, log = FALSE)
  
  
  thisSubjectY = thisSubjectY/sum(thisSubjectY)/binWidth
  maxY = max(maxY, max(thisSubjectY))
  
  if(mated[pairNumber]=="Mated")
  {
    theColor = matedColor
    
    MatedSum = MatedSum + thisSubjectY 
    
  }
  else
  {
    theColor = nonMatedColor
    
    NonmatedSum = NonmatedSum + thisSubjectY   
    
  }
  if (pairNumber == 1)
  {
    # plot.new()
    plot(AllLatentX, thisSubjectY, type = "l", col = theColor, xlim = c(plottingBottomXValue, plottingTopXValue),
         main = theTitle, font.main = 4, font.lab = 2,
         #log="y", #uncomment to see on log scale
         #ylim = c(1e-9, 1), 
         ylim = c(0, maxY),
         lwd = lineWidth, xlab="Latent Dimension (support for different or same sources)", ylab = "Probability Density (relative likelihood)")
    minor.tick(nx=5, ny=0, tick.ratio = 0.5, x.args = list(), y.args = list())
    
  }
  else
  {
    lines(AllLatentX, thisSubjectY, type = "l", col = theColor, xlim = c(plottingBottomXValue, plottingTopXValue),
          ylim = c(0, maxY), 
          lwd = lineWidth)
    
  }
}

#normalize both to create probability density functions
MatedSum = MatedSum/sum(MatedSum)/binWidth
NonmatedSum = NonmatedSum/sum(NonmatedSum)/binWidth

lines(AllLatentX, NonmatedSum, type = "l", col = bigNonMatedColor, lwd = bigLineWidth)
lines(AllLatentX, MatedSum, type = "l", col=bigMatedColor, lwd = bigLineWidth)

xlab="Latent Dimension (support for different or source sources)"

if (doBlackBox==TRUE)
{
  if (doingT==FALSE)
  {  saveGraph( file="CombinedPlotBlackBox", type=graphFileType)}
  else{saveGraph( file="CombinedPlotBlackBoxTDist", type=graphFileType)}
}else
{
  if (doingT==FALSE)
  {saveGraph( file="CombinedPlot", type=graphFileType)}
  else{saveGraph( file="CombinedPlotTDist", type=graphFileType)}
}


openGraph(width=7.0,height=7.0)
likelihood = MatedSum/NonmatedSum

if (doBlackBox)
{
  theTitle = "Likelihood Ratio (Black Box Data)"
}else
{
  theTitle = "Likelihood Ratio (Busey et al. Data)"
}

plot(AllLatentX, likelihood, type = "l",  xlim = c(plottingBottomXValue, plottingTopXValue),
     xlab="Latent Dimension (support for different or same sources)",
     lwd = 3,col = bigMatedColor,
     log="y", ylim = c(.01, 1e10), ylab = "Likelihood Ratio (log scale)",
     main = theTitle, font.main = 4, font.lab = 2, yaxt="n")

aty <- axTicks(2)
labels <- sapply(aty,function(i)
  as.expression(bquote( .(i)))
)
axis(2,at=aty,labels=labels)

grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
minor.tick(nx=5, ny=0, tick.ratio = 0.5, x.args = list(), y.args = list())


if (doBlackBox==TRUE)
{
  if (doingT==FALSE)
  {  saveGraph( file="LikelihoodRatioBlackBox", type=graphFileType)}
  else{saveGraph( file="LikelihoodRatioBlackBoxTDist", type=graphFileType)}
}else
{
  if (doingT==FALSE)
  {saveGraph( file="LikelihoodRatio", type=graphFileType)}
  else{saveGraph( file="LikelihoodRatioTDist", type=graphFileType)}
}




LikelihoodData = data.frame(MatedSum)
LikelihoodData$NonmatedSum = NonmatedSum
LikelihoodData$AllLatentX = AllLatentX
LikelihoodData$likelihood = likelihood

if (doBlackBox==TRUE)
{
  if (doingT==FALSE)
  {  save(LikelihoodData, file= "LikelihoodDataBlackBox.rdata")}
  else{save(LikelihoodData, file= "LikelihoodDataBlackBoxTDist.rdata")}
}else
{
  if (doingT==FALSE)
  {save(LikelihoodData, file= "LikelihoodData.rdata")}
  else{save(LikelihoodData, file= "LikelihoodDataTDist.rdata")}
}

#now write out data table if doing black box
if ((doBlackBox) & (doingT == FALSE))
{
  
  #now accumulate the counts and LRs
  Table1Colnames =  c('pairID', 'mu', 'sigma', 'LR', 'TradExcl', 'TradInc', 'TradID', 'TradNV')
  Table1DF <- data.frame(matrix(ncol = length(Table1Colnames), nrow = 0))
  colnames(Table1DF) = Table1Colnames
  
  for (pairNumber in 1:length(OrdModelResults$NForEachPairNumber))
  {
    
    
    medianLocation = OrdModelResults$OrdSummaryMat[paste0("mu[",pairNumber, "]"), "Median"]
    medianSigma = OrdModelResults$OrdSummaryMat[paste0("sigma[",pairNumber, "]"), "Median"]
    
    EstimatedLikelihood = approx(LikelihoodData$AllLatentX, LikelihoodData$likelihood, medianLocation)
    
    
    #now update dataframe
    Table1DF[pairNumber, 'pairID'] = OrdModelResults$fingerprints_summary_all$Pair_ID[pairNumber]
    Table1DF[pairNumber, 'mu'] = medianLocation
    Table1DF[pairNumber, 'sigma'] = medianSigma
    Table1DF[pairNumber, 'LR'] = EstimatedLikelihood$y
    
    Table1DF[pairNumber, 'TradExcl'] = OrdModelResults$fingerprints_summary_all$Exclusion[pairNumber]
    Table1DF[pairNumber, 'TradInc'] = OrdModelResults$fingerprints_summary_all$Inconclusive[pairNumber]
    Table1DF[pairNumber, 'TradID'] = OrdModelResults$fingerprints_summary_all$Individualization[pairNumber]
    Table1DF[pairNumber, 'TradNV'] = OrdModelResults$fingerprints_summary_all$NV[pairNumber]
    
  }
  Table1DF = Table1DF[order(Table1DF$LR),]
  write.table(Table1DF, "BlackBoxDataForTable2.csv", sep = ",", col.names = TRUE, row.names = FALSE)
  
  decimatedRows = seq(1,nrow(Table1DF), 10)
  
  Table1DFDecimated = Table1DF %>% slice(decimatedRows)
  
  write.table(Table1DFDecimated, "BlackBoxDataForTable2Decimated.csv", sep = ",", col.names = TRUE, row.names = FALSE)
}