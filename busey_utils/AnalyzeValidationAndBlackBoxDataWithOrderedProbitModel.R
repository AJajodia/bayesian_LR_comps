
#AnalyzeValidationAndBlackBoxDataWithOrderedProbitModel.R
# Thomas A. Busey, June 2022

#based on 
# OrderedProbitModel-Example.R
# John K. Kruschke, January-August 2018.
# Part of Doing Bayesian Data Analysis
# https://sites.google.com/site/doingbayesiandataanalysis/




# To use this script, you must have in R's working directory the following four
# files:
# (1) OrderedProbitModel.R
# (2) DBDA2E-utilities.R (included, but part of the site above)
# (3) results2.csv, which is the data from Busey et al (2022)
# (4) FBINoblisBlackBoxTestResponses.csv, which is the data downloaded from the FBI/Noblis site linked in the main manuscript

#make sure to set the working directory, and you may need to install jags from this link:
#https://sourceforge.net/projects/mcmc-jags/


#
# See important comments at the top of OrderedProbitModel.R for more information
# about the function.
#
#-----------------------------------------------------------------------------------
# Optional housekeeping:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#-----------------------------------------------------------------------------------
source("OrderedProbitModel.R") #need to source this first; commented out so we can debug

#read in data

library(tidyverse)
library(brms)

doBlackBox = TRUE #set to false to do the Busey et al (2022) data, TRUE to do the FBI/Noblis black box data
MinimumNumberOfExaminersInBB = 16

doingT = FALSE #if true, use a t distribution instead of normal distribution
tDF = 5 #not actually used; hard coded in orderedprobitmodel.r

#To turn off shrinkage, see this line near line 201: hierarchSD = TRUE

loadFromFile = FALSE #after you run the MCMC once for each dataset, set this to true to speed things up

set.seed(47408)#mainly just so that the plots the same order when doing individual cells 

if (loadFromFile == FALSE)
{
  
  if (doBlackBox==FALSE)
  {
    
    for (doFiltering in c(TRUE, FALSE))
    {
      fingerprints <- read.csv("results2.csv", header = TRUE, sep = ";") %>%
        # create a new column indicating if the prints were mated or nonmated
        # even stimid numbers are mated pairs, odd numbers are nonmated pairs
        mutate(mated = if_else(stimid %% 2 == 0, 1, 0))
      
      fingerprintsUnfilteredByValue = fingerprints
      #delete no-value trials if necessary
      if (doFiltering == TRUE)
        fingerprints = fingerprints %>%
        filter(valuedecision>0)
      
      #compute N with and without no value trials
      
      
      #make results one-based
      fingerprints$decision=fingerprints$decision+1
      
      # summary table
      fingerprints_summary <- fingerprints %>%
        #filter(conclusionType == 1) %>%
        filter(condition == 3) %>%
        group_by(stimid, conclusionType, condition, decision) %>%  #group_by(condition, conclusionType, mated, decision) %>%
        summarize(n = n()) %>% 
        pivot_wider(names_from = decision, values_from = n) %>%
        select(condition, conclusionType, `1`, `3`, `5`)
      #select( `0`, `1`, `2`, `3`, `4`)
      #select(condition, conclusionType, mated, `0`, `1`, `2`, `3`, `4`)
      
      fingerprints_summary$`2` = rep.int(0, nrow(fingerprints_summary))
      fingerprints_summary$`4` = rep.int(0, nrow(fingerprints_summary))
      fingerprints_summary$isThree = rep.int(1, nrow(fingerprints_summary))
      fingerprints_summary$condition = rep.int(3, nrow(fingerprints_summary))
      
      fingerprints_summary_trad = fingerprints_summary %>%
        select(stimid, condition, conclusionType, isThree,  `1`, `2`, `3`, `4`, `5`)
      
      
      fingerprints_summary_expanded_sos <- fingerprints %>%
        filter(conclusionType == 1) %>%
        filter(condition == 5) %>%
        group_by(stimid, conclusionType, condition, decision) %>%  #group_by(condition, conclusionType, mated, decision) %>%
        summarize(n = n()) %>% 
        pivot_wider(names_from = decision, values_from = n) %>%
        select(condition, conclusionType, `1`, `2`, `3`, `4`, `5`)
      fingerprints_summary_expanded_sos$isThree = rep.int(0, nrow(fingerprints_summary_expanded_sos))
      fingerprints_summary_expanded_sos = fingerprints_summary_expanded_sos %>%
        select(condition, conclusionType, isThree, `1`, `2`, `3`, `4`, `5`)
      
      
      
      
      fingerprints_summary_expanded_trad <- fingerprints %>%
        filter(conclusionType == 0) %>%
        filter(condition == 5) %>%
        group_by(stimid, conclusionType, condition, decision) %>%  #group_by(condition, conclusionType, mated, decision) %>%
        summarize(n = n()) %>% 
        pivot_wider(names_from = decision, values_from = n) %>%
        select(condition, conclusionType, `1`, `2`, `3`, `4`, `5`)
      fingerprints_summary_expanded_trad$isThree = rep.int(0, nrow(fingerprints_summary_expanded_trad))
      fingerprints_summary_expanded_trad = fingerprints_summary_expanded_trad %>%
        select(condition, conclusionType, isThree, `1`, `2`, `3`, `4`, `5`)
      
      if (doFiltering == TRUE)
      {
        fingerprints_summary_all=rbind(fingerprints_summary_trad,fingerprints_summary_expanded_sos)
        fingerprints_summary_all=rbind(fingerprints_summary_all,fingerprints_summary_expanded_trad)
        
        fingerprints_summary_all = fingerprints_summary_all %>% tidyr::replace_na(list( '1' = 0, '2' = 0, '3' = 0, '4' = 0, '5' = 0))
      }else
      {
        fingerprints_summary_allUnfiltered=rbind(fingerprints_summary_trad,fingerprints_summary_expanded_sos)
        fingerprints_summary_allUnfiltered=rbind(fingerprints_summary_allUnfiltered,fingerprints_summary_expanded_trad)
        
        fingerprints_summary_allUnfiltered = fingerprints_summary_allUnfiltered %>% tidyr::replace_na(list( '1' = 0, '2' = 0, '3' = 0, '4' = 0, '5' = 0))
        
      }
      yColNames = c("1","2","3","4","5")
      fingerprints_summary_all$Pair_ID = fingerprints_summary_all$stimid
    }
  }else
  {
    #FBINoblisBlackBoxTestResponses.csv
    fingerprints <- read.csv("FBINoblisBlackBoxTestResponses.csv", header = TRUE, sep = ",") 
    # create a new column indicating if the prints were mated or nonmated
    # even stimid numbers are mated pairs, odd numbers are nonmated pairs
    
    fingerprints_summary <- fingerprints %>%
      #filter(conclusionType == 1) %>%
      group_by(Pair_ID, Compare_Value, Mating) %>%  #group_by(condition, conclusionType, mated, decision) %>%
      summarize(n = n()) %>% 
      pivot_wider(names_from = Compare_Value, values_from = n)%>%
      select(`Pair_ID`, `Mating`, `Exclusion`, `Inconclusive`, `Individualization`)
    
    fingerprints_summary = fingerprints_summary %>% tidyr::replace_na(list( 'Exclusion' = 0, 'Inconclusive' = 0, 'Individualization' = 0))
   
    fingerprints_summary_nv <- fingerprints %>%
      #filter(conclusionType == 1) %>%
      group_by(Pair_ID, Latent_Value) %>%  #group_by(condition, conclusionType, mated, decision) %>%
      summarize(n = n()) %>% 
      pivot_wider(names_from = Latent_Value, values_from = n)#%>%
     # select(`Pair_ID`, `Latent_Value`)
    
    fingerprints_summary_nv = fingerprints_summary_nv %>% tidyr::replace_na( list('NV' = 0))
    
    
    fingerprints_summary$TotalCompared = fingerprints_summary$Exclusion+
      fingerprints_summary$Inconclusive+
      fingerprints_summary$Individualization
    
    fingerprints_summary$NV = fingerprints_summary_nv$NV
    
    fingerprints_summary <-fingerprints_summary %>%
      filter(TotalCompared >= MinimumNumberOfExaminersInBB)
    
    #rename columns 
    fingerprints_summary$`1` = fingerprints_summary$Exclusion
    fingerprints_summary$`3` = fingerprints_summary$Inconclusive
    fingerprints_summary$`5` = fingerprints_summary$Individualization
    
    fingerprints_summary$`stimid` = as.numeric(factor(fingerprints_summary$Pair_ID))
    fingerprints_summary$`mated` = fingerprints_summary$Mating == "Mates"
    
    fingerprints_summary$`2` = rep.int(0, nrow(fingerprints_summary))
    fingerprints_summary$`4` = rep.int(0, nrow(fingerprints_summary))
    
    fingerprints_summary$`isThree` = rep.int(1, nrow(fingerprints_summary))
    fingerprints_summary$`conclusionType` = rep.int(0, nrow(fingerprints_summary))
    fingerprints_summary$`condition` = rep.int(3, nrow(fingerprints_summary))
    
    yColNames = c("1","2","3","4","5")
    fingerprints_summary_all = fingerprints_summary
    fingerprints_summary_allUnfiltered = fingerprints_summary_all #deal with later
  }
  
  # Example with the movies data.
  OrdModelResultsList = ordinalAnalysis( 
    datFrm = fingerprints_summary_all ,
    yColNames , # column names in data file
    caseIDcolName = "stimid" ,
    compareCases = NULL , # close parenthesis of compareCases list
    hierarchSD = TRUE, #change to true to turn on shrinkage
    doBlackBox = doBlackBox,
    doingT = doingT,
    tDF = tDF
  ) # close parenthesis of function arguments
  
  OrdModelResults = OrdModelResultsList
  OrdModelResults$fingerprints_summary_all=fingerprints_summary_all
  OrdModelResults$fingerprints_summary_allUnfiltered=fingerprints_summary_allUnfiltered
  OrdModelResults$yColNames = yColNames
  #...................................................................................
  
}else
{
  
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
  
  fingerprints_summary_all = OrdModelResults$fingerprints_summary_all
  fingerprints_summary_allUnfiltered = OrdModelResults$fingerprints_summary_allUnfiltered
  yColNames = OrdModelResults$yColNames
  
  graphFileType=c("pdf","png","eps","jpg")[1]
  OrdModelResults$graphFileType = graphFileType
  
}

NForEachPairNumberList = plotParametersAfterMCMC(OrdModelResults, fingerprints_summary_all, 
                                                 fingerprints_summary_allUnfiltered, yColNames)
NForEachPairNumber = NForEachPairNumberList$NForEachPairNumber
NForEachPairNumberUnfiltered = NForEachPairNumberList$NForEachPairNumberUnfiltered

#...................................................................................

OrdModelResults$NForEachPairNumber = NForEachPairNumber
OrdModelResults$NForEachPairNumberUnfiltered =NForEachPairNumberUnfiltered

OrdModelResults$doBlackBox = doBlackBox
OrdModelResults$doingT = doingT

# The following lines are optional, presented as examples of how to do further
# exploration of the output.

# Display the first few lines of the parameter summary matrix:
OrdModelResults$OrdSummaryMat[1:6,]

# Display the top left of the MCMC chain matrix:
OrdModelResults$OrdMcmcMat[1:6,1:10]


if (doBlackBox==TRUE)
{
  if (doingT==FALSE)
  {
  save(OrdModelResults, file= "MCMCParameterSummaryBlackBox.rdata")
  write.csv(OrdModelResults$OrdSummaryMat, file = "MCMCParameterSummaryBlackBox.csv");
  }
  else
  {
    save(OrdModelResults, file= "MCMCParameterSummaryBlackBoxTDist.rdata")
    write.csv(OrdModelResults$OrdSummaryMat, file = "MCMCParameterSummaryBlackBoxTDist.csv");
  }
}else
{
  if (doingT==FALSE)
  {
  save(OrdModelResults, file= "MCMCParameterSummary.rdata")
    write.csv(OrdModelResults$OrdSummaryMat, file = "MCMCParameterSummary.csv");
  }
  else
  {
    save(OrdModelResults, file= "MCMCParameterSummaryTDist.rdata")
    write.csv(OrdModelResults$OrdSummaryMat, file = "MCMCParameterSummaryTDist.csv");
  }
}