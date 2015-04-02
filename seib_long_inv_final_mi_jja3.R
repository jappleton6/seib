############################################################################################*
# this file provides longitudinal invariance code for the SEI-B data as implemented within R
#       file created on 2014.01.11 by James Appleton
#       last modified on 2014.02.02 by James Appleton
#############################################################################################*

# bracket all code - to stop all execution at point stopifnot() is false
if (1 == 1) {
  
  
  
# setup #######################################################################
  rm(list = ls())
  
  # functions ##########################################################*
  # returns string w/o leading or trailing whitespace
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  # change variable case; df name in quotations to be accepted
  case.cols <- function(x) {
    x.df <- get(x)
    colnames(x.df) <- tolower(names(x.df))
    assign(x,x.df, env = .GlobalEnv)
  }
  
  ################################# switches ###########################*
  # only set ffs to 1 if you need to convert file formats - see just below
  swnm      <- c("install", "import", "efa", "cfa", "invariance", "graph")
  sw        <- c( 0,         1,        0,     0,     1,            0     )
  switches  <- cbind(swnm,sw)
  rm(swnm,sw)
  ######################################################################*
  
  # load packages - 
  # will install if selected under switches above
  
  pkgs <- c("ggplot2", "plyr", "qgraph", "semPlot", "DMwR", "lavaan", "semTools", "mi")
      if (switches[1, 2] == 1) {
        lapply(pkgs, install.packages, dependencies = T)
      }
  lapply(pkgs, require, character.only = T)
  rm(pkgs)
  
  # set working directory
  setwd("E:\\Pinzone_Data")
####################

# import and format
###################
  if (switches[2, 2] == 1) {
    
    sei.b <- read.csv("contrived_id_2013.02.26_seib.csv", header = T)
    demgr <- read.csv("contrived_id_2013.02.26_nm.demgr.eocts.csv", header = T)
    ids <- ids[, 1]
      # import data
        
        # should be no duplicates here
        stopifnot(anyDuplicated(demgr$contrivid)==0)
    
        # merge sei.b w/ demgr to keep only those with legitimate IDs
        sei.b <- merge(sei.b, demgr[,c(1, 3)], by.x = "contrivid", by.y = "contrivid")
          sei.b <- sei.b[,-30] # drop ethnicity which is duplicated
          sei.b <- sei.b[sei.b$contrivid %in% ids, ]
    
        # check dates and create 3 timepoints
        date.distr <- table(sei.b$date)
          
          # save output
          sink("date.distr.txt", append = FALSE)
            date.distr
            round(prop.table(date.distr)*100, 1)
          sink()
    
            date.distr
            round(prop.table(date.distr)*100, 1)
            # main administrations at: 2011.03.24, 2011.04.20, 2011.05.11
              # in interval notation: time 1 = [2011.02.24, 2011.04.06]
              #                       time 2 = [2011.04.07, 2011.04.30]
              #                       time 3 = [2011.05.01, 2011.06.11]
    
    
    
    
    
    
    
        # recode dates; drop cases before 2011.02.24 and after 2011.06.11
        sei.b <- sei.b[sei.b$date >= 20110224 & sei.b$date <= 20110611, ]
          # check new date range 
            table(sei.b$date) # looks good
    
        sei.b$time <- ""
        sei.b[sei.b$date >= 20110224 & sei.b$date <= 20110406, 30] <- "Time 1"
        sei.b[sei.b$date >= 20110407 & sei.b$date <= 20110430, 30] <- "Time 2"
        sei.b[sei.b$date >= 20110501 & sei.b$date <= 20110611, 30] <- "Time 3"
    
        # check recoding
        table(sei.b$date, sei.b$time) # looks good
        class(sei.b$time) # think character is ok but may need factor
    


    
  }

###################
# efa
###################
if (switches[3, 2] == 1) {
    
    
  }

###################
# cfa
###################
if (switches[4, 2] == 1) {
  
      seib.cfa.60 <- read.csv("cfa.csv", header = F)
  
    # name variables
    colnames(seib.cfa.60) <- paste0(rep("q", 27), 1:27)
  
    # assert no NAs (replaced by 9s)
    stopifnot(sum(apply(seib.cfa.60, 1, function(x) {sum(is.na(x))} )) == 0)

    
  cfa.60 <- ' tsr  =~ q16 + q2 + q4 + q9 + q12 + q17 + q20 + q24
              crsw =~ q26 + q11 + q19 + q21 + q25 + q27
              pss  =~ q5 + q3 + q6 + q10 + q18
              fga  =~ q14 + q7 + q13 + q23
              fsl  =~ q15 + q1 + q8 + q22 '
  
  fit.cfa.60 <- cfa(cfa.60, data = seib.cfa.60 ,
                    ordered = c("q1", "q2", "q3", "q4", "q5",
                                "q6", "q7", "q8", "q9", "q10",
                                "q11", "q12", "q13", "q14",
                                "q15", "q16", "q17", "q18",
                                "q19", "q20", "q21", "q22",
                                "q23", "q24", "q25", "q26",
                                "q27"))
  
  # send table of parameters to be estimated to output
  sink("parTableCFA60.txt", append=FALSE, split=FALSE)
  parTable(fit.cfa.60)
  sink()
  

  
  # initial summary of fit
  sink("initialFitCFA60.txt", append=FALSE, split=TRUE)
  summary(fit.cfa.60, fit.measures = TRUE)
  sink()
  
  # modification indices
  MI <- modificationIndices(fit.cfa.60)
  subset(MI, mi > 10)
  
  # standardized coefficients
  parameterEstimates(fit.cfa.60, standardized = TRUE)
    

    
    
  }

###################
# invariance
###################
if (switches[5, 2] == 1) {
    
    
    seib2 <- sei.b[, 3:30]
    
      #subset to those with at least one item completed at all 3 times
      all.3 <- ddply(sei.b[, c(1, 30)], .(contrivid), summarise, nrow(piece))
        all.3 <- all.3[all.3[, 2] == 3, ]
    
        sei.b.mi <- merge(sei.b, all.3, by.x = "contrivid", by.y = "contrivid")
          names(sei.b.mi)[31] <- "responses"
    
            sei.b.mi <- sei.b.mi[, c(1, 30, 3:29)]
    
                names(sei.b.mi)[3:29] <- paste(rep("q", 27), 1:27, sep="")
    
      ##Amelia II multiple imputation of missing data, appears to converge after the 2nd chain?
  
                    ##set.seed to replicate results
              set.seed(4612)  
                    ##Prediction of missing data
              a.out <- amelia(sei.b.mi, m = 5, p2s = 2, idvars = c("contrivid", "time")) 
    
                    ##Check convergence, produce a missmap, check densities of variables -- So many variables make graphs maybe not so useful
              disperse(a.out, dims = 1, m = 5)
              missmap(a.out)
              plot(a.out) ##mean imputations appear to fall within our bounds (1-5) for each question
                    
                    ##save imputations as separate .csv files to visually inspect, some values are over 5 -- does this matter?
                  write.amelia(obj=a.out, file.stem = "SEI_B_mi")
  
    ##old way of multiple imputation begins here
  
#                  info <- mi.info(sei.b.mi)  ## creates an info matrix for the mi
#                  info
#    
#                  # mark contrivid as ID
#                  info <- update(info, "is.ID", list("contrivid" = TRUE))
#                  # check boolean that it worked - should now be TRUE
#                  info$is.ID
#    
                    # missingness map
#                    missing.pattern.plot(sei.b.mi, y.order = TRUE, x.order = TRUE, gray.scale = TRUE)
#    
#
#
#    
#                  sei.b.mimp <- mi(sei.b.mi, info=info, n.imp = 5, n.iter = 30,          ##multiple imputation, performing iterations for 30 minutes
#                                   R.hat = 1.1, max.minutes = 30, rand.imp.method = "bootstrap",
#                                   run.past.convergence = FALSE,
#                                   seed = 123456789, check.coef.convergence = TRUE,
#                                   add.noise = FALSE)
#                  pdf("plots.pdf")
#                  plot(sei.b.mimp)  ##visual analysis (I wasn't sure how to interpret these plots, they are different than the pdf)           
#                      dev.off()
#                      dev.off()
#                      dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off()
#                       dev.off() 
#                       dev.off()
#     
#                   # retrieve imputed datasets
#                   for (i in 1:5) {
#                   assign(paste("mi", i, ".dat", sep = ""), mi.data.frame(sei.b.mimp, m = 1))
#                   }
#     
#               ##  sei.b.mi2 <- sei.b.mi[complete.cases(sei.b.mi),] # remove incomplete cases
#     
#     
#     
#           sei.b.mi[, names(sei.b.mi)[3:29]] <- 
#                   lapply( sei.b.mi[, names(sei.b.mi)[3:29]], ordered)
    
 

          seib.model <- ' tsr  =~ q16 + q2 + q4 + q9 + q12 + q17 + q20 + q24
                          crsw =~ q26 + q11 + q19 + q21 + q25 + q27
                          pss  =~ q5 + q3 + q6 + q10 + q18
                          fga  =~ q14 + q7 + q13 + q23
                          fsl  =~ q15 + q1 + q8 + q22 '
    

            ## calls in each imputation csv and places them into a list
        sei.b.amelia = list()               
        for (i in 1:5) {
          oname <- paste("sei.b.amelia", i, sep="")
         sei.b.amelia[[i]] <- assign(oname, read.csv(paste("SEI_B_mi",i,".csv", sep = ""), header=T))
        }
          

          x <- 1     ##initial loop counter
        while (x != 6) {        ##loop to do invariance analysis on each imputation
              inv.fit <- measurementInvariance(seib.model, data = sei.b.amelia[[x]], ## returns the data frame at each imputation 1-5
                                     group = "time")
              x <- x + 1  ##loop incrementation
        }

  }
  

###################
  # graph
###################
if (switches[6, 2] == 1) {
    
    semPaths(
      fit.cfa.60,
      layout="circle",
      groups=NULL,
      vsize.man=3,
      vsize.lat=6,
      filename="qgraph",
      filetype="pdf",
      residuals=TRUE,
      include=1:12,
      curve=0,
      residSize=0.2,
      onefile=TRUE,
      width=12,
      height=8,
      titles=TRUE)
    
    pdf("semPaths.pdf")
    semPaths(fit.cfa.60, what = "paths", style = "lisrel", layout = "tree2", curve = 1, sizeMan = 3, sizeLat = 7, 
             sizeInt = 3, intercepts = FALSE, residuals = TRUE, thresholds = TRUE, edge.color = "gray", color = "gray", 
             rotation = 2)
    dev.off()
    
    pdf("semPathscoeffs.pdf")
    semPaths(fit.cfa.60, what = "stand", layout = "circle", curve = 1, sizeMan = 5, sizeLat = 7, sizeInt = 3, 
             intercepts = FALSE, residuals = FALSE, thresholds = TRUE, edge.color = "gray", color = "gray", 
             ndigits = 3)
    dev.off()

    

    
    
  }
  
}
  
