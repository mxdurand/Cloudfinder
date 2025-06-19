library(data.table)
source("C:/Users/Localadmin_durandma/Dropbox/Work/R/NikyOvercast/NOV_MAESPA_Utils.R")
source("C:/Users/Localadmin_durandma/Dropbox/Work/R/findflecks/findflecks.R")
source("C:/Users/Localadmin_durandma/Dropbox/Work/R/findflecks/findzeroes.R")


# CSthres = 0.7
# OVthres = 0.6
# CStimelimit = 45
# OVtimelimit = 30
# leftovers = TRUE


classifySky <- function(data, CSthres, OVthres, CStimelimit, OVtimelimit, leftovers)
{
  ## 0. Return NA is no data logged
  if(mean(data$PARmes, na.rm = T) < 10){
    verdict = NA
    return(verdict)
    
  }
  
  ## 1. check if overcast
  ovc_flag <- 0 ; time_count <- 0
  for(i in 1:nrow(data))
  {
    TF <- data[i,"PARmes"] < data[i,"PARscl"] * OVthres
    
    # Count duration of overcast period
    if(TF == TRUE & ovc_flag == 0){
      # If new overcast period, turn on flag and initialize count
      ovc_flag <- 1
      time_count <- 0
      prev_time <- data[i,"tsec"]
    } else if (TF == TRUE & ovc_flag == 1){
      # else if continuation, add up time
      current_time <- data[i,"tsec"]
      time_count <- time_count + (current_time - prev_time)
      prev_time <- data[i,"tsec"]
    } else if (TF == FALSE){
      # If period ends of is not overcast, turn off flag
      ovc_flag <- 0
    }
    
    # Check if duration > OVtimelimit, in which case verdict is overcast
    if(time_count / 60 > OVtimelimit){
      verdict = "OV"
      return(verdict)
    }
  }
  
  ## 2. Check clear sky
  csk_flag <- 0 ; time_count <- 0
  for(i in 1:nrow(data))
  {
    TF <- data[i,"PARmes"] > data[i,"PARscl"] * CSthres
    
    # Count duration of clear sky period
    if(TF == TRUE & csk_flag == 0){
      # If new clear sky period, turn on flag and initialize count
      csk_flag <- 1
      time_count <- 0
      prev_time <- data[i,"tsec"]
    } else if (TF == TRUE & csk_flag == 1){
      # else if continuation, add up time
      current_time <- data[i,"tsec"]
      time_count <- time_count + (current_time - prev_time)
      prev_time <- data[i,"tsec"]
    } else if (TF == FALSE){
      # If period ends and is not clear, turn off flag
      csk_flag <- 0
    }
    
    # Check if duration > CStimelimit, in which case verdict is clear sky
    if(time_count / 60 > CStimelimit){
      verdict = "CS"
      return(verdict)
    }
  }
  
  
  ## 3. If not overcast or clear, check partly cloudy
  data$tsec0 <- data$tsec - data$tsec[1]
  z <- findZeros(time = data$tsec0, var = data$PARmes, timeSplit = 10, lim = 1e-4)
  
  dfSf <- NULL ; try({
    dfSf <- findFlecks(time = data$tsec0, var = data$PARmes, zeroes = z,
                       minTime = 0, minAmp = 5, minPdiff = 0, asymmetry = 1/4, trimCV = 0.025,
                       shadeflecks = FALSE, asmMethod = "max", timeSplit = 10, verbose = FALSE)
  }, silent = TRUE)
  dfSh <- NULL ; try({
    dfSh <- findFlecks(time = data$tsec0, var = data$PARmes, zeroes = z,
                       minTime = 0, minAmp = 5, minPdiff = 0, asymmetry = 1/4, trimCV = 0.025,
                       shadeflecks = TRUE, asmMethod = "max", timeSplit = 10, verbose = FALSE)
  }, silent = TRUE)
  
  # Try to detect if PPFD jumps from CS to OV using sunfleck detected
  ppfdIsCS <- FALSE ; ppfdIsOV <- FALSE
  if(!is.null(dfSf) & !is.null(dfSh)){
    
    # If found both sun and shadeflecks:
    for(isf in 1:nrow(dfSf))
    {
      # Select  each sunfleck
      ipeak <- dfSf[isf,"peak"]
      ipeakTime <- dfSf[isf, "peakTime"]
      
      # Check if any peak counts as clear sky
      if(ipeak > data[data$tsec0 == ipeakTime,"PARscl"] * CSthres){
        ppfdIsCS <- TRUE
      }
    }
    for(isf in 1:nrow(dfSh))
    {
      # Select  each shadefleck
      ipeak <- dfSh[isf,"peak"]
      ipeakTime <- dfSh[isf, "peakTime"]
      
      # Check if any peak counts as overcast
      if(ipeak < data[data$tsec0 == ipeakTime,"PARscl"] * OVthres){
        ppfdIsOV <- TRUE
      }
    }
    
    if(ppfdIsCS == TRUE & ppfdIsOV == TRUE){
      verdict = "PC"
      return(verdict)
    }
    
  } else if (!is.null(dfSf)){
    
    # If found only sunflecks:
    for(isf in 1:nrow(dfSf))
    {
      # Select  each sunfleck
      ipeak <- dfSf[isf,"peak"]
      ipeakTime <- dfSf[isf, "peakTime"]
      ibase <- unname(unlist(dfSf[isf,c("baseline1", "baseline2")][which.min(dfSf[isf,c("baseline1", "baseline2")])]))
      ibaseTime <- unname(unlist(dfSf[isf,c("baselineTime1", "baselineTime2")][which.min(dfSf[isf,c("baseline1", "baseline2")])]))
      
      # Check if any peak counts as clear sky and baseline as overcast
      if(ipeak > data[data$tsec0 == ipeakTime,"PARscl"] * CSthres){
        ppfdIsCS <- TRUE
      }
      if(ibase < data[data$tsec0 == ibaseTime,"PARscl"] * OVthres){
        ppfdIsOV <- TRUE
      }
    }
    
    if(ppfdIsCS == TRUE & ppfdIsOV == TRUE){
      verdict = "PC"
      return(verdict)
    }
    
  } else if (!is.null(dfSh)){
    
    # If found only shadeflecks:
    for(isf in 1:nrow(dfSh))
    {
      # Select  each shadefleck (peak and base are inverted)
      ipeak <- dfSh[isf,"peak"]
      ipeakTime <- dfSh[isf, "peakTime"]
      ibase <- unname(unlist(dfSh[isf,c("baseline1", "baseline2")][which.max(dfSh[isf,c("baseline1", "baseline2")])]))
      ibaseTime <- unname(unlist(dfSh[isf,c("baselineTime1", "baselineTime2")][which.max(dfSh[isf,c("baseline1", "baseline2")])]))
      
      # Check if any peak counts as clear sky and baseline as overcast
      if(ibase > data[data$tsec0 == ibaseTime,"PARscl"] * CSthres){
        ppfdIsCS <- TRUE
      }
      if(ipeak < data[data$tsec0 == ipeakTime,"PARscl"] * OVthres){
        ppfdIsOV <- TRUE
      }
    }
    
    if(ppfdIsCS == TRUE & ppfdIsOV == TRUE){
      verdict = "PC"
      return(verdict)
    }
  }
  
  ## 4. Try to decide remaining points
  if(leftovers == TRUE){
    medPARmes <- median(data$PARmes, na.rm = T)
    medPARscl <- median(data$PARscl, na.rm = T)
    if(medPARmes >  medPARscl * CSthres){
      verdict = "CS"
      return(verdict)
    } else if(medPARmes <  medPARscl * OVthres){
      verdict = "OV"
      return(verdict)
    } else if(max(data$PARmes, na.rm = T) > medPARscl * CSthres & min(data$PARmes, na.rm = T) < medPARscl * OVthres){
      verdict = "PC"
      return(verdict)
    }
  }
  
  ## 5. If none the above set as ambiguous
  verdict = "AMB"
  return(verdict)
}

procClassif <- function(file, CSthres, OVthres, CStimelimit, OVtimelimit, leftovers){
  df <- fread(file) ; setDF(df)
  alltimes <- paste(trunc(seq(1,1439,1) / 60, 0), round(((seq(0,1439,1) / 60) - trunc(seq(0,1439,1) / 60, 0)) * 60,0), "00", sep = ":")
  
  # Calculate clear sky 
  idate <- df$date[1]
  dfPAR <- clearpar(localtime = alltimes, date = idate)
  dfPAR$datetime <- as.POSIXct(paste(dfPAR$date, dfPAR$localtime))
  dfPAR$tsec <- tsec(dfPAR$datetime)
  
  # Interpolate on measured data time
  intrp <- approx(x = dfPAR$tsec, y = dfPAR$PARtot, xout = df$tsec)
  PARintrp <- intrp$y
  
  # Find PAR at noon
  solarnoon <- df[df$tsec == intrp$x[which.max(PARintrp)],"Time"]
  PARatnoon <- max(PARintrp, na.rm = T)
  PARatnoon_corrected <- (PARatnoon * 0.9612421) - 153.19768              # 2012
  # PARatnoon_corrected <- (PARatnoon * 0.9612421) - 153.19768 - 241.55339  # 2013
  
  # Scale all calculated PAR based on corrected PAR at noon
  PARscaled <- PARatnoon_corrected * PARintrp / max(PARintrp, na.rm = T)
  
  df <- data.frame("date" = idate, "time" = df$Time, "solarnoon" = solarnoon, "PARatnoon" = PARatnoon_corrected, "PARmes" = df$PAR, "PARmod" = PARintrp, "PARscl" = PARscaled)
  df$tsec <- tsec(df$time)
  rm(dfPAR, intrp, solarnoon, PARatnoon, PARatnoon_corrected, PARintrp, PARscaled)
  
  mor = 25200
  noo = 46800
  aft = 68400
  dfClass <- data.frame()
  
  lday <- list("mor" = df[df$tsec >= mor & df$tsec < noo,], "aft" = df[df$tsec >= noo & df$tsec < aft,])
  for (idt in c("mor", "aft"))
  {
    ilist <- lday[[idt]]
    setDF(ilist)
    
    i <- 1
    while(i <= nrow(ilist))
    {
      # Find start and end of each hour
      start = i
      if(ilist[i,"tsec"] + 60 ** 2 > tail(ilist$tsec,1)){
        end <- which(ilist$tsec == tail(ilist$tsec,1)) # If last hour isnt complete, get what is remaining
      } else {
        end <- which(ilist$tsec == ilist[i,"tsec"] + 60 ** 2) - 1
      }
      
      # In case there is am offset of +/- 1 sec
      s = 0
      while(length(end) == 0)
      {
        end <- which(ilist$tsec == ilist[i,"tsec"] + (60 ** 2) + s) - 1
        s = s - 1
      }
      
      # Subset one hour
      dfH <- ilist[start:end,]
      
      # find classification
      verdict <- classifySky(dfH, CSthres, OVthres, CStimelimit, OVtimelimit, leftovers)
      
      # Fill into data frame
      ihour <- paste0(substring(dfH[trunc(nrow(dfH)/2),"time"], 12, 13), ":00:00")
      df0 <- data.frame("date" = idate, "hour" = ihour, "period" = idt, "sky" = verdict)
      dfClass <- rbind(dfClass, df0)
      
      i <- end + 1 # to start next hour
    }
  }
  return(list(df, dfClass))
}









