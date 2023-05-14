setwd("C:/Rafael/BH/Piecewise")
install.packages("reshape")
require("plyr")
require("dplyr")
require(tidyverse)
require(segmented)
require(formattable)
require(reshape)
require(viridis)
require(wbs)
require(ecp)

######################################################################
###### create picewise linear functions depending on cumulative events
######################################################################
# Create Piecewise Linear models of cumulative sum of your data to identify breakpoints.
{
#'
#' @param dates dd/mm/yyyy format, the dates on which events occurred
#' @param counts The count for each date to generate cumulative sum
#' @param max_pieces The maximum number of breakpoints to consider in analysis, be wary of using too many break points
#' @return A list containing the Aggregated dataset, breakpoints found, gradients between break points and the piecewise linear models generated.
#'
#' @export

piecewise_linear <- function(dates, #The dates that events occurred
                             counts = NULL, #Count of each variable
                             max_pieces = 3){ #max number of breakpoints
  
  if(is.null(counts)){
    counts = rep(1, length(dates)) #If no counts supplied assume counting events
  }
  
  #Now going to run piecewise linear regression for the desired numbers of break points
  
  dates <- as.Date(dates, format='%d/%m/%Y', origin ="01/01/1900") #convert to datetime type
  
  df <- data.frame(dates = dates, counts = counts)
  
  df <- df[order(df$dates),] #order by date
  
  df.agg <- aggregate(df$counts, by=list(dates=df$dates), FUN=sum) #add together counts for events with same date
  names(df.agg)[2] <- 'counts'
  df.agg <- df.agg %>% #this adds a datapoint for each day and takes the cumulative sum
    #complete(dates = seq.Date(as.Date("01/01/1997", format='%d/%m/%Y', origin ="01/01/1900"),
    complete(dates = seq.Date(as.Date("01/01/2010", format='%d/%m/%Y', origin ="01/01/1900"),
                              #min(dates), 
                              #max(dates), 
                              as.Date("20/03/2020", format='%d/%m/%Y', origin ="01/01/1900"),
                              by="day"), fill=list(counts=0))
  df.agg$cum <- cumsum(df.agg$counts)
  
  df.agg$index <- 1:length(df.agg$dates) #Add index for the regression
  
  #Define linear model for counts versus index
  df.lm<-lm(cum~index,data=df.agg)
  
  #Define lists to hold the resulting models and summary statistics
  pwl.models <- vector(mode='list', length=max_pieces)
  summary <- vector(mode='list', length=max_pieces)
  gradients <- vector(mode='list', length=max_pieces)
  break_points <- vector(mode='list', length=max_pieces)
  
  for(pieces in 1:max_pieces){
    o<-segmented(df.lm,
                 seg.Z=~index,
                 model = TRUE,
                 npsi=pieces)
    
    pwl.models[[pieces]] <- o
    summary[[pieces]] <- summary(o)
    
    model = o
    break_point_indices <- c(1) %>%
      append(round(model$psi[,2])) %>%
      append(length(model$fitted.values))
    
    break_points[[pieces]] = df.agg$dates[break_point_indices]
    
    gradients_vec <- vector(length=pieces+1)
    for(piece in 1:length(gradients_vec)){
      x0 <- break_point_indices[piece]
      x1 <- break_point_indices[piece+1]
      y0 <- model$fitted.values[break_point_indices][piece]
      y1 <- model$fitted.values[break_point_indices][piece+1]
      gradients_vec[piece] <- (y1-y0)/(x1-x0)
    }
    gradients[[pieces]] <- gradients_vec
    
  }
  
  get_dates <- function(i){
    ldply(break_points[i], as.Date.character, format="%Y-%m-%d")
  }
  
  columns <- c('One Breakpoint', 'Two Breakpoints', 'Three Breakpoints', 'Four Breakpoints', 'Five Breakpoints')
  
  table_break_points <- ldply(1:max_pieces, get_dates)
  names(table_break_points) <- NULL
  table_break_points <- t(table_break_points)
  table_break_points <- as.data.frame(table_break_points)
  names(table_break_points) <- columns[1:max_pieces]
  null_formatter <- formatter("span", style = x ~ style(color = ifelse(is.na(x), 'white', 'black')))
  table_break_points <- formattable(data.frame(table_break_points),
                                    list(area(col=names(table_break_points)) ~ null_formatter),
                                    col.names=columns[1:max_pieces])
  
  
  gradients.copy <- gradients
  null_formatter <- formatter("span", style = x ~ style(color = ifelse(is.na(x), 'white', 'black')))
  for(piece in 1:max_pieces){
    gradients.copy[[piece]] <- gradients.copy[[piece]][1:(max_pieces+1)]
  }
  gradients.df <- as.data.frame(gradients.copy)
  names(gradients.df) <- columns[1:max_pieces]
  table.gradients <- formattable(gradients.df,
                                 list(area(col=names(gradients.df)) ~ null_formatter),
                                 col.names=columns[1:max_pieces])
  
  output = list(
    'aggregated.data' = df.agg,
    'break.points' = break_points,
    'table.break.points' = table_break_points,
    'model' = pwl.models,
    'gradients' = gradients,
    'table.gradients' = table.gradients
  )
  
  return(output)
  
}
}

# A function to visulaise output of piecewise_linear().
{
#' @param pwlf The dictionary output from piecewise_linear().
#' @param breakpoints_to_plot A vector of numbers of breakpoints to plot.
#' @param legend.position A valid legend position for ggplot2, use 'none' to have no legend.
#' @param breakpoint_lines Bool, whether or not to include vertical lines corresponding to identified break points.
#'
#' @return A ggplot2 object showing the cumulative sum of your data with pieceise linear models fitted.
#'
#' @export

piecewise_linear_plotter <- function(pwlf, #this is the output from my_piecewise_linear()
                                     breakpoints_to_plot=NULL, #Which numbers of breakpoints to plot, vector
                                     legend.position=c(.85,.25), #where to position legend, 'none' for no legend
                                     breakpoint_lines=TRUE){ #Whether to plot the lines of breakpoints or not
  
  colours <- c('#332288', '#44AA99', '#117733', '#CC6677', '#AA4499') #Chosen to be distinct non-primary colours
  
  model <- pwlf$model
  breakpoints <- pwlf$break.points
  data <- pwlf$aggregated.data[,c('dates', 'cum')]
  
  for(pieces in breakpoints_to_plot){
    newcol <- as.data.frame(predict(model[[pieces]]))
    names(newcol) <- paste(toString(pieces),'breakpoints',sep=' ')
    data <- cbind(data,newcol)
  }
  names(data)[2] <- 'True Curve'
  
  data <- melt(data, id.vars='dates')
  
  plot <- ggplot(data=data, aes(x=dates, y=value, color=variable))+
    scale_color_manual(values=colours)+
    geom_line(size=1.2)+#, linetype='dashed', alpha=1)+
    labs(x='',
         y='',
         title='',
         color='legend')+
    theme(legend.position=legend.position,
          legend.background=element_blank(),
          legend.title=element_blank(),
          legend.key = element_blank())
  if(breakpoint_lines){
    color_number=2
    for(number_breakpoints in breakpoints_to_plot){
      for(breakpoint in breakpoints[[number_breakpoints]][-length(breakpoints[[number_breakpoints]])][-1]){
        plot <- plot+geom_vline(xintercept=breakpoint, color=colours[color_number], size=0.9, linetype='dashed')
      }
      color_number <- color_number+1
    }
    
  }
  
  return(plot)
  
}
}

#### Use data for Boko HAram and analyse directionalit
BH <- read.csv("BH_data.csv")
BH <- BH[BH$dateNum >= 40179, ]
BHResponsable <- BH$CLASSIFIED %in% c("Boko Haram", "Battle")
FiltDF <- BH[BHResponsable, ] ### events by Boko Haram
FiltDF <- BH[!BHResponsable, ] ### events againts BH
DF1 <- FiltDF
for (pieces in 1:20){
  max_pieces = pieces
  M <- try(piecewise_linear(DF1$date, max_pieces = max_pieces), silent = TRUE)
  Correct <- try(M$break.points[[1]][1]>1, silent = TRUE)
  if(Correct == TRUE){
    rateMatrix <- matrix(rep(0, max_pieces * length(M$aggregated.data$dates)), ncol = max_pieces)
    for (iter in 1:max_pieces){
      BreakPoints <- iter
      M$rate <- rep(M$gradients[[BreakPoints]][1], length(M$aggregated.data))
      for (k in (BreakPoints + 2):2){
        u <- which(M$aggregated.data$dates <= M$break.points[[BreakPoints]][k])
        M$rate[u] <- M$gradients[[BreakPoints]][k-1]
      }
      rateMatrix[, iter] <- M$rate
    }
    rateMatrix <- as.data.frame(rateMatrix)
    quantile1 <- function(x){quantile(x, probs = 0.25)}
    quantile9 <- function(x){quantile(x, probs = 0.75)}
    rateMatrix$Q.1 <- apply(rateMatrix, 1, quantile1)
    rateMatrix$Q.9 <- apply(rateMatrix, 1, quantile9)
    rateMatrix$MeanR <- apply(rateMatrix, 1, mean)
    rateMatrix$date <- M$aggregated.data$dates
    #FileName <- paste("BHDirectionality/Results_BH_Directionality_pieces_", max_pieces ,".RData", sep = "")
    FileName <- paste("BHDirectionality/Results_BHMNJTF_Directionality_pieces_", max_pieces ,".RData", sep = "")
    save(rateMatrix, file = FileName)
    cat(pieces, "\n")
  }
}
write.csv(rateMatrix, file = "BHDirectionality.csv", row.names = FALSE)

#### SPECIALISATION
{
#### divide by type of event
BH <- read.csv("BH_data.csv")
BH <- BH[BH$dateNum >= 40179, ]
ToDrop <- c("Peaceful protest",
            "Protest with intervention",
            "Protests",
            "Agreement",
            "Headquarters or base established",
            "Non-violent transfer of territory",
            "Other")
ToDropFilt <- BH$sub_event_type %in% ToDrop
BH <- BH[!ToDropFilt,]

#### armed clash
filt <- BH$sub_event_type == "Armed clash"
for (pieces in 1:20){
  max_pieces = pieces
  M <- try(piecewise_linear(BH$date[filt], max_pieces = max_pieces), silent = TRUE)
  Correct <- try(M$break.points[[1]][1]>1, silent = TRUE)
  if(Correct == TRUE){
    rateMatrix <- matrix(rep(0, max_pieces * length(M$aggregated.data$dates)), ncol = max_pieces)
    for (iter in 1:max_pieces){
      BreakPoints <- iter
      M$rate <- rep(M$gradients[[BreakPoints]][1], length(M$aggregated.data))
      for (k in (BreakPoints + 2):2){
        u <- which(M$aggregated.data$dates <= M$break.points[[BreakPoints]][k])
        M$rate[u] <- M$gradients[[BreakPoints]][k-1]
      }
      rateMatrix[, iter] <- M$rate
    }
    rateMatrix <- as.data.frame(rateMatrix)
    quantile1 <- function(x){quantile(x, probs = 0.25)}
    quantile9 <- function(x){quantile(x, probs = 0.75)}
    rateMatrix$Q.1 <- apply(rateMatrix, 1, quantile1)
    rateMatrix$Q.9 <- apply(rateMatrix, 1, quantile9)
    rateMatrix$MeanR <- apply(rateMatrix, 1, mean)
    rateMatrix$date <- M$aggregated.data$dates
    FileName <- paste("ArmedClash_", max_pieces ,".RData", sep = "")
    save(rateMatrix, file = FileName)
    cat(pieces, "\n")
  }
}

#### attack
filt <- BH$sub_event_type == "Attack"
for (pieces in 1:20){
  max_pieces = pieces
  M <- try(piecewise_linear(BH$date[filt], max_pieces = max_pieces), silent = TRUE)
  Correct <- try(M$break.points[[1]][1]>1, silent = TRUE)
  if(Correct == TRUE){
    rateMatrix <- matrix(rep(0, max_pieces * length(M$aggregated.data$dates)), ncol = max_pieces)
    for (iter in 1:max_pieces){
      BreakPoints <- iter
      M$rate <- rep(M$gradients[[BreakPoints]][1], length(M$aggregated.data))
      for (k in (BreakPoints + 2):2){
        u <- which(M$aggregated.data$dates <= M$break.points[[BreakPoints]][k])
        M$rate[u] <- M$gradients[[BreakPoints]][k-1]
      }
      rateMatrix[, iter] <- M$rate
    }
    rateMatrix <- as.data.frame(rateMatrix)
    quantile1 <- function(x){quantile(x, probs = 0.25)}
    quantile9 <- function(x){quantile(x, probs = 0.75)}
    rateMatrix$Q.1 <- apply(rateMatrix, 1, quantile1)
    rateMatrix$Q.9 <- apply(rateMatrix, 1, quantile9)
    rateMatrix$MeanR <- apply(rateMatrix, 1, mean)
    rateMatrix$date <- M$aggregated.data$dates
    FileName <- paste("Attack_", max_pieces ,".RData", sep = "")
    save(rateMatrix, file = FileName)
    cat(pieces, "\n")
  }
}

#### remote attacks
filt <- BH$sub_event_type %in% c("Remote explosive/landmine/IED",
                               "Air/drone strike",
                               "Shelling/artillery/missile attack", 
                               "Grenade")
for (pieces in 1:20){
  max_pieces = pieces
  M <- try(piecewise_linear(BH$date[filt], max_pieces = max_pieces), silent = TRUE)
  Correct <- try(M$break.points[[1]][1]>1, silent = TRUE)
  if(Correct == TRUE){
    rateMatrix <- matrix(rep(0, max_pieces * length(M$aggregated.data$dates)), ncol = max_pieces)
    for (iter in 1:max_pieces){
      BreakPoints <- iter
      M$rate <- rep(M$gradients[[BreakPoints]][1], length(M$aggregated.data))
      for (k in (BreakPoints + 2):2){
        u <- which(M$aggregated.data$dates <= M$break.points[[BreakPoints]][k])
        M$rate[u] <- M$gradients[[BreakPoints]][k-1]
      }
      rateMatrix[, iter] <- M$rate
    }
    rateMatrix <- as.data.frame(rateMatrix)
    quantile1 <- function(x){quantile(x, probs = 0.25)}
    quantile9 <- function(x){quantile(x, probs = 0.75)}
    rateMatrix$Q.1 <- apply(rateMatrix, 1, quantile1)
    rateMatrix$Q.9 <- apply(rateMatrix, 1, quantile9)
    rateMatrix$MeanR <- apply(rateMatrix, 1, mean)
    rateMatrix$date <- M$aggregated.data$dates
    FileName <- paste("RemoteAttack_", max_pieces ,".RData", sep = "")
    save(rateMatrix, file = FileName)
    cat(pieces, "\n")
  }
}

#### suicide bombs
filt <- BH$sub_event_type %in% c("Suicide bomb")
for (pieces in 1:20){
  max_pieces = pieces
  M <- try(piecewise_linear(BH$date[filt], max_pieces = max_pieces), silent = TRUE)
  Correct <- try(M$break.points[[1]][1]>1, silent = TRUE)
  if(Correct == TRUE){
    rateMatrix <- matrix(rep(0, max_pieces * length(M$aggregated.data$dates)), ncol = max_pieces)
    for (iter in 1:max_pieces){
      BreakPoints <- iter
      M$rate <- rep(M$gradients[[BreakPoints]][1], length(M$aggregated.data))
      for (k in (BreakPoints + 2):2){
        u <- which(M$aggregated.data$dates <= M$break.points[[BreakPoints]][k])
        M$rate[u] <- M$gradients[[BreakPoints]][k-1]
      }
      rateMatrix[, iter] <- M$rate
    }
    rateMatrix <- as.data.frame(rateMatrix)
    quantile1 <- function(x){quantile(x, probs = 0.25)}
    quantile9 <- function(x){quantile(x, probs = 0.75)}
    rateMatrix$Q.1 <- apply(rateMatrix, 1, quantile1)
    rateMatrix$Q.9 <- apply(rateMatrix, 1, quantile9)
    rateMatrix$MeanR <- apply(rateMatrix, 1, mean)
    rateMatrix$date <- M$aggregated.data$dates
    FileName <- paste("SuicideBombs_", max_pieces ,".RData", sep = "")
    save(rateMatrix, file = FileName)
    cat(pieces, "\n")
  }
}
}

#### create figure of specialisation By BH
{
  load(paste("BySuicideBombs_", 20,".RData", sep = ""))
  SB <- rateMatrix
  load(paste("ByRemoteAttack_", 14 ,".RData", sep = ""))
  RA <- rateMatrix
  load(paste("ByAttack_", 15 ,".RData", sep = ""))
  Att <- rateMatrix
  load(paste("ByArmedClash_", 18 ,".RData", sep = ""))
  AC <- rateMatrix
  
  colSB <- rev(plasma(20, alpha = 0.5))[12]
  colRA <- rev(plasma(20, alpha = 0.5))[9]
  colAtt <- rev(plasma(20, alpha = 0.5))[5]
  colAC <- rev(plasma(20, alpha = 0.5))[18]
  
  colSBb <- rev(plasma(20, alpha = 1))[12]
  colRAb <- rev(plasma(20, alpha = 1))[9]
  colAttb <- rev(plasma(20, alpha = 1))[5]
  colACb <- rev(plasma(20, alpha = 1))[18]
  cols <- gray.colors(20, start = 0.5, alpha = 0.71)
  png("SpecialisationBHDirectionalityBy.PNG", width = 1200, pointsize = 20)
  par(mar = c(0,0,0,0))
  plot(SB$date, SB$V1, type = "l",
       xlab = "Date", ylab = "Daily event rate",
       xlim = as.Date(c(40179, 43911), format='%d/%m/%Y', origin ="01/01/1900"),
       ylim = c(0, 1.3), col = NA)
  for (k in 1:10){points(SB$date, rep(0, length(SB$date)) + k,
                         type = "l", col = "gray")}
  
  #### SB
  {
    for (iter in 1:20){
      points(SB$date,
             SB[, iter],
             type = "l", lwd = 1, 
             col = cols[iter])
    }
    polygon(c(SB$date, rev(SB$date)),
            c(SB$Q.1, rev(SB$Q.9)),
            col = colSB,
            border = NA)
    points(SB$date, SB$MeanR, type = "l", 
           col = colSBb, 
           lwd = 4)
  }
  
  #### RA
  {
    for (iter in 1:14){
      points(RA$date,
             RA[, iter],
             type = "l", lwd = 1, 
             col = cols[iter])
    }
    polygon(c(RA$date, rev(RA$date)),
            c(RA$Q.1, rev(RA$Q.9)),
            col = colRA,
            border = NA)
    points(RA$date, RA$MeanR, type = "l", 
           col = colRAb, 
           lwd = 4)
  }
  
  #### Att
  {
    for (iter in 1:15){
      points(Att$date,
             Att[, iter],
             type = "l", lwd = 1, 
             col = cols[iter])
    }
    polygon(c(Att$date, rev(Att$date)),
            c(Att$Q.1, rev(Att$Q.9)),
            col = colAtt,
            border = NA)
    points(Att$date, Att$MeanR, type = "l", 
           col = colAttb, 
           lwd = 4)
  }
  
  #### AC
  {
    for (iter in 1:18){
      points(AC$date,
             AC[, iter],
             type = "l", lwd = 1, 
             col = cols[iter])
    }
    polygon(c(AC$date, rev(AC$date)),
            c(AC$Q.1, rev(AC$Q.9)),
            col = colAC,
            border = NA)
    points(AC$date, AC$MeanR, type = "l", 
           col = colACb, 
           lwd = 4)
  }
  dev.off()
}


