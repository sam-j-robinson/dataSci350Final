#Get the log file name
get_log_filename = function(){
  #Prefer to use epoch time.
  epoch_time = toString(as.integer(Sys.time()))
  log_file_name = paste("HW6_", epoch_time, ".log", sep="")
  return(log_file_name)
}

#Create regressive feature from N units ago. Finding specific value
prev_value_n <- function(atomic_table, lag_value){
  regress_data <- sapply(1:length(atomic_table), function(x){
    if(x <= lag_value){
      return(atomic_table[x])
    }
    else {
      return(atomic_table[x-lag_value])
    }
  })
  return(regress_data)
}

#Create regressive features of the mean of previous time amounts
prev_mean_n <- function(atomic_table, lag_value){
  regress_data <- sapply(1:length(atomic_table), function(x){
    if(x <= lag_value+1){
      return(atomic_table[x])
    }
    else {
      return(mean(atomic_table[x-1-lag_value:x-1]))
    }
  })
  return(regress_data)
}

#Create regressive feature from the delta since previous features
prev_delta_n <- function(atomic_table, lag_value){
  regress_data <- sapply(1:length(atomic_table), function(x){
    if(x <= lag_value+1){
      return(atomic_table[x])
    }
    else {
      return(atomic_table[x-lag_value-1] - atomic_table[x-1])
    }
  })
  return(regress_data)
}

#Create a table of regressive features for our predictive model
create_arma_table <- function(start_table, n_back_values, compare_value){
  loginfo("Creating ARMA table")
  nameList = names(start_table)
  for(name in nameList){
    for(i in hours_vector){
      start_table[,paste0(name,"_value_",as.character(i))] <- prev_value_n(start_table[,name], i)
      start_table[,paste0(name,"_mean_",as.character(i))] <- prev_mean_n(start_table[,name], i)
      start_table[,paste0(name,"_delta_",as.character(i))] <- prev_delta_n(start_table[,name], i)
    }
  }
  
  loginfo("Cleaning out non-predictive sale values")
  non_regressive_columns <- sapply(names(start_table), function(x){
    if(startsWith(x, "buy_")){
      return(x)
    }
  })
  return(start_table)
}

#Taken from here
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

test_arma_table <- function(){
  loginfo("Testing arma data creation")
  test_df <- import_and_clean_df("test_data.csv")
  test_arma_df <- create_arma_table(test_df, 1, "spent")
  stopifnot(length(test_arma_df) > 0)
  stopifnot(length(names(test_arma_df)) > length(names(test_df)))
  stopifnot(typeof(test_arma_df) == "list")
}

#Finds an optimal number of pca columns to use for our model.
get_best_pca_model <- function(atomic_comparison, model_pc){
  compareValues <- as.data.frame(t(c(0,0,0,0,0)))
  colnames(compareValues) <- c("pca", "adjRSquared", "RSME", "median_residuals", "AIC")
  
  best_model <- lm(atomic_comparison ~ model_pc$x[,0:1])
  
  for(i in seq(1:length(model_pc$sdev))){
    #Get New Model
    cur_model <- lm(atomic_comparison ~ model_pc$x[,0:i])
    cur_summary <- summary(cur_model)
    
    #Get the units we want to look at for evaluation
    temp_compareValues <-  as.data.frame(t(
      c(i, cur_summary$adj.r.squared, 
        sqrt(mean(cur_summary$residuals^2)), 
        median(cur_summary$residuals), 
        abs(AIC(cur_model)))))
    colnames(temp_compareValues) <- c("pca", "adjRSquared", "RSME", "median_residuals", "AIC")
    
    compareValues <- rbind(compareValues, temp_compareValues)
    
    #Establish multiple benchmarks we want to hit
    if(abs(compareValues[nrow(compareValues), 4]) > abs(compareValues[nrow(compareValues)-1, 4]) &&
       abs(compareValues[nrow(compareValues), 3]) < 100 &&
       compareValues[nrow(compareValues), 2] > .99){
      bes_pcValues = i-1
      best_model <- cur_model
      break
    }
  }
  
  compareValues <- compareValues[2:nrow(compareValues), ]
  p_AdjstRSquared <- ggplot(data=compareValues) +
    geom_line(aes(pca, adjRSquared)) +
    ggtitle("PCA Adjusted R Squared Values") +
    xlab("PCA Variables") + ylab("adj r squared")
  
  p_RSME <- ggplot() +
    geom_line(data=compareValues, aes(pca, RSME)) +
    ggtitle("PCA RSME Values") +
    xlab("PCA Variables") + ylab("RSME")
  
  p_median_residuals <- ggplot() +
    geom_line(data=compareValues, aes(pca, median_residuals)) +
    ggtitle("PCA median residuals") +
    xlab("PCA Variables") + ylab("median residuals")
  
  p_AIC <- ggplot() +
    geom_line(data=compareValues, aes(pca, AIC)) +
    ggtitle("PCA AIC") +
    xlab("PCA Variables") + ylab("AIC")
  
  multiplot(p_AdjstRSquared, p_RSME, p_median_residuals, p_AIC, cols=2)
  
  return(best_model)
}

#Clean output for how long the script has been running
hour_minute_output <- function(seconds){
  cur_seconds = seconds
  hour = 3600
  minute = 60
  duration_hours = 0
  duration_minutes = 0
  
  while(cur_seconds > 60){
    if(cur_seconds >= 3600){
      duration_hours = duration_hours + 1
      cur_seconds = cur_seconds - hour
    }else if(cur_seconds >= 60){
      duration_minutes = duration_minutes + 1
      cur_seconds = cur_seconds - minute
    } else {
      break
    }
  }
  loginfo(paste0("Run Time - ", duration_hours, ":", duration_minutes, ":", cur_seconds))
}

test_import_and_clean <- function(){
  loginfo("Test importing and cleaning data")
  test_df <- import_and_clean_df("test_data.csv")
  stopifnot(length(test_df) > 0)
  stopifnot(typeof(test_df) == "list")
}

#Take our data table and clean it up to be used for the rest of the script
import_and_clean_df <- function(file_name){
  cur_df <- read_csv(file_name)
  loginfo("Cleaning data table")
  cur_df <- cur_df[,sapply(cur_df, function(col) length(unique(col)) > 1)]  
  names(cur_df) <- sapply(names(cur_df), function(name){gsub('-','_', name)})
  #turn time values to factors
  cur_df$saleHour <- factor(cur_df$saleHour)
  cur_df$eventHour <- factor(cur_df$eventHour)  
  cur_df <- cur_df[, -which(names(cur_df) %in% c("spenders", "rpp", "rps", "conversion"))]
  return(cur_df)
}


if(interactive()){
  require(readr)
  require(ggplot2)
  require(stringr)
  require(reshape)
  require(logging)
  setwd("~/Desktop/personal_git/schoolwork/finalProject/")
  
  log_file_name = get_log_filename()
  basicConfig()
  addHandler(writeToFile, file=log_file_name, level='INFO')
  start_time <- proc.time() #Track how long the script takes to run
  hour_minute_output(start_time[['elapsed']])
  
  loginfo("#######################################")
  loginfo("###   Final Project Data Science 350")
  loginfo("###   Predicting Hourly Sales based on previous player behaviour for Empire Z mobile app")
  loginfo("###   By Sam Robinson")
  loginfo("#######################################")
  loginfo("")
  
  #Retrieving and cleaning data
  loginfo("Unit Testing import and clean function")
  test_import_and_clean()
  loginfo("Finish unit import and clean unit test")
  all_items <- import_and_clean_df("sales_empirez_2015-12-1-16_2017-2-25-16.csv")
  
  #How many hours back do we want to look?
  compare_value <- "spent"
  hours_vector <- c(72,96,120,168)
  
  #Create table with regressive features.
  loginfo("Unit testing creating arma data")
  test_arma_table()
  loginfo("Finish arma data creation unit test")
  regressive_values <- create_arma_table(all_items, hours_vector, compare_value)
  
  #Double check that all our factors have multiple levels
  regressive_values <-  regressive_values[, sapply(regressive_values, function(col) length(unique(col))) > 1]
  hour_minute_output((proc.time()-start_time)[['elapsed']])
  
  spent_regressive_matrix <- model.matrix(spent~., data = regressive_values)
  spent_regressive_model <- prcomp(spent_regressive_matrix)
  hour_minute_output((proc.time()-start_time)[['elapsed']])
  
  best_pca_model <- get_best_pca_model(all_items$spent, spent_regressive_model)
  
  best_summary <- summary(best_pca_model)
  best_aic <- AIC(best_pca_model)
  
  loginfo(paste0("Final adjusted rsquared value:", best_summary$adj.r.squared))
  loginfo(paste0("Final RSME Value: ", best_summary$sigma))
  loginfo(paste0("Final median residuals: ", median(best_summary$residuals)))
  loginfo(paste0("Best AIC Value: ", best_aic))
  loginfo("finished running script")
  endTime <- proc.time()-start_time
  hour_minute_output(endTime[['elapsed']])
}