#get difference of time of transmission
#for pair of IDs identified as false positives with phyloscanner analysis

library(stringr)
library(dplyr)
library(ggplot2)


get_difference <- function(phylopair){

  timeTrans <- str_split(phylopair$year_trans, pattern = "-")[[1]]



  if(length(timeTrans) == 1){
    #print(phylopair$byGroup)
    View(phylopair)
  }



  if (length(timeTrans) == 2){

    difference_results <- data.frame(difference_timeTrans = as.numeric(timeTrans[2]) - as.numeric(timeTrans[1]))

  }

  if(length(timeTrans) > 2){

    difference_results <- data.frame(difference_timeTrans = as.numeric(timeTrans[length(timeTrans)]) - as.numeric(timeTrans[1]))

  }

  all_results <- cbind(difference_results, phylopair)

  all_results["difference_donor_recip"] <- all_results$st_donor_ID - all_results$st_recip_ID

  return(all_results)

}


#read pairs in which filtering was W >= 80%

FP_W80 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_FP_results_W80_final.RDS")
FP_W80["param_mig"] <- paste(FP_W80$param, FP_W80$mig, sep = "_")
FP_W80["byGroup"] <- c(1:nrow(FP_W80))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)


differences <- FP_W80 %>%
  group_by(byGroup) %>%
  group_modify(~get_difference(.x))

differences$difference_timeTrans <- abs(differences$difference_timeTrans)

#get summary from total chain
summary_by_params <- differences %>%
  group_by(param_mig) %>%
  mutate(mean = mean(difference_timeTrans),
         median = median(difference_timeTrans),
         Min = summary(difference_timeTrans)[1],
         Max = summary(difference_timeTrans)[6],
         quantile1 = summary(difference_timeTrans)[2],
         median1 = summary(difference_timeTrans)[3],
         quantile3 = summary(difference_timeTrans)[5],
         total_rows = n()) %>%
  select(param_mig, mean, median, median1, total_rows, Min, Max, quantile1, quantile3) %>%
  distinct()

#table for order of transmissions
order_table_FP_W80 <- as.data.frame.matrix(table(FP_W80$param_mig,FP_W80$order))
order_table_FP_W80["total"] <- unname(apply(order_table_FP_W80, 1, sum))



differences %>%
  ggplot( aes(x=param_mig, y=difference_timeTrans, fill=param_mig)) +
  geom_boxplot()




#read pairs in which filtering was W >= 1%

FP_W0.01 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_FP_results_W0.01_final.RDS")
FP_W0.01["param_mig"] <- paste(FP_W0.01$param, FP_W0.01$mig, sep = "_")
FP_W0.01["byGroup"] <- c(1:nrow(FP_W0.01))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)

differences_W0.01 <- FP_W0.01 %>%
  group_by(byGroup) %>%
  group_modify(~ get_difference(.x))

differences_W0.01$difference_timeTrans <- abs(differences_W0.01$difference_timeTrans)

#get summary from total chain
summary_by_params_W0.01 <- differences_W0.01 %>%
  group_by(param_mig) %>%
  mutate(mean = mean(difference_timeTrans),
         median = median(difference_timeTrans),
         Min = summary(difference_timeTrans)[1],
         Max = summary(difference_timeTrans)[6],
         quantile1 = summary(difference_timeTrans)[2],
         median1 = summary(difference_timeTrans)[3],
         quantile3 = summary(difference_timeTrans)[5],
         total_rows = n()) %>%
  select(param_mig, mean, median, median1, total_rows, Min, Max, quantile1, quantile3) %>%
  distinct()


#table for order of transmissions
order_table_FP_W0.01 <- as.data.frame.matrix(table(FP_W0.01$param_mig,FP_W0.01$order))
order_table_FP_W0.01["total"] <- unname(apply(order_table_FP_W0.01, 1, sum))






differences_W0.01 %>%
  ggplot( aes(x=param_mig, y=difference_timeTrans, fill=param_mig)) +
  geom_boxplot()




