#get difference of time of transmission
#for pair of IDs identified as false positives with phyloscanner analysis

library(stringr)
library(dplyr)
library(ggplot2)


get_difference <- function(phylopair){

  timeTrans <- str_split(phylopair$year_trans, pattern = "-")[[1]]



  if(length(timeTrans) == 1){
    difference_results <- data.frame(difference_timeTrans = NA)
    print(count)

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

TN_W80 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_TN_results_W80_final.RDS")
TN_W80["param_mig"] <- paste(TN_W80$param, TN_W80$mig, sep = "_")
TN_W80["byGroup"] <- c(1:nrow(TN_W80))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)


differences <- TN_W80 %>%
  group_by(byGroup) %>%
  group_modify(~get_difference(.x))

differences$difference_timeTrans <- abs(differences$difference_timeTrans)

#get summary from total chain
summary_by_params <- differences %>%
  group_by(param_mig) %>%
  mutate(mean = mean(difference_timeTrans, na.rm =TRUE),
         median = median(difference_timeTrans, na.rm =TRUE),
         Min = summary(difference_timeTrans, na.rm =TRUE)[1],
         Max = summary(difference_timeTrans, na.rm =TRUE)[6],
         quantile1 = summary(difference_timeTrans, na.rm =TRUE)[2],
         median1 = summary(difference_timeTrans, na.rm =TRUE)[3],
         quantile3 = summary(difference_timeTrans, na.rm =TRUE)[5],
         total_rows = n()) %>%
  select(param_mig, mean, median, median1, total_rows, Min, Max, quantile1, quantile3) %>%
  distinct()


#table for order of transmissions
order_table_TN_W80 <- as.data.frame.matrix(table(TN_W80$param_mig,TN_W80$order))
order_table_TN_W80["total"] <- unname(apply(order_table_TN_W80, 1, sum))


#check difference of transmission times for order 1-2_2-3 etc
order1223_W80 <- differences[differences$order == "1-2_2-3"|
                               differences$order == "1-2_2-3_3-4"|
                               differences$order == "1-2_2-3_3-4_4-5",]

summary_by_params_order1223_W80 <- order1223_W80 %>%
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



differences %>%
  ggplot( aes(x=param_mig, y=difference_timeTrans, fill=param_mig)) +
  geom_boxplot()




#read pairs in which filtering was W >= 1%

TN_W0.01 <- readRDS("Analyses/deep_sequencing/Result_analyses/Results/all_TN_results_W0.01_final.RDS")
TN_W0.01["param_mig"] <- paste(TN_W0.01$param, TN_W0.01$mig, sep = "_")
TN_W0.01["byGroup"] <- c(1:nrow(TN_W0.01))
#get maximum number of columns to separate IDs involved in a chain
#from inf to sus individual (pair analysed with phyloscanner)

differences_W0.01 <- TN_W0.01 %>%
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
order_table_TN_W0.01 <- as.data.frame.matrix(table(TN_W0.01$param_mig,TN_W0.01$order))
order_table_TN_W0.01["total"] <- unname(apply(order_table_TN_W0.01, 1, sum))



#check difference of transmission times for order 1-2_2-3 etc
order1223_W0.01 <- differences_W0.01[differences_W0.01$order == "1-2_2-3"|
                                       differences_W0.01$order == "1-2_2-3_3-4"|
                                       differences_W0.01$order == "1-2_2-3_3-4_4-5"|
                                       differences_W0.01$order == "1-2_2-3_3-4_4-5_5-6" |
                                       differences_W0.01$order == "1-2_2-3_3-4_4-5_5-6_6-7" |
                                       differences_W0.01$order == "1-2_2-3_3-4_4-5_5-6_6-7_7-8" |
                                       differences_W0.01$order == "1-2_2-3_3-4_4-5_5-6_6-7_7-8_8-9",]

summary_by_params_order1223_W0.01 <- order1223_W0.01 %>%
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




differences_W0.01 %>%
  ggplot( aes(x=param_mig, y=difference_timeTrans, fill=param_mig)) +
  geom_boxplot()




