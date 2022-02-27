#script to create windows to use as argument in the script
#to generate alignments using Phyloscanner

window_width <- 250
window_overlap <- 125
start_pos <- 800
end_pos <- 9400

# define windows
windows_start <- seq(start_pos, end_pos, window_overlap)
windows_start <- windows_start[windows_start < end_pos]
windows_end <- windows_start + (window_width - 1)
windows_end <- windows_end[windows_end < end_pos]

#get difference between windows (start and end)
difference <- length(windows_start) - length(windows_end)

#if difference > 0 then remove last elements of windows_end vector
final_wstart <- head(windows_start, -difference)

#create dataframe with start and end position

positions <- data.frame(start = final_wstart, end = windows_end)

write.table(positions, file = "phyloscanner_windows.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)


#to rerun jobs that did not have enough walltime
windows2 <- as.data.frame(windows[c(14:16, 23:32,34:38, 41,43,45:47,50:52,54:67),])
colnames(windows2) <- "V1"
write.table(windows2, file = "phyloscanner_windows_rerun1Oct2021.csv",
            sep = ",", row.names = FALSE, col.names = FALSE)
