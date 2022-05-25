#create dag file to run the script to process network analysis in the cluster

#total number of jobs to add in the DAG file
n = 10
params = 1067


jobs <- as.list(rep(1:n))
dagfilename <- "seed_network_10jobs.dag"

#part1 <- unlist(lapply(lines, function(x) paste("JOB SIM", x, " hiv_epimodel_simulation2.sub", sep = "")))
#part2 <- unlist(lapply(lines, function(x) paste("VARS SIM", x, " line_number=", "\"", x, "\"",
#                                                " seed=", "\"", floor(runif(1, min=10000, max=100001)), "\"",
#                                                " dirname=\"params_", x, "/rep_1\"", sep = "")))

part1 <- unlist(lapply(jobs, function(x) paste("JOB SIM", x, " 500mig_1067_deepseq.sub", sep = "")))
part2 <- unlist(lapply(jobs, function(x) paste("VARS SIM", x, " line_number=", "\"", x, "\"",
                                               " seed=", "\"", floor(runif(1, min=10000, max=100001)), "\"",
                                               paste(" dirname=\"params_", params, sep = ""), "/rep_", x, "\"",
                                               paste(" inputname=\"params_", params, sep = ""), "/rep_", x, "/sim_results.tar.gz\"", sep = "")))


write.table(part1, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

write.table(part2, file = dagfilename, quote = FALSE,
            row.names = FALSE, col.names = FALSE,
            append =  TRUE)
