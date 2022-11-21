for(t in 1:length(all_trans)){

  if(all_trans[[t]]$sus[1] == df1$inf){

    chain <- paste(all_trans[[t]]$inf[1], all_trans[[t]]$sus[1], sep = "-")

  }

  if(do.call(paste0, all_trans[[t + 2]]) %in% do.call(paste0, all_trans[[t]]))

  for(i in 1:nrow(all_trans[[t]])){

    if(all_trans[[t]]$inf[i] == all_trans[[t]]$sus[1]){
      if(all_trans[[t]]$sus[i] %in% all_trans[[t+1]]$inf){
        chain <- paste(chain, all_trans[[t]]$sus[i], sep = "-")
      }else if(all_trans[[t]]$sus[i] %in% all_trans[[t+1]]$sus){
        chain <- paste(chain, all_trans[[t]]$sus[i], sep = "-")
      }

    }
  }

    if(t <= length(all_trans) )

    for(j in 1:nrow(all_trans[[t+1]])){

      if(all_trans[[t]]$sus[i] == all_trans[[t+1]]$inf[j]){
        chain <- paste(chain, all_trans[[t+1]]$sus[j], sep = "-")
        if(all_trans[[t+1]]$sus[j] == df1$sus){
          break
        }
      }else if(all_trans[[t]]$sus[i] == all_trans[[t+1]]$sus[j]){
        #check whether last element in chain in already from all_trans[[t+1]]
        elements <- str_split(chain, pattern = "-")
        last_element <- elements[[1]][length(elements[[1]])]

        if(last_element == all_trans[[t+1]]$sus[j]){
          break
        }
      }
    }
  }




}


