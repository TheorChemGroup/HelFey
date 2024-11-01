vanish_calculation = function(data, force){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        distance = numeric())

  # iterate over the unique pairs of method and molecule
  for (method_ in unique(data$method)) {
    for (molecule in unique(data$moleculeName)) {
      # filter the data by method and molecule
      subset <- subset(data, method == method_ & moleculeName == molecule)
      # set distance to NA, if there is no force vanishing
      distance <- NA
      for (i in 1:nrow(subset)) {
        if (subset[[force]][i] <= 0) {
          if (i==1){
            distance <- 2
            break}
          else{
          # found a value below zero, 
          # so interpolate linearly the corresponding distance
          x1 <- subset$R[i-1]
          x2 <- subset$R[i]
          y1 <- subset[[force]][i-1]
          y2 <- subset[[force]][i]
          distance <- x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
          break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           distance = distance))
    }
  }
  dist_name <- paste0(force, '_vanish_dist')
  results_renamed <- results %>% rename(!!dist_name := distance)
  return(results_renamed)
}



