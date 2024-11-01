sensing_distance_calculation <- function(data, force, percent){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        distance = numeric())
  
  QCF_min_dist_CCSD <- new.env()

  for (molecule in unique(data$moleculeName)) {
    CCSD_molecule_data = filter(data,
                                moleculeName == molecule, 
                                method == 'CCSD-full')
    QCF_min_dist_CCSD$molecule = min(CCSD_molecule_data$QCF)
    #print(QCF_min_dist_CCSD$molecule)
    for (method_ in unique(data$method)) {
      subset <- subset(data, method == method_ & moleculeName == molecule)
      distance <- 2
      for (i in nrow(subset):1) {
        if (subset[[force]][i] <= QCF_min_dist_CCSD$molecule * percent){
          if (i==nrow(subset)){
            distance <- 7
            break}
          else{
            # found a value above needed force, 
            # so interpolate linearly the corresponding distance
            x1 <- subset$R[i]
            x2 <- subset$R[i+1]
            y1 <- subset[[force]][i]
            y2 <- subset[[force]][i+1]
            distance <- x1 + (QCF_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
            break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           distance = distance))
    }
  }
  colnames(results)[colnames(results) == "distance"] <- paste0(force, '_sensing_distance_', percent * 100)
  return(results)
}



cumulative_sensing_distance_calculation <- function(data, force, percent, is_CP){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        distance = numeric())
  
  energy_min_dist_CCSD <- new.env()
  
  for (molecule in unique(data$moleculeName)) {
    CCSD_molecule_data = filter(data,
                                moleculeName == molecule, 
                                method == 'CCSD-full')
    
    energy_min_dist_CCSD$molecule = min(CCSD_molecule_data$energy_interaction_CP_Ha)
    for (method_ in unique(data$method)) {
      subset <- subset(data, method == method_ & moleculeName == molecule)
      subset <- subset %>%
        mutate(integrative_force = rev(cumtrapz(subset$R, rev(subset[[force]]))))
      #print(subset$integrative_force)
      distance <- 2
      for (i in nrow(subset):1) {
        if (subset$integrative_force[i] <= energy_min_dist_CCSD$molecule * percent){
          if (i==nrow(subset)){
            distance <- 7
            break}
          else{
            # found a value above needed force, 
            # so interpolate linearly the corresponding distance
            x1 <- subset$R[i]
            x2 <- subset$R[i+1]
            y1 <- subset$integrative_force[i]
            y2 <- subset$integrative_force[i+1]
            distance <- x1 + (energy_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
            break}
      #distance <- 2
      #for (i in 1:nrow(subset)) {
      #  if (subset$integrative_force[i] <= energy_min_dist_CCSD$molecule * percent){
      #    if (i==nrow(subset)){
      #      distance <- 7
      #      break}
      #    else{
      #      # found a value above needed force, 
      #      # so interpolate linearly the corresponding distance
      #      x1 <- rev(subset$R)[i-1]
      #      x2 <- rev(subset$R)[i]
      #      y1 <- subset$integrative_force[i-1]
      #      y2 <- subset$integrative_force[i]
      #      distance <- x1 + (energy_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
      #      break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           distance = distance))
    }
  }
  colnames(results)[colnames(results) == "distance"] <- paste0(force, '_cumulative_sensing_distance_', percent * 100)
  return(results)
}


driven_energy_contrib <- function(data, data_overall, force, Angs){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        value = numeric())
  
  energy_min_dist_CCSD <- new.env()
  vanish_dist_CCSD <- new.env()
  
  for (molecule in unique(data$moleculeName)) {
    CCSD_molecule_data = filter(data,
                                moleculeName == molecule, 
                                method == 'CCSD-full')
    CCSD_molecule_data_overall = filter(data_overall,
                                 moleculeName == molecule, 
                                 method == 'CCSD-full')
    
    # add if new variable is True then energy_interaction_NoCP_Ha else
    # energy_interaction_NoCP_Ha

    #energy_min_dist_CCSD$molecule = min(CCSD_molecule_data$energy_interaction_CP_Ha)
    #energy_min_dist_CCSD$molecule = filter(CCSD_molecule_data, R == Angs)
    vanish_dist_CCSD$molecule = CCSD_molecule_data_overall$QCF_CP_vanish_dist

    #print(vanish_dist_CCSD$molecule)
    
    for (method_ in unique(data$method)) {
      subset <- subset(data, method == method_ & moleculeName == molecule)
      subset <- subset %>%
        mutate(integrative_force = rev(cumtrapz(subset$R, rev(subset[[force]]))))
      #print(subset$integrative_force)
      for (i in nrow(subset):1) {
        if (subset$R[i] <= vanish_dist_CCSD$molecule){
          if (i==nrow(subset)){
            value <- subset$integrative_force[1]
            break}
          else{
            # found a value above needed force, 
            # so interpolate linearly the corresponding distance
            x1 <- subset$R[i]
            x2 <- subset$R[i+1]
            y1 <- subset$integrative_force[i]
            y2 <- subset$integrative_force[i+1]
            value <- y1 + (y2-y1)*(vanish_dist_CCSD$molecule-x1)/(x2-x1)
            #distance <- x1 + (energy_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
            break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           value = value/energy_min_dist_CCSD$molecule*100))
    }
  }
  colnames(results)[colnames(results) == "value"] <- paste0(force, '_driven_energy')
  return(results)
}


driven_energy_contrib_all <- function(data, data_overall, force, angs=NULL){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        value = numeric())
  
  energy_min_dist_CCSD <- new.env()
  vanish_dist <- new.env()
  
  for (molecule in unique(data$moleculeName)) {
    CCSD_molecule_data = filter(data,
                                moleculeName == molecule, 
                                method == 'CCSD-full')
    
    
    # add if new variable is True then energy_interaction_NoCP_Ha else
    # energy_interaction_NoCP_Ha
    
    energy_min_dist_CCSD$molecule = filter(CCSD_molecule_data, R==angs)$energy_interaction_CP_Ha
    

    #print(vanish_dist_CCSD$molecule)
    
    for (method_ in unique(data$method)) {
      molecule_data_overall = filter(data_overall,
                                      moleculeName == molecule, 
                                      method == method_)
      method_molecule_data = filter(data,
                                  moleculeName == molecule, 
                                  method == method_)
      vanish_dist$molecule = molecule_data_overall$QCF_CP_vanish_dist
      #method_energy_min$molecule = min(method_molecule_data$energy_interaction_CP_Ha)
      
      if (!is.null(angs)) {
        vanish_dist$molecule <- angs
      }
      if (is.na(vanish_dist$molecule)){
        results <- rbind(results, data.frame(method = method_,
                                             moleculeName = molecule,
                                             value = NA))
        next}
      subset <- subset(data, method == method_ & moleculeName == molecule)
      subset <- subset %>%
        mutate(integrative_force = rev(cumtrapz(subset$R, rev(subset[[force]]))))
      #print(subset$integrative_force)
      for (i in nrow(subset):1) {
        if (subset$R[i] <= vanish_dist$molecule){
          if (i==nrow(subset)){
            value <- subset$integrative_force[nrow(subset)]
            break}
          else{
            # found a value above needed force, 
            # so interpolate linearly the corresponding distance
            x1 <- subset$R[i]
            x2 <- subset$R[i+1]
            y1 <- subset$integrative_force[i]
            y2 <- subset$integrative_force[i+1]
            value <- y1 + (y2-y1)*(vanish_dist$molecule-x1)/(x2-x1)
            #distance <- x1 + (energy_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
            break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           value = value/energy_min_dist_CCSD$molecule*100))
    }
  }
  colnames(results)[colnames(results) == "value"] <- paste0(force, '_driven_energy')
  return(results)
}


driven_energy_contrib_all_pure <- function(data, data_overall, force, angs=NULL){
  
  # create an empty data frame to store the results
  results <- data.frame(method = character(),
                        moleculeName = character(),
                        value = numeric())
  
  energy_min_dist_CCSD <- new.env()
  vanish_dist <- new.env()
  
  for (molecule in unique(data$moleculeName)) {
    for (method_ in unique(data$method)) {
      molecule_data_overall = filter(data_overall,
                                     moleculeName == molecule, 
                                     method == method_)
      method_molecule_data = filter(data,
                                    moleculeName == molecule, 
                                    method == method_)
      vanish_dist$molecule = molecule_data_overall$QCF_CP_vanish_dist
      #method_energy_min$molecule = min(method_molecule_data$energy_interaction_CP_Ha)
      
      if (!is.null(angs)) {
        vanish_dist$molecule <- angs
      }
      if (is.na(vanish_dist$molecule)){
        results <- rbind(results, data.frame(method = method_,
                                             moleculeName = molecule,
                                             value = NA))
        next}
      subset <- subset(data, method == method_ & moleculeName == molecule)
      subset <- subset %>%
        mutate(integrative_force = rev(cumtrapz(subset$R, rev(subset[[force]]))))
      #print(subset$integrative_force)
      for (i in nrow(subset):1) {
        if (subset$R[i] <= vanish_dist$molecule){
          if (i==nrow(subset)){
            value <- subset$integrative_force[nrow(subset)]
            break}
          else{
            # found a value above needed force, 
            # so interpolate linearly the corresponding distance
            x1 <- subset$R[i]
            x2 <- subset$R[i+1]
            y1 <- subset$integrative_force[i]
            y2 <- subset$integrative_force[i+1]
            value <- y1 + (y2-y1)*(vanish_dist$molecule-x1)/(x2-x1)
            #distance <- x1 + (energy_min_dist_CCSD$molecule * percent - y1) * (x2 - x1) / (y2 - y1)
            break}
        }
      }
      # add the results to the data frame
      results <- rbind(results, data.frame(method = method_,
                                           moleculeName = molecule,
                                           value = value))
    }
  }
  colnames(results)[colnames(results) == "value"] <- paste0(force, '_driven_energy')
  return(results)
}