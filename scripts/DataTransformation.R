library(readxl)
library(tidyverse)
library(pracma)
library(openxlsx)


# import "vanish_calculation" function,
# which calculates force vanishing
# (minima of QCF/HelFey)
#source("Vanishing.R")

# import "sensing_distance_calculation" function
# which calculates distance of force excess 
# compared with some percent of CCSD QCF maximum
source('SensingDistance.R')

# set data dir
SI_dir <- normalizePath("..")
data_dir <- file.path(SI_dir, "data")
setwd(data_dir)

# set up Minnesota names
Minnesota_names <- c('GAM', 'M05', 'M052X', 'M06', 'M062X', 'M06HF', 'M06L', 'M08-HX',
                     'M08-SO', 'M11', 'M11L', 'MN12L', 'MN12SX', 'MN15',
                     'MN15L', 'N12', 'N12SX', 'revM06', 'revM06L', 'revM11',
                     'SOGGA11X', 'HLE16', 'HLE17', 'MPWLYP1W')

# read data
# dataframe with 2-7A points for each functional for each system
Functional_points <- read_delim('Functional_points.csv', del = ';') %>%
  filter(!method%in%c('MN12L',"HLE16","HLE17")) %>%
  mutate( # interpolate counterpoise correction for MN12SX in HeNe 5.6A QCF
    QCF_Bq = if_else(method == 'MN12SX' & moleculeName == 'HeNe' & R == 5.6, 
                     (QCF_Bq[method == 'MN12SX' & moleculeName == 'HeNe' & R == 5.5] + 
                        QCF_Bq[method == 'MN12SX' & moleculeName == 'HeNe' & R == 5.7]) / 2,
                     QCF_Bq))

Functionals_overall <- read_delim('Functional_overall.csv', del = ';') %>%
  filter(!method%in%c('MN12L',"HLE16","HLE17"))

# additional information about funcs (year, rung)
rung_year_data <- read_delim('FuncData.csv', del=';') %>%
  select(Functional, Rung, Year) %>%
  rename(method=Functional, rung=Rung, year=Year)


# energy_Bq is He + Bq and Bq + He (where Bq is a dummy atom with basis)
# CP - is Counter-Poise corrected descriptor
Functional_points_corrected <- Functional_points %>%
  mutate(Pulay = QCF - HelFey) %>%
  mutate(HelFey_CP = HelFey - HelFey_Bq) %>%
  mutate(QCF_CP = QCF - QCF_Bq) %>%
  mutate(Pulay_CP = QCF_CP - HelFey_CP) %>%
  # introduce Pulay(CP*) (or Pulay delta) = QCF (NoCP) - HelFey (CP)
  mutate(Pulay_delta = QCF - HelFey_CP) %>%
  # translate energy to kcal/mol
  mutate(energy_interaction_CP = (energy - energy_Bq)*627.5095,
         energy_interaction_CP_Ha = energy - energy_Bq) %>% 
  left_join(select(Functionals_overall, method, moleculeName, zero_energy),
             by=c('method', 'moleculeName')) %>%
  mutate(energy_interaction_NoCP = (energy - zero_energy)*627.5095,
         energy_interaction_NoCP_Ha = energy - zero_energy)


DrivenEnergyContrib_A_pure <- function(Functional_points_corrected, 
                                       Functionals_overall_vanished, Angstrom){
  PulayDrivenEnergyContrib_all <- driven_energy_contrib_all_pure(Functional_points_corrected, 
                                                                 Functionals_overall_vanished, 
                                                                 'Pulay_CP', Angstrom)
  HFFDrivenEnergyContrib_all <- driven_energy_contrib_all_pure(Functional_points_corrected, 
                                                               Functionals_overall_vanished, 
                                                               'HelFey_CP', Angstrom) 
  QCFDrivenEnergyContrib_all <- driven_energy_contrib_all_pure(Functional_points_corrected, 
                                                               Functionals_overall_vanished, 
                                                               'QCF_CP', Angstrom)
  DrivenEnergyContrib_Angs <- PulayDrivenEnergyContrib_all %>% 
    left_join(HFFDrivenEnergyContrib_all, by=c('method', 'moleculeName')) %>% 
    left_join(QCFDrivenEnergyContrib_all, by=c('method', 'moleculeName')) %>%
    mutate(is_Minnesota=ifelse(method %in% Minnesota_names,'M','not_M'))
  return (DrivenEnergyContrib_Angs)}

