## install these packages before

library(readxl)
library(tidyverse)
library(pracma)
library(openxlsx)
library(ggrepel)
library(svglite)
library(glue)


SI_dir <- normalizePath("..")
scripts_dir <- file.path(SI_dir, "scripts")
setwd(scripts_dir)
source("DataTransformation.R")
plots_dir <- file.path(SI_dir, "plots")
setwd(plots_dir)


Ha_to_ev = 27.2114
Ha_bohr_to_ev_A = 51.42208619083232

FuncSpread <- tibble(method=c("PBE0", "B3LYPV1R", "M11", "M11L", "revM11", "M06L", 
                              "revM06L", "revTPSS", "VSXC", "M062X", "SLATER", "MN15"),
                     column=c(1,1,1,1,1,2,2,1,2,2,1,2),
                     row=c(1,1,2,2,2,1,1,1,2,2,1,2))

SelFuncData <- Functional_points_corrected %>% 
  select(method,moleculeName,R,QCF_CP,HelFey_CP,Pulay_CP) %>% 
  filter(moleculeName=="NeNe",
         method%in%FuncSpread$method) %>% 
  gather("Force","ForceVal",QCF_CP,HelFey_CP,Pulay_CP) %>% 
  mutate(Force=factor(Force,levels=c("QCF_CP","HelFey_CP","Pulay_CP")),
         method=factor(method, levels=FuncSpread$method))

SelFuncData %>% distinct(method)

SelCCSDData <- Functional_points_corrected %>% 
  select(method,moleculeName,R,QCF_CP,HelFey_CP,Pulay_CP) %>% 
  filter(moleculeName=="NeNe", method%in%c("CCSD-full")) %>% 
  gather("Force","ForceVal",QCF_CP,HelFey_CP,Pulay_CP) %>% 
  mutate(Force=factor(Force,levels=c("QCF_CP","HelFey_CP","Pulay_CP")))

Figure1 <- SelFuncData %>% 
  left_join(FuncSpread) %>% 
  filter(Force!="Pulay_CP",
         method!="SLATER") %>% 
  ggplot(aes(y=ForceVal* Ha_to_ev,x=R,color=method, linetype=Force)) +
  geom_hline(yintercept = 0,color="grey80",size=1.2) +
  geom_line(data=SelCCSDData %>% filter(Force!="Pulay_CP"),color="black",size=1.25) +
  geom_line(size=1.1) +
  facet_grid(column~row,scales="free_y",space="free_y") +
  coord_cartesian(xlim=c(2.8,4.5),
                  ylim=c(NA,0.00035*Ha_to_ev)) +
  scale_x_continuous(name="Interatomic distance, Å") +
  scale_y_continuous(name="Force value, eV/Å") +
  scale_color_manual(values=c("#e31a1c", "#fdbf6f", "#ff7f00","#fb9a99",
                              "#a6cee3", "#1f78b4","#349AB6","#68E52D", "#2FBE2F",
                              "#984ea3", "#cab2d6", "#FF90FF", "#ffff99","6a3d9a"),
                     breaks=c("PBE0", "B3LYPV1R","revTPSS", "SLATER", 
                              "M11L", "M11", "revM11", 
                              "M06L", "revM06L", 
                              "M062X", "VSXC", "MN15"),
                     labels=c("PBE0", "B3LYP", "revTPSS", "SLATER", 
                              "M11L", "M11", "revM11",
                              "M06-L", "revM06-L", 
                              "M062X", "VSXC","MN15"),
                     name="DFT functional") +
  scale_linetype(name="Force type",
                 labels=c("Quantum\nChemical\nForce (QCF)","Hellmann-\nFeynman\nForce (HFF)"),
                 #labels=c("QCF","HFF"),
                 limits=c("QCF_CP","HelFey_CP")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid=element_blank())
Figure1
ggsave("Figure1.png",Figure1, width=18, heigh=12, units="cm",dpi=2000)


PureDrivenE_many <- tibble(coff=seq(2,4,0.2)) %>% 
  group_by(coff) %>% 
  do(DrivenEnergyContrib_A_pure(Functional_points_corrected, 
                                Functionals_overall, .$coff)) %>% 
  unnest()

PureDrivenE_many

RefEn_many <- PureDrivenE_many %>% 
  filter(method=="CCSD-full") %>% 
  select(moleculeName,coff,QCF_CP_driven_energy) %>% 
  rename(QCF_CP_ref=QCF_CP_driven_energy)
RefEn_many

QCFerr_many <- PureDrivenE_many %>% 
  filter(method!="CCSD-full") %>% 
  select(coff,method,moleculeName,QCF_CP_driven_energy) %>% 
  left_join(RefEn_many) %>% 
  rowwise() %>% 
  mutate(QCFerr0 = abs(QCF_CP_driven_energy-QCF_CP_ref)) %>% 
  group_by(coff,moleculeName) %>% 
  mutate(QCFerr=-QCFerr0/QCF_CP_ref) %>% 
  ungroup() %>% 
  distinct(coff,method,moleculeName,QCFerr) %>% 
  group_by(coff,method) %>% 
  mutate(maxQCFerr=max(QCFerr),
         IsMinnesota=ifelse(method%in%Minnesota_names,"Minnesota","Other")) %>%
  arrange(desc(maxQCFerr)) 
QCFerr_many

PulayHelp_many <- PureDrivenE_many %>% 
  filter(method!="CCSD-full") %>% 
  select(coff,method,moleculeName,QCF_CP_driven_energy,HelFey_CP_driven_energy) %>% 
  left_join(RefEn_many) %>% 
  rowwise() %>% 
  mutate(PulayHelp0 = abs(HelFey_CP_driven_energy-QCF_CP_ref)-abs(QCF_CP_driven_energy-QCF_CP_ref)) %>% 
  rowwise() %>% 
  mutate(PulayHelp=PulayHelp0/abs(QCF_CP_ref)) %>% 
  ungroup() %>% 
  distinct(coff,method,moleculeName,PulayHelp) %>% 
  group_by(coff,method) %>% 
  mutate(maxPH=max(PulayHelp),
         maxSystem=moleculeName[which.max(PulayHelp)],
         IsMinnesota=ifelse(method%in%Minnesota_names,"Minnesota","Other")) %>%
  arrange(desc(maxPH))
PulayHelp_many
PulayHelp_many %>% arrange(desc(abs(maxPH)))
PulayHelp_many %>% spread(moleculeName,PulayHelp) %>% filter(coff==3.2,
                                                             method%in%c("M11L","M06L","B3LYPV1R","M052X"))


PulHelpPlot_MS <- left_join(PulayHelp_many %>% select(coff,method,maxPH,IsMinnesota,maxSystem) %>% 
                              distinct() %>% rename(moleculeName=maxSystem),
                            QCFerr_many %>% select(coff,method,moleculeName,QCFerr,IsMinnesota))  %>% 
  filter(coff==3.2) %>% 
  ggplot(aes(x=100*QCFerr,y=100*maxPH,color=IsMinnesota,fill=IsMinnesota)) +
  geom_vline(xintercept = 1,color="grey80") +
  geom_hline(yintercept = 0,color="grey80") +
  geom_hline(yintercept = 10,color="pink") +
  geom_hline(yintercept = -10,color="pink") +
  geom_point(size=1) +
  geom_text_repel(aes(label=method),size=2.5,force_pull = 40,force = 0.5) +
  scale_color_manual(values=c('black', '#7CAE00',"#e41a1c", "#377eb8"),
                     breaks=c(0, 1, "Minnesota", "Other"),
                     labels=c('CCSD-full', 'Dispersionless', 'Minnesota', 'Other DFT')) +
  scale_fill_manual(values=c('black', '#7CAE00',"#e41a1c", "#377eb8"),
                    breaks=c(0, 1, "Minnesota", "Other"),
                    labels=c('CCSD-full', 'Dispersionless', 'Minnesota', 'Other DFT')) +
  scale_x_log10(name="Absolute error in QCF-driven energy at 3.2 Å, %") +
  scale_y_continuous(name="Maximum Pulay helpness at 3.2 Å, %")+
  #                breaks=c(1,5,10,25,50,100)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid=element_blank())
PulHelpPlot_MS
ggsave("Figure2_PulHelpPlot.png",PulHelpPlot_MS,width=9,heigh=12,units="cm",dpi=1000)
ggsave("Figure2_PulHelpPlot.svg",PulHelpPlot_MS,width=9,heigh=12,units="cm")
