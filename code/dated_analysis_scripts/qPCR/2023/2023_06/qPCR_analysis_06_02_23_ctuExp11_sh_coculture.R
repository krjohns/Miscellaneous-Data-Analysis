#read in formatted csv file and perform delta delta cq analysis then make plots
library(data.table)
library(fs)
library(tidyverse)
library(ggpubr) #ggarrange
library(scales) #log2 plotting
library(cowplot) #save_plot
library(monochromeR)

# Plot settings ----
#compress PDFs
tools::compactPDF('mypdfs/',  gs_quality='screen')
###PLOTTING THEME ###
theme_update(text = element_text(size = rel(2)),
             legend.key.size = unit(.1, 'inches'),
             legend.text = element_text(size = rel(2.5)),
             panel.background = element_rect(fill = "white", colour = "black"),
             axis.line = element_blank(),
             axis.text.x = element_text(size=rel(2)),
             axis.text.y = element_text(size=rel(2)),
             axis.title.y = element_text(margin=margin(t = 0, r = 0.3, b = 0, l = 0.5), size=rel(2.5)),
             axis.title.x = element_text(size=rel(2.5)),
             plot.title = element_text(size = rel(2.5), face= "bold", hjust=0.5),
             plot.margin = unit(c(0, 0.1, 0, 0.1), "inches"),
             strip.text.x = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.text.y = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.background = element_rect(size=0.35, color="black"),#trbl
             panel.spacing = unit(0, "lines"))
 

#set directories
wd = getwd()
input_dir = path(wd, "input/data/2023/2023_06")
output_dir = path(wd, "output/2023/2023_06/2023_06_02")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "2023_06_02_ctuExp11_formatted.csv")) %>%
           data.frame() %>%
           select(Well, Gene, Target, Oligo, Bio_Rep, Cq) %>%
           #filter(!(grepl("std", Target))) %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2"),
                  #str combine target & oligo cols
                  Target_oligo = gsub("_$", "",str_c(Target, "_", Oligo))) %>%
           #drop nan and outliers
           filter(!is.nan(Cq)) %>%
           filter(!(Cq > 32)) %>%
           filter(!(Gene == "Gapdh" & Cq > 30))
           
            
dat = qPCR_dat %>%  
           dplyr::select(-Well) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, 
                       values_from = Cq, 
                       names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Gapdh = dat %>%
        filter(Gene == "Gapdh") %>%
        select(Bio_Rep, Target_oligo, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging

dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Tech_Rep_1, Tech_Rep_2,
                  Target_oligo, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Bio_Rep", "Target_oligo")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Gapdh_Avg_Cq` = 2^(-Gapdh_Avg_Cq)) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

# map this to make it more conscise
means = dCq_dat %>%
        group_by(Gene, Target_oligo) %>%
        summarise(rep_means_Cq = mean(Avg_Cq),
                  rep_se_Cq = sd(Avg_Cq)/sqrt(2),
                  rep_means_2_Cq = 2^-(rep_means_Cq),
                  rep_means_dCq = mean(dCq),
                  rep_se_2_dCq = sd(2^-dCq)/sqrt(2),
                  rep_means_2_dCq = 2^-(rep_means_dCq),
                  rep_means_Gapdh = mean(Gapdh_Avg_Cq),
                  rep_se_Gapdh = sd(Gapdh_Avg_Cq)/sqrt(2),
                  rep_means_2_Gapdh = 2^-(rep_means_Gapdh))

dCq_dat = dCq_dat %>%
          full_join(means, by = c("Gene", "Target_oligo")) %>%
          separate(Target_oligo, into = c("Target", "Oligo"), sep = "_", remove = F)
          

dCq_dat = dCq_dat %>%
  mutate(system = "shRNA")

write.csv(as.matrix(dCq_dat), path(wd, "output/2023/2023_06/2023_06_02/2023_06_02_ctuExp11_sh_coculture.csv"))

# plotting fxns  ----
make_Cq_plot_bar = function(gene, dat, colors = Target, dCq = F){
  
  colors = enquo(colors)
  
  #make relevant plot and add mean line
  if(gene == "Gapdh"){
    p = ggplot(dat, 
               aes(y = rep_means_Gapdh, x = Target_oligo)) +
        geom_bar(aes(fill = !!colors),
                 color = "black",
                 position=position_dodge(width=0.5),
                 stat= "identity") +
      geom_point(aes(y = `Gapdh_Avg_Cq`, shape = as.character(Bio_Rep)),
                 color = "black", size = 1.5) +
      geom_text(aes(y = rep_means_Gapdh + 2.5,
                    x=Target_oligo,
                    label = round(rep_means_Gapdh,1)),
                    size = rel(2.5)) +
     geom_errorbar(dat, mapping =aes(ymin= rep_means_Gapdh - rep_se_Gapdh,
                                      ymax= rep_means_Gapdh + rep_se_Gapdh), 
                    width=0.4, position = position_dodge(0.5), 
                    stat = "identity", 
                    colour="black", alpha=0.9, size=0.5) 
  } else if(gene != "Gapdh" & dCq == F){
    p = ggplot(dat %>% filter(Gene == gene), 
               aes(y = rep_means_Cq, x = Target_oligo)) +
      geom_bar(aes(fill = !!colors), 
               color = "black",
               position=position_dodge(width=0.5),
               stat= "identity") +
      geom_point(aes(y = `Avg_Cq`, shape = as.character(Bio_Rep)),
                 color = "black", size = 1.5) +
      geom_text(aes(y = rep_means_Cq +2.5, x=Target_oligo, 
                    label = round(rep_means_Cq,1)), size = rel(2.5)) +
      geom_errorbar(dat, mapping =aes(ymin= rep_means_Cq - rep_se_Cq,
                                      ymax= rep_means_Cq + rep_se_Cq), 
                    width=0.4, position = position_dodge(0.5), 
                    stat = "identity", 
                    colour="black", alpha=0.9, size=0.5)
  } else if(gene != "Gapdh" & dCq == T) {
    p = ggplot(dat %>% filter(Gene == gene), 
               aes(y = rep_means_2_dCq, x = Target_oligo)) +
      geom_bar(aes(fill = !!colors), 
               color = "black",
               position=position_dodge(width=0.5),
               stat= "identity") +
      geom_point(aes(y = `2_dCq`, shape = as.character(Bio_Rep)),
                 color = "black", size = 1.5) +
      geom_errorbar(dat, mapping =aes(ymin= rep_means_2_dCq - (rep_se_2_dCq),
                                      ymax= rep_means_2_dCq + (rep_se_2_dCq)), 
                    width=0.4, position = position_dodge(0.5), 
                    stat = "identity", 
                    colour="black", alpha=0.9, size=0.5) 
  }
  
  p = p +
  labs(x = "", 
       y = ifelse(dCq == F, "Cq", expression(2^-dCq))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(3)),
        strip.text.x = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        plot.title = element_text(size = rel(2.5)),
        text = element_text(size = rel(4)),
        legend.position = "bottom",
        legend.title = element_text(size = rel(3))) 
  

  
  return(p)
  
  }

dat = dCq_dat %>% 
           filter(Gene != "Myd88") %>%
           filter(Gene != "Ctu2") %>%
           filter(Target %in% c("Empty","Scr","Ctu1")) %>%
           mutate(Target = factor(Target, levels = c("Empty",
                                                     "Scr",
                                                     "Myd88",
                                                     "Ctu1", 
                                                     "Ctu2")),
           Target_oligo = factor(Target_oligo, levels = c("Empty_sh11",
                                        str_c("Scr_sh", c("14","15")),
                                        "Myd88_sh21",
                                        str_c("Ctu1_sh", c(1,2,3)),
                                        str_c("Ctu2_sh", c(6,8,17)))))

p1<-make_Cq_plot_bar("Ctu1", dat) + 
     scale_fill_manual(values =c( Scr = "grey",
                                  Untransduced = "white",
                                  Empty = "grey90",
                                  Myd88 = "#68C3A5",
                                  Ctu1 = "#8C9FCA",
                                  Ctu2 = "#F68064")) +
     labs(title = "Ctu1 Expression") 

p2<-make_Cq_plot_bar("Gapdh", dat) + 
      scale_fill_manual(values =c( Scr = "grey",
                                   Untransduced = "white",
                                   Empty = "grey90",
                                   Myd88 = "#68C3A5",
                                   Ctu1 = "#8C9FCA",
                                   Ctu2 = "#F68064")) +
      labs(title = "Gapdh Expression") 

p3 <- make_Cq_plot_bar("Ctu1", dat, dCq = T) + 
      scale_fill_manual(values =c( Scr = "grey",
                               Untransduced = "white",
                               Empty = "grey90",
                               Myd88 = "#68C3A5",
                               Ctu1 = "#8C9FCA",
                               Ctu2 = "#F68064")) +
      labs(title = "Delta Cq (Ctu1 - Gapdh)") 

save_plot(path(output_dir, "Ctu1_KD/Ctu1_cq_only.pdf"), p1)
save_plot(path(output_dir, "Ctu1_KD/Gapdh_cq-Ctu1KD_only.pdf"), p2)
save_plot(path(output_dir, "Ctu1_KD/Ctu1-Gapdh_dcq_only.pdf"), p3)

