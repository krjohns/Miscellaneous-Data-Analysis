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
             axis.title.y = element_text(margin=margin(t = 0, r = 0.8, b = 0, l = 0.8), size=rel(2.5)),
             axis.title.x = element_text(size=rel(2.5)),
             plot.title = element_text(size = rel(2.5), face= "bold", hjust=0.5),
             plot.margin = unit(c(0, 0.05, 0, 0.05), "inches"),
             strip.text.x = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.text.y = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.background = element_rect(linewidth=0.35, color="black"),#trbl
             panel.spacing = unit(0.5, "lines"))
 

#set directories
wd = getwd()
input_dir = path(wd, "input/data/2024/2024_04")
output_dir = path(wd, "output/2024/2024_04/20240423")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "20240418_snoRNA_bmdc_test1.csv")) %>%
           data.frame() %>%
           dplyr::select(-Type) %>%
           #combine primer info with gene info and drop primer set column
           mutate(Gene = ifelse(Gene == "Gapdh", Gene, str_c(Gene, "-", Primer_set))) %>%
           dplyr::select(-Primer_set) %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2")) %>%
           #drop nan and outliers
           filter(!is.nan(Cq)) %>%
           #filter(!(Cq > 32)) %>%
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
        select(Bio_Rep,Treatment, Time, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging


dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Tech_Rep_1, Tech_Rep_2,
                  Treatment,Time, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Bio_Rep",  "Treatment", "Time")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Gapdh_Avg_Cq` = 2^(-Gapdh_Avg_Cq)) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

# map this to make it more conscise
means = dCq_dat %>%
        group_by(Gene, Treatment) %>%
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
          full_join(means, by = c("Gene",  "Treatment")) 
          

write.csv(as.matrix(dCq_dat), path(output_dir, "20240423_snoRNA_bmdc_dCqTable.csv"))

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


dat2 = dCq_dat %>% 
      select(Gene, Treatment, Time, rep_means_2_dCq, rep_se_2_dCq) %>%
      mutate(Treatment = factor(Treatment, levels = c("Control", "L100", "L10_S10_G20")),
             #recreate primer set column
             Primer_set = gsub(".*-", "", Gene),
             Gene = gsub("-.*", "", Gene)) %>%
      unique()
     

p<-ggplot(dat2, aes(x = Treatment, y = rep_means_2_dCq)) +
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  geom_errorbar(mapping =aes(ymin= rep_means_2_dCq - rep_se_2_dCq,
                             ymax= rep_means_2_dCq + rep_se_2_dCq), 
                width=0.4, position = position_dodge(0.5), 
                stat = "identity", 
                colour="black", alpha=0.9, size=0.5) +
  ggh4x::facet_grid2(Primer_set~Gene, independent = "y",scales = "free") +
  scale_fill_manual(values = c("grey", "#8C9FCA", "#F68064")) +
  labs(y = "2^dCq (Target - Gapdh)") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = rel(5)),
        strip.text = element_text(size = rel(1)),
        legend.position = "bottom")
  
save_plot(path(output_dir, "20240423_snora_bmdc_test1.pdf"), p, base_aspect_ratio = 1.5)

