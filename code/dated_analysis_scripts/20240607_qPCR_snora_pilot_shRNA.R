#read in formatted csv file and perform delta delta cq analysis then make plots
library(data.table)
library(fs)
library(tidyverse)
library(ggpubr) #ggarrange
library(scales) #log2 plotting
library(cowplot) #save_plot
library(monochromeR)
 

#set directories
wd = getwd()
input_dir = path(wd, "input/data/2024/2024_06")
output_dir = path(wd, "output/2024/2024_06/20240607")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "20240607_shSnora_pilot_repeat.csv")) %>%
           data.frame() %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2"),
                  #add tech reps for repeat of controls
                  Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) == 21 &
                           gsub("21","", Well) %in% c("A", "C", "E", "G", "I", "K") , "3", 
                           ifelse(as.numeric(gsub("[A-P]","", Well)) == 21 &
                                    gsub("21","", Well) %in% c("B", "D", "F", "H", "J", "L") ,
                                  "4", Tech_Rep))) %>%
           #drop nan and outliers
           #filter(!is.nan(Cq)) %>%
           #filter(!(Cq > 32)) %>%
           #filter(!(Gene == "Actb" & Cq > 30))
            
dat = qPCR_dat %>%  
           dplyr::select(-Well) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, 
                       values_from = Cq, 
                       names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2, Tech_Rep_3, Tech_Rep_4), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Actb = dat %>%
        filter(Gene == "Actb") %>%
        select(Sample, Bio_Rep, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Actb_", .))#record Gapdh genes for merging


dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Tech_Rep_1, Tech_Rep_2,
                  Sample,Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Actb") %>% # remove Gapdh rows
           merge(Actb, by = c("Bio_Rep", "Sample")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Actb_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Actb_Avg_Cq` = 2^(-Actb_Avg_Cq)) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

# map this to make it more conscise
means = dCq_dat %>%
        group_by(Gene, Sample) %>%
        summarise(rep_means_Cq = mean(Avg_Cq),
                  rep_se_Cq = sd(Avg_Cq)/sqrt(2),
                  rep_means_2_Cq = 2^-(rep_means_Cq),
                  rep_means_dCq = mean(dCq),
                  rep_se_2_dCq = sd(2^-dCq)/sqrt(2),
                  rep_means_2_dCq = 2^-(rep_means_dCq),
                  rep_means_Actb = mean(Actb_Avg_Cq),
                  rep_se_Actb = sd(Actb_Avg_Cq)/sqrt(2),
                  rep_means_2_Actb = 2^-(rep_means_Actb))

dCq_dat = dCq_dat %>%
          full_join(means, by = c("Gene",  "Sample")) 
          

write.csv(as.matrix(dCq_dat), path(output_dir, "20240607_shSnora_pilot_dCqTable.csv"))


#plot
dat2 = dCq_dat %>% 
       select(Gene, Sample,  rep_means_2_dCq, rep_se_2_dCq) %>%
       #remove individual replicate rows
       unique() %>%
       #Add columns with additional info
       mutate(Sample = gsub("73ab_pool", "73pool_ab", Sample)) %>%
       separate(Sample, into = c("Target_type", "shRNA"), sep = "-", remove = F) %>%
       mutate(Pool = ifelse(grepl("pool", Sample), "Pooled shRNA", "Single shRNA"),
              shRNA = ifelse(Pool == "Pooled shRNA" , "Pool", shRNA),
              Target = gsub("_.*", "", gsub("sh","",Target_type)),
              Type = ifelse(grepl("OE", Sample), "OE",
                            ifelse(grepl("-", Sample), "Single KD", "Pooled KD"))) 

     
p2<-ggplot(dat2 %>% filter(Gene == "Snora73b" & (Pool == "Pooled shRNA" | Target == "Ctrl")), 
       aes(x = Sample, y = rep_means_2_dCq)) +
  geom_bar(stat = "identity", 
           color = "black",
           aes(fill = Target, alpha = Type)) +
  geom_errorbar(mapping =aes(ymin= rep_means_2_dCq - rep_se_2_dCq,
                             ymax= rep_means_2_dCq + rep_se_2_dCq), 
                width=0.4, position = position_dodge(0.5), 
                stat = "identity", 
                colour="black", alpha=0.9, size=0.5) +
  #facet_wrap(~Pool, scales = "free") +
  scale_fill_manual(values =c("Ctrl" ="white",
                              "Snora57" = "#7CAE00", 
                              "Snora73a" = "#00BFC4",
                              "Snora73b" = "#F8766D",
                              "Snora73pool" = "#C77CFF")) +
  scale_alpha_manual(values = c("OE" =1,
                                "Pooled KD" = 1,
                                "Single KD" = 1)) +
  labs(y = NULL, x = NULL, 
       title = "Pooled shRNA") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, vjust = 0.5))


p1<-ggplot(dat2 %>% filter(Gene == "Snora73b" & Pool == "Single shRNA" & Sample != "shSnora73a-2"), 
       aes(x = Sample, y = rep_means_2_dCq)) +
  geom_bar(stat = "identity", 
           color = "black",
           aes(fill = Target, alpha = Type)) +
  geom_errorbar(mapping =aes(ymin= rep_means_2_dCq - rep_se_2_dCq,
                             ymax= rep_means_2_dCq + rep_se_2_dCq), 
                width=0.4, position = position_dodge(0.5), 
                stat = "identity", 
                colour="black", alpha=0.9, size=0.5) +
  #facet_wrap(~Pool, scales = "free") +
  scale_fill_manual(values =c("Ctrl" ="white",
                              "Snora57" = "#7CAE00", 
                              "Snora73a" = "#00BFC4",
                              "Snora73b" = "#F8766D",
                              "Snora73pool" = "#C77CFF")) +
  scale_alpha_manual(values = c("OE" =1,
                                "Pooled KD" = 1,
                                "Single KD" = 1)) +
  labs(y = "Relative Expression \n (Target - Actb)", x = NULL, 
       title = "Single shRNA") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, vjust = 0.5)) 

  

p<-plot_grid(plotlist = list(p1,p2), nrow = 1, align = "hv", rel_widths =c(1,1))
p

save_plot(path(output_dir, "KD/20240607_snora73b_KD.pdf"), p, base_aspect_ratio = 2.5)

