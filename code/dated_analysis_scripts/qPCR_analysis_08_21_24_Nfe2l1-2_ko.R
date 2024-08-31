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
input_dir = path(wd, "input/data/2024/2024_08")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "20240821_nfe2l1-2ko_qPCR.csv")) %>%
           data.frame() %>%
           select(Well, Sample, Gene, Bio_Rep, Cq) %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2")) %>%
           #drop unreliable data
           filter(!(is.na(Cq))) 
            
dat = qPCR_dat %>%  
           dplyr::select(-c(Well)) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, values_from = Cq, names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Gapdh = dat %>%
        filter(Gene == "Gapdh") %>%
        select(Bio_Rep, Sample, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging

dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Sample, Tech_Rep_1, Tech_Rep_2, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, 
                 by = c("Bio_Rep", "Sample")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Gapdh_Avg_Cq` = 2^(-Gapdh_Avg_Cq)) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

means = dCq_dat %>%
        group_by(Gene, Sample) %>%
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
          full_join(means, by = c("Gene", "Sample")) 


write.csv(as.matrix(dCq_dat), path(wd, "output/2024/2024_08/20240821/20240821_Nfe2l1-2ko_qPCR_dat.csv"))

plot_me <- dCq_dat %>%
           mutate(Target = gsub("-.*", "", Sample)) %>%
           #filter(!((Sample %in% c("sgCtrl-3", "sgCtrl-2","sgNfe2l2-2", "sgNfe2l1-1", "sgNfe2l1-2")))) %>%
           select(Gene, Sample, Target, rep_means_2_dCq, rep_se_2_dCq) %>%
           #remove individual replicate rows
           unique() %>%
           mutate(Gene = factor(Gene, levels = c("Slc7a11", "Il6", "Il1b")),
                  Target = factor(Target, levels =c("sgCtrl", "sgNfe2l1", "sgNfe2l2"))) %>%
           arrange(Target, Sample)

plot_me2 <- plot_me %>%
            mutate(Sample = factor(Sample, levels = unique(plot_me$Sample))) #%>%
            filter(Target != "sgNfe2l1")
    
#PLOTS
p<-ggplot(plot_me2, aes(x = Sample, 
                y = rep_means_2_dCq)) +
  geom_bar(stat = "identity", 
           color = "black",
           aes(fill = Target)) +
  geom_errorbar(mapping =aes(ymin= rep_means_2_dCq - rep_se_2_dCq,
                             ymax= rep_means_2_dCq + rep_se_2_dCq), 
                width=0.4, position = position_dodge(0.5), 
                stat = "identity", 
                colour="black", alpha=0.9, size=0.5) +
  facet_wrap(~Gene, scales = "free_y", nrow=1) +
  scale_fill_manual(values =c("sgCtrl" ="white",
                              "sgNfe2l1" = "#7CAE00", 
                              "sgNfe2l2" = "#00BFC4")) +
  labs(y = "Relative Expression \n (TARGET - Gapdh)", x = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 8),
        text = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, vjust = 0.5)) 


save_plot(path(wd, "output/2024/2024_08/20240821/20240821_Nfe2l1ko_Nfe2l2ko.pdf"), p,
          base_aspect_ratio = 3.8)



##OLD __________________________________________ 
# ddCq calcs  ----           
controls = il6_plot %>% 
           filter(Treatment == "None") %>%
           group_by(Gene, Virus_amount, Vector, Target_oligo) %>%
           summarise(Avg_Control_dCq = mean(dCq)) # compute average delta Cq value for controls from each target

ddCq_dat = il6_plot %>%
           filter(Treatment != "None") %>% #remove controls before merging
           merge(controls, by = c("Gene", "Vector", "Virus_amount", "Target_oligo")) %>%
           mutate(ddCq = dCq - Avg_Control_dCq,  # compute ddCq (sample Cq - avg dCq for control)
                  expression = 2^-ddCq) %>% # compute expression (2^ddCq)
           #dplyr::select(-Sample) %>%
           pivot_wider(id_cols = c("Treatment", "Target_oligo", "Vector", "Virus_amount"),
                       names_from = Bio_Rep, 
                       values_from = expression, 
                       names_prefix = "Bio_Rep_") %>%
           rowwise() %>%
          filter(!(grepl("Untrans", Target_oligo))) %>%
          mutate(Avg_expression = mean(c(Bio_Rep_1, Bio_Rep_2), na.rm=T),
                 se_expression = sd(c(Bio_Rep_1, Bio_Rep_2)) /sqrt(2)) # average bio reps and compute error (but keep other cols)

#plot ddCq
ggplot(ddCq_dat, aes(x = Target_oligo, y = Avg_expression, fill = Vector)) +
       geom_bar(position= position_dodge(0.9), stat = "identity") +
       geom_errorbar(ddCq_dat, mapping = aes(ymin= Avg_expression-se_expression,
                                        ymax= Avg_expression+se_expression), 
                     width=0.4,
                     position = position_dodge(0.9), 
                stat = "identity",
                colour="black", alpha=0.9, size=0.5) +
       #scale_color_manual(values = c("#F68064")) +
       #scale_fill_manual(values = c("#F68064")) + 
       facet_wrap(~Virus_amount, ncol = 5) +
       labs(y = expression(2^-ddCq), x = NULL) +
       theme(legend.position = "bottom",
            axis.text.x = element_text(angle =45, h = 1, size = rel(3)),
            strip.text.y = element_text(angle = 360),
            legend.text= element_text(size = rel(2)),
            text = element_text(size = rel(5)))
       
save_plot(path(wd, "output/2023/2023_03_22/il6_myd88_viral_titration_relativeExpression_ddCq.pdf"), p)
