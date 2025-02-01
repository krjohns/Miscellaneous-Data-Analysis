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
input_dir = path(wd, "input/data/2024/2024_05")
output_dir = path(wd, "output/2024/2024_05/20240509")

#read in data
raw_dat = read_csv(fs::path(input_dir, "20240509_bmdc_snora_formatted.csv")) 
           
#separate out gene list info
gene_list = raw_dat %>%
            select(Column, Gene) %>%
            na.omit() %>%
            mutate(Gene = gsub("_.*", "", Gene))
  
qPCR_dat = raw_dat %>%
           #drop gene / column info before merging back in
           select(-c(Column, Gene)) %>%
           filter(!is.nan(Cq)) %>%
           separate(Well, into = c("Row", "Column"), sep = 1, remove = T) %>%
           #remove trailing zero to match with gene_list numbering format
           #also make numeric for adding target #s for I-P rows (targets # 25-30)
           mutate(Column = as.numeric(gsub("^0", "", Column)),
                  #add + 24 to I-P wells to match target numbering
                  Column = ifelse(grepl("[I-P]", Row), Column + 24, Column),
                  Tech_Rep = ifelse(as.numeric(factor(Row)) %% 2 != 0, "1", "2"),
                  Treatment = ifelse(grepl("[A-D, I-L]", Row), "Control", "LPS (4h)"),
                  Bio_Rep = ifelse(grepl("[A-B, E-F, I-J, M-N]", Row), 1,2),
                  Pass = ifelse(Cq < 30, T, F)) %>%
           merge(gene_list, by = "Column") %>%
           arrange(Gene, Row)
           
            
dat = qPCR_dat %>%  
           dplyr::select(-c(Row, Column, Pass)) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, 
                       values_from = Cq, 
                       names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Gapdh = dat %>%
        filter(Gene == "Gapdh") %>%
        select(Bio_Rep, Treatment, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging


dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Tech_Rep_1, Tech_Rep_2,
                  Treatment, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Bio_Rep",  "Treatment")) %>% # merge with Gapdh data
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
          

write.csv(as.matrix(dCq_dat), path(output_dir, "20240509_snoRNA_bmdc_dCqTable.csv"))

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
      select(Gene, Treatment, rep_means_2_dCq, rep_se_2_dCq) %>%
      unique() %>%
      mutate(human = ifelse(grepl("[()]", Gene), gsub(".*[(]", "", Gene), str_to_upper(Gene)),
             human = ifelse(human == "SNORA7A)", "SNORA7", human)) %>%
      arrange(human)

dat3 = dat2 %>%
       mutate(Gene = factor(Gene, levels = unique(dat2$Gene)))
     

p<-ggplot(dat3, aes(x = Treatment, y = rep_means_2_dCq)) +
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  geom_errorbar(mapping =aes(ymin= rep_means_2_dCq - rep_se_2_dCq,
                             ymax= rep_means_2_dCq + rep_se_2_dCq), 
                width=0.4, position = position_dodge(0.5), 
                stat = "identity", 
                colour="black", alpha=0.9, linewidth=0.5) +
  facet_wrap(~Gene ,scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("grey", "#8C9FCA", "#F68064")) +
  labs(y = "2^dCq (Target - Gapdh)", x = NULL) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.text.y = element_text(size = rel(1.5)),
        text = element_text(size = 4.5),
        strip.text.x = element_text(size = rel(1.5)),
        legend.position = "none")
  
save_plot(path(output_dir, "20240509_snora_bmdc_primerset1.pdf"), p, base_aspect_ratio = 1.9)

