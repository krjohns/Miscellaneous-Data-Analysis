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
             plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "inches"),
             strip.text.x = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.text.y = element_text(size = rel(2), face = "bold", margin = margin(0.05,0.05,0.05,0.05, "cm")),
             strip.background = element_rect(size=0.35, color="black")) #trbl

#set directories
wd = getwd()
input_dir = path(wd, "input/data")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "kj_myd88_titration_il6_032223_formatted.csv")) %>%
           data.frame() %>%
           select(Well, Gene, Target, Virus_origin, Virus_amount, Vector, Oligo, Bio_Rep, Cq, Treatment) %>%
           filter(!(grepl("std", Target))) %>%
           #filter(Treatment == "Unstim") %>%
           #select(-Treatment) %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2"),
                  #str combine target & oligo cols
                  Target_oligo = gsub("_$", "",str_c(Target, "_", Oligo)),
                  Vector = ifelse(Virus_amount != "0",str_c(Vector, " (", Virus_origin, ")"), "")) %>%
           select(-Virus_origin)
                  
            
dat = qPCR_dat %>%  
           dplyr::select(-Well) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, values_from = Cq, names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Gapdh = dat %>%
        filter(Gene == "Gapdh") %>%
        select(Bio_Rep, Target_oligo, Vector, Virus_amount, Avg_Cq, Treatment) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging

dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Vector, Virus_amount, Tech_Rep_1, Tech_Rep_2,
                  Treatment, Target_oligo, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Bio_Rep", "Target_oligo", "Treatment",
                               "Vector", "Virus_amount")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Gapdh_Avg_Cq` = 2^(-Gapdh_Avg_Cq)) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

means = dCq_dat %>%
        group_by(Gene, Target_oligo, Treatment, Vector, Virus_amount) %>%
        summarise(rep_means_Cq = mean(Avg_Cq),
                  rep_means_2_Cq = 2^-(rep_means_Cq),
                  rep_means_dCq = mean(dCq),
                  rep_means_2_dCq = 2^-(rep_means_dCq),
                  rep_means_Gapdh = mean(Gapdh_Avg_Cq),
                  rep_means_2_Gapdh = 2^-(rep_means_Gapdh))

dCq_dat = dCq_dat %>%
          full_join(means, by = c("Gene", "Target_oligo", "Treatment", "Vector", "Virus_amount")) %>%
          mutate(Target_oligo = factor(Target_oligo, levels = 
                                         c("Untrans. (-puro)", "Untrans. (+puro)", "Scr_g11", 
                                           "Scr_1-1 (g11)", "Scr_g12",
                                           "Myd88_g14", "Myd88_1-2 (g14)",
                                           "Myd88_1-3")))

dCq_dat = dCq_dat %>%
 # filter(Target_oligo == "None") %>%
  mutate(system = "gRNA")

write.csv(as.matrix(dCq_dat), path(wd, "output/2023/2023_03_22/2023_03_22_myd88_virus_titration.csv"))

#make plots ----
make_Cq_plot = function(gene, dat, colors = Target, dCq = F){
  
  colors = enquo(colors)
  
  #make relevant plot and add mean line
  if(gene == "Gapdh"){
    p = ggplot(dat, 
               aes(y = `2_Gapdh_Avg_Cq`, x = Target_oligo)) +
      geom_point(aes(fill = !!colors),
                 pch = 21,
                 size = 2,
                 color = "black",
                 position = position_dodge(width = 0.5)) +
      geom_point(aes(y= rep_means_2_Gapdh), shape=95, size=5, position=position_dodge(width=0.5))
  } else if(gene != "Gapdh" & dCq == F){
    p = ggplot(dat %>% filter(Gene == gene), 
               aes(y = `2_Avg_Cq`, x = Target_oligo)) +
        geom_point(aes(fill = !!colors),
                pch = 21,
                size = 2,
                color = "black",
                position = position_dodge(width = 0.5)) +
        geom_point(aes(y= rep_means_2_Cq), shape=95, size=5, position=position_dodge(width=0.5))
  } else if(gene != "Gapdh" & dCq == T) {
    p = ggplot(dat %>% filter(Gene == gene), 
               aes(y = `2_dCq`, x = Target_oligo)) +
      geom_point(aes(fill = !!colors),
                 pch = 21,
                 size = 2,
                 color = "black",
                 position = position_dodge(width = 0.5)) +
      geom_point(aes(y= rep_means_2_dCq), shape=95, size=5, position=position_dodge(width=0.5))
  }
  p = p +
  facet_wrap(~system, scales = "free") +
  labs(x = "", 
       y = ifelse(dCq == F, expression(2^-Cq), expression(2^-dCq))) +
  scale_y_log10() +
  #scale_y_continuous(trans = log2_trans())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(3)),
        axis.text.y = element_text(size = rel(2)),
        plot.title = element_text(size = rel(3.5)),
        text = element_text(size = rel(3.5))) 
  
  if(gene == "Gapdh" | dCq == T){
    p = p +
      theme(strip.text.x= element_blank())
  } 
  
  if (dCq == F) {
    p = p +
      theme(axis.text.x = element_blank(),
            legend.position = "none")
  }
  
  return(p)
  
  }


# il6 plot ----
il6_plot <- dCq_dat %>%
            filter(grepl("Myd|Scr|Unt", Target_oligo)) %>%
            mutate(Treatment = ifelse(Treatment == "Pam", "PAM", Treatment),
                  #Sample = ifelse(Treatment == "None", as.character(Target_oligo),
                                 # str_c(Target_oligo, " (+ ", Treatment, ")")),
                  #Sample = factor(Sample, 
                                 # levels = c("Untrans. (-puro)", "Untrans. (-puro) (+ LPS)","Untrans. (+puro)", "Untrans. (+puro) (+LPS)",
                                          #   "Scr_g11", 
                                          #   "Scr_1-1 (g11)", "Scr_g12",
                                           #  "Myd88_g14", "Myd88_1-2 (g14)",
                                           #  "Myd88_1-3")),
                 Target_oligo = ifelse(Treatment == "None", str_c(as.character(Target_oligo), ),
                                              str_c(Target_oligo, " (+ ", Treatment, ")")),
                  Treatment = factor(Treatment, 
                                     levels = c("None", "LPS")))
          

p<-ggarrange(plotlist=list(
  ggarrange(plotlist = list(make_Cq_plot("Il6", il6_plot,
                                         colors = Treatment) + 
                              scale_fill_manual(values =c("white", "#F68064")) +
                              facet_wrap(~Virus_amount, nrow =1),
                           make_Cq_plot("Gapdh", il6_plot, colors = Treatment) +
                           scale_fill_manual(values =c("white", "#F68064")) +
                          facet_wrap(~Virus_amount, nrow =1)),
            ncol =1, 
            align = "hv" 
            ),
            make_Cq_plot("Il6", 
                         il6_plot, 
                         colors = Treatment,
                         dCq=T) +
                      scale_fill_manual(values = c("white", "#F68064")) +
                      facet_wrap(~Virus_amount, nrow=1)),
            ncol = 1, 
            labels = NULL,
            common.legend = T,
            legend = "bottom",
            heights = c(0.5,0.5),
            #widths = c(50,1),
            align= "hv")
  
   #8C9FCA" #B9CEE6"

#'None' = "white",
#"PAM" = "#8C9FCA",
#"LPS_100" = "#F68064",
#"LPS_10" = "#F9B2A2",
#"LPS_1" = "#FDE5E0"

save_me <- annotate_figure(p, top = text_grob("Il6 Expression"))

save_plot(path(wd, "output/2023/2023_03_22/il6_myd88_viral_titration.pdf"), save_me)

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
p <- ggplot(ddCq_dat, aes(x = Target_oligo, y = Avg_expression, fill = Vector)) +
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
