#read in formatted csv file and perform delta delta cq analysis then make plots
library(data.table)
library(fs)
library(tidyverse)
library(ggpubr) #ggarrange
library(scales) #log2 plotting
library(cowplot) #save_plot

#set directories
wd = getwd()
input_dir = path(wd, "input/data")

#read in data
qPCR_dat = read.csv(fs::path(input_dir, "kj_ctu_kdko_exp2_lpsStim_primer2_formatted_editedGapdh.csv")) %>%
           data.frame() %>%
           select(Well, Gene, Target, Oligo, Bio_Rep, Cq, Repeat_values) %>%
           #add columns for qpcr gene, biol rep, tech rep
           mutate(Tech_Rep = ifelse(as.numeric(gsub("[A-P]","", Well)) %% 2 != 0, "1", "2"),
                  #str combine target & oligo cols
                  Target_oligo = gsub("_$", "",str_c(Target, "_", Oligo)),
                  Cq = ifelse(!(is.na(Repeat_values)), Repeat_values, Cq)) %>%
           select(-Repeat_values)
            
dat = qPCR_dat %>%  
           dplyr::select(-Well) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, values_from = Cq, names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rowwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T),
                  `2_Avg_Cq` = 2^(-Avg_Cq))  #%>%# average Cq over technical replicates

Gapdh = dat %>%
        filter(Gene == "Gapdh") %>%
        select(Bio_Rep, Target_oligo, Avg_Cq) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging

dCq_dat =  dat %>%
           select(Gene, Bio_Rep, Tech_Rep_1, Tech_Rep_2, Target_oligo, Avg_Cq, `2_Avg_Cq`) %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Bio_Rep", "Target_oligo")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq,
                  `2_dCq` = 2^(-dCq),
                  `2_Gapdh_Avg_Cq` = 2^(-Gapdh_Avg_Cq),
                  system = ifelse(grepl("g", Target_oligo), "gRNA", "shRNA")) %>% # compute delta Cq (gene Cq - gapdh Cq)
           arrange(Gene)

means = dCq_dat %>%
        group_by(Gene, Target_oligo) %>%
        summarise(rep_means_Cq = mean(Avg_Cq),
                  rep_means_2_Cq = 2^-(rep_means_Cq),
                  rep_means_dCq = mean(dCq),
                  rep_means_2_dCq = 2^-(rep_means_dCq),
                  rep_means_Gapdh = mean(Gapdh_Avg_Cq),
                  rep_means_2_Gapdh = 2^-(rep_means_Gapdh))

dCq_dat = dCq_dat %>%
          full_join(means, by = c("Gene", "Target_oligo")) %>%
          mutate(Target_oligo = factor(Target_oligo, levels = 
                                         c("None","Scr_g1", "Scr_g2", "Myd88_g1",
                                           "Empty_sh1", "Scr_sh2", "Scr_sh3",
                                           "Ctu1_g1", "Ctu1_g2", "Ctu1_g4", "Ctu1_g5",
                                           "Ctu1_sh1", "Ctu1_sh2", "Ctu1_sh3", "Ctu1_sh4",
                                           "Ctu2_g1", "Ctu2_g3", "Ctu2_g4", "Ctu2_g5",
                                           "Ctu2_sh1", "Ctu2_sh2", "Ctu2_sh3", "Ctu2_sh4")),
                 Target = ifelse(grepl("Ctu1", Target_oligo),
                                 "Ctu1",
                                 ifelse(grepl('Ctu2', Target_oligo),
                                        "Ctu2", "Control")))

no_target = dCq_dat %>%
  filter(Target_oligo == "None") %>%
  mutate(system = "gRNA")

dCq_dat = dCq_dat %>%
          rbind(no_target)

write.csv(as.matrix(dCq_dat), path(wd, "output/2023/2023_02_23_dCqTable_editedGapdh.csv"))

#make plots ----
make_Cq_plot = function(gene, dCq = F){
  
  #make relevant plot and add mean line
  if(gene == "Gapdh"){
    p = ggplot(dCq_dat, 
               aes(y = `2_Gapdh_Avg_Cq`, x = Target_oligo)) +
      geom_point(aes(color = Target), 
                 position = position_dodge(width = 0.5)) +
      geom_point(aes(y= rep_means_2_Gapdh), shape=95, size=5, position=position_dodge(width=0.5))
  } else if(gene != "Gapdh" & dCq == F){
    p = ggplot(dCq_dat %>% filter(Gene == gene), 
               aes(y = `2_Avg_Cq`, x = Target_oligo)) +
        geom_point(aes(color = Target), 
                 position = position_dodge(width = 0.5)) +
        geom_point(aes(y= rep_means_2_Cq), shape=95, size=5, position=position_dodge(width=0.5))
  } else if(gene != "Gapdh" & dCq == T) {
    p = ggplot(dCq_dat %>% filter(Gene == gene), 
               aes(y = `2_dCq`, x = Target_oligo)) +
      geom_point(aes(color = Target), 
                 position = position_dodge(width = 0.5)) +
      geom_point(aes(y= rep_means_2_dCq), shape=95, size=5, position=position_dodge(width=0.5))
  }
    
      
  p = p +
  facet_wrap(~system, scales = "free") +
  labs(x = "", 
       y = ifelse(dCq == F, expression(2^-Cq), expression(2^-dCq))) +
  scale_y_log10() +
  #scale_y_continuous(trans = log2_trans())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(2)),
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
    
make_Cq_plot("Ctu1")

p <- ggarrange(plotlist=list(
  ggarrange(plotlist = list(make_Cq_plot("Ctu1"),
                           make_Cq_plot("Gapdh")),
            ncol =1, 
            align = "hv"
            ),
            make_Cq_plot("Ctu1", dCq=T)),
            ncol = 1, 
            labels = NULL,
            common.legend = T,
            legend = "bottom",
            heights = c(0.6,0.4),
            align= "hv") 

save_me <- annotate_figure(p, top = text_grob("Ctu1 Expression (Edited)"))

save_plot(path(wd, "output/2023/ctu2_lps_treated.pdf"), save_me)

# ddCq calcs  ----           
controls = ddCq_dat %>% 
           filter(Treatment == "Control") %>%
           group_by(Gene) %>%
           summarise(Avg_Control_dCq = mean(dCq)) # compute average delta Cq value for controls from each target

ddCq_dat = dCq_dat %>%
           filter(Treatment != "Control") %>% #remove controls before merging
           merge(controls, by = c("Gene")) %>%
           mutate(ddCq = dCq - Avg_Control_dCq,  # compute ddCq (sample Cq - avg dCq for control)
                  expression = 2^-ddCq) %>% # compute expression (2^ddCq)
           dplyr::select(-Sample) %>%
           pivot_wider(id_cols = c("Treatment", "Dose", "Time", "Gene"),names_from = Biol_Rep, 
                       values_from = expression, names_prefix = "Biol_Rep_") %>%
           rowwise() %>%
           mutate(Avg_expression = mean(c(Biol_Rep_1, Biol_Rep_2), na.rm=T),
                  se_expression = sd(c(Biol_Rep_1, Biol_Rep_2)) /sqrt(2)) # average bio reps and compute error (but keep other cols)

dat = ddCq_dat %>%
      mutate(Type = ifelse(nchar(Treatment) ==2, "Pair", "Single"),
             Treatment = factor(Treatment, levels = c("L", "S", "LS")),
             Dose = factor(Dose, levels = c("low", "med", "high")),
             Time = str_c(Time, " hours"),
             Time = factor(Time, levels = c("2 hours", "6 hours")))

#plot by dose
ggplot(dat, aes(x = Treatment, y = Avg_expression, alpha = Dose, fill = Type)) +
       geom_bar(position= position_dodge(0.9), stat = "identity") +
       geom_errorbar(dat, mapping = aes(ymin= Avg_expression-se_expression,
                                   ymax= Avg_expression+se_expression, fill = Dose), width=0.4, position = position_dodge(0.9), 
                stat = "identity",
                colour="black", alpha=0.9, size=0.5) +
       scale_color_manual(values = c( Single = "#707FA1",
                                      Pair = "#3E7563")) +
       scale_fill_manual(values = c( Single = "#707FA1",
                                      Pair = "#3E7563")) +
       scale_alpha_manual(values = c(low = 0.2, med = 0.5, high = 1)) +
      
       ggh4x::facet_grid2(Time ~ Gene, scales = "free", independent = "y") +
       labs(y = "Relative Expression (2^ddCq)", x = NULL) +
       theme(legend.position = "bottom",
            strip.text.y = element_text(angle = 360),
            legend.text= element_text(size = rel(0.5)),
            text = element_text(size = 12))
       

#plot by time

#plot by treatment