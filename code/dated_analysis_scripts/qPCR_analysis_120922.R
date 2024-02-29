#read in formatted csv file and perform delta delta cq analysis then make plots
library(data.table)
library(dplyr)
library(fs)
library(tidyverse)
library(ggplot2)

#set directories
wd = "/Users/Kate/Documents/Chevrier_Lab/Projects/U01/Data/qPCR_data/2022_1208_BMDC_LS_cytokines"

#read in data
qPCR_dat = read.csv(fs::path(wd, "BMDC_qPCR_formatted_120822.csv")) %>%
           data.frame() %>%
           mutate(Sample = str_c(Treatment, "_", Dose, "_", Time)) %>% #record sample info
           dplyr::select(-Well) %>%  #drop well info for pivoting
           pivot_wider(names_from = Tech_Rep, values_from = Cq, names_prefix = "Tech_Rep_") %>% #pivot to have technical reps in different cols
           rowwise() %>% #prepare to compute rpwwise average 
           mutate(Avg_Cq = mean(c(Tech_Rep_1, Tech_Rep_2), na.rm=T)) # average Cq over technical replicates

Gapdh = qPCR_dat %>%
        filter(Gene == "Gapdh") %>%
        dplyr::select(-c(Gene, contains('Tech_Rep'))) %>%
        rename_at(vars(Avg_Cq), ~str_c("Gapdh_", .))#record Gapdh genes for merging

dCq_dat = qPCR_dat %>%
           filter(Gene != "Gapdh") %>% # remove Gapdh rows
           merge(Gapdh, by = c("Biol_Rep", "Treatment", "Time", "Dose", "Sample")) %>% # merge with Gapdh data
           mutate(dCq = Avg_Cq - Gapdh_Avg_Cq)  # compute delta Cq (gene Cq - gapdh Cq)
           
             
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