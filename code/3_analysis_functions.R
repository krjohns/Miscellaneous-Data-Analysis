#read in formatted csv file and perform delta delta cq analysis then make plots
library(fs)
library(tidyverse)

#assumes data input with folling columns:
  # Well : well id from plate
  # Biol_Rep : biological replicate info
  # Tech_Rep : technival replicate info
  # Gene : target info (assumes Gapdh is housekeeping control with mouse gene spelling)
  # Cq : Cq values for each well (numeric)
  
# -- project spefific assumptions you might need to change
  # Treatment: adjuvant treatment ("L", "S", "LS", or "Control")
  # Time: timepoint collected (datatype doesnt matter for this but might matter for plotting)
  # Dose: what dose of adjuvant was given ("low", "med", "high")


# note you'll probably need to re-factor your data after running this before plotting!

analyze_qPCR_dat = function(qPCR_dat, save_dir = NULL, save_name = NULL, save =F){
      
    qPCR_dat = qPCR_dat %>%
                data.frame() %>%
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
      merge(Gapdh, by = c("Biol_Rep", "Treatment", "Time", "Dose")) %>% # merge with Gapdh data
      mutate(dCq = Avg_Cq - Gapdh_Avg_Cq)  # compute delta Cq (gene Cq - gapdh Cq)
    
    
    controls = dCq_dat %>% 
      filter(Treatment == "Control") %>%
      group_by(Gene) %>%
      summarise(Avg_Control_dCq = mean(dCq)) # compute average delta Cq value for controls from each target
    
    ddCq_dat = dCq_dat %>%
      filter(Treatment != "Control") %>% #remove controls before merging
      merge(controls, by = c("Gene")) %>%
      mutate(ddCq = dCq - Avg_Control_dCq,  # compute ddCq (sample Cq - avg dCq for control)
             expression = 2^-ddCq) %>% # compute expression (2^ddCq)
      pivot_wider(id_cols = c("Treatment", "Dose", "Time", "Gene"),names_from = Biol_Rep, 
                  values_from = expression, names_prefix = "Biol_Rep_") %>%
      rowwise() %>%
      mutate(Avg_expression = mean(c(Biol_Rep_1, Biol_Rep_2), na.rm=T),
             se_expression = sd(c(Biol_Rep_1, Biol_Rep_2)) /sqrt(2)) # average bio reps and compute error (but keep other cols)
    
    #make the output folder structure if it doesnt already exist
    if(!(dir.exists(save_dir))){
      dir.create(save_dir)
    }
    
    sub_folder = fs::path(save_dir, save_name)
    
    if(!(dir.exists(sub_folder))){
      dir.create(fs::path(save_dir, save_name))
    }
    
    #save csv files if desired
    if(save==T){
    write.csv(as.matrix(ddCq_dat), fs::path(sub_folder, str_c ("Expression_Table.csv")))
    write.csv(as.matrix(dCq_dat), fs::path(sub_folder, str_c("dCq_Table.csv")))
    write.csv(as.matrix(Gapdh), fs::path(sub_folder, str_c ("Gapdh_Table.csv")))
    write.csv(as.matrix(qPCR_dat), fs::path(sub_folder, str_c("Raw_data_formatted.csv")))
    }
    
    return(list(ddCq_dat = ddCq_dat, dCq_dat= dCq_dat, Gapdh_dat = Gapdh, Raw_dat = qPCR_dat))
    
    } 
