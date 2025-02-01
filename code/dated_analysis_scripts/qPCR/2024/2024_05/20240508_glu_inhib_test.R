#libraries
library(tidyverse)
library(fs)
library(cowplot)


# directories ----
wd = getwd()
input_dir = fs::path(wd, "input/data")
output_dir = fs::path(wd, "output/2024/2024_05")


#read in data
dat = read_csv(path(input_dir, "2024/2024_05/glu_inhib_test.csv")) #%>%

dat2 = dat %>%
      dplyr::select(-Sample) %>%
      mutate(Treatment = factor(Treatment, levels = c("OVA", "L", "LSG")),
             Duration = factor(Duration, levels = c("None", "Stim" ,"Co-Culture", "Stim+Co-culture"))) %>%
      arrange(Duration, Inhihbitor, Treatment, Rep) %>%
      mutate(Treatment_rep = str_c(Treatment, "_", Rep)) %>%
      dplyr::select(-c(Treatment, Rep)) %>%
      pivot_wider(names_from = "Treatment_rep", values_from = "% Proliferation") %>%
      dplyr::select(Duration, Inhihbitor, starts_with("OVA"), starts_with("L_"), starts_with("LS"))

write.csv(as.matrix(dat2), path(output_dir, "glu_inhib_percent_proliferation_formatted.csv"))
