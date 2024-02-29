### Title: qPCR Plotting Fxns
### Author: Kate Johnson
### Date: 12/10/22
### Add'tl Notes: Project specific so far - bar plots for combos related qPCR data

# Plot settings ----
#compress PDFs
tools::compactPDF('mypdfs/',  gs_quality='screen')
###PLOTTING THEME ###
theme_update(text = element_text(size = rel(2)),
             legend.key.size = unit(.1, 'inches'),
             legend.text = element_text(size = rel(2.5)),
             panel.background = element_rect(fill = "white", colour = "black"),
             axis.line = element_blank(),
             axis.text.x = element_text(size=rel(3)),
             axis.text.y = element_text(size=rel(3)),
             axis.title.y = element_text(margin=margin(t = 0, r = 1, b = 0, l = 1), size=rel(3)),
             axis.title.x = element_text(size=rel(3)),
             plot.title = element_text(size = rel(2.5), face= "bold", hjust=0.5),
             plot.margin = unit(c(0.1, 0.1, 0.1, 0.2), "inches"),
             strip.text.x = element_text(size = rel(4), face = "bold", margin = margin(0.05,0.05,0.2,0.05, "cm")),
             strip.text.y = element_text(size = rel(4), face = "bold", margin = margin(0.5,0.05,0.05,0.2, "cm")),
             strip.background = element_blank()) #trbl


# Load libraries ----
library(tidyverse)
library(cowplot)

reference_colors =  c(Control = "grey72", L = "#B9CEE6",S = "#8C9FCA", G = "#707FA1", 
                      LS = "#68C3A5" , LG = "#539C84" , SG = "#3E7563",
                      LSG = "#F68064" )

#fill -- how to color bars
#alpha -- how to shade bars
#grouping -- which samples are grouped together (ie dose, time, treatment; 
             #used for placing error bars)
#facet_x -- which variable to place per row 
#facet_y -- which variable to place per column
  
make_qPCR_plot = function(dat, x, y, se, fill, alpha, grouping, facet_x, facet_y, 
                          colors = c(reference_colors["Control"],
                                     reference_colors["G"], 
                                     reference_colors["LG"]), 
                          title = NULL,
                          save_dir=NULL, 
                          save_name = NULL,
                          save=F) {
                 
                 #enquote arguments for mapping
                 x = enquo(x) #they get unquoted for piping into aes with !! later on
                 y = enquo(y)
                 se = enquo(se)
                 alpha = enquo(alpha)
                 fill = enquo(fill)
                 grouping = enquo(grouping)
                 facet_x = enquo(facet_x)
                 facet_y = enquo(facet_y)
                 
                 #make plot!
                 p=ggplot(dat, aes(x = !!x, y = !!y, alpha = !!alpha, fill = !!fill)) +
                                  geom_bar(color = "black", 
                                           position= position_dodge(0.9), stat = "identity") +
                                  geom_errorbar(dat, mapping =aes(ymin= !!y - !!se,
                                                                  ymax= !!y + !!se,
                                                                  fill = !!grouping), 
                                                width=0.4, position = position_dodge(0.9), 
                                                stat = "identity", 
                                                colour="black", alpha=0.9, size=0.5) +
                                  scale_fill_manual(values = c(None = colors[[1]],
                                                               Single = colors[[2]],
                                                               Pair = colors[[3]])) +
                                  scale_alpha_manual(values = c(None = 0, Low = 0.2, Med = 0.5, High = 1)) +
                                  
                                  ggh4x::facet_grid2(c(facet_x, facet_y), scales = "free", independent = "y") +
                                  labs(y = expression(paste("Relative Expression (", 2^{~Delta~Delta~Cq}, ")")), x = NULL, title = title) +
                                  theme(legend.position = "bottom",
                                        text = element_text(size = rel(3)),
                                        strip.text.y = element_text(angle = 0, face = "plain"))
                 
                 #save plot if needed
                 if(save ==T) {
                 save_plot(fs::path(save_dir, str_c(save_name, "qPCR_plot.pdf")), p, dpi=90)
                 }
                 
                 return(p)
          }
