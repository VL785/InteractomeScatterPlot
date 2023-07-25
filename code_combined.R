# install packages --------------------------------------------------------
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("here")
install.packages("fuzzyjoin")
install.packages("gridExtra") # tableGrob
install.packages("ggrepel")
install.packages("grid)")
install.packages("formattable")
install.packages("patchwork")
install.packages("ggforce")
install.packages("magrittr")
# naming and general variables ------------------------------------------------------------------
naming<-"combined_adjusted_x-axis_2"
#qqSIZEshape <- "1"
# load packages ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(here)
library(fuzzyjoin)
library(gridExtra) 
library(ggrepel) 
library(grid)
library(formattable)
library(patchwork)
library(ggforce)
library(magrittr)

# load mitocarta and any other datasets to append -------------------------

mitocarta <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/code/databases/mitocarta/Mouse.MitoCarta3.0.txt")
mitocarta <- mitocarta %>%
  select(Symbol, MitoCarta3.0_MitoPathways, MitoCarta3.0_SubMitoLocalization, UniProt, MitoCarta3.0_List) %>%
  # Make emtpy (0) UniProt entries NA, otherwise fuzzy_join will misbehave
  mutate(UniProt = ifelse(`UniProt` == 0, NA, UniProt))
mitocarta

# kidney upload -----------------------------------------------------------
my.data_tissue_kidney <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/kidney.txt")
colnames(my.data_tissue_kidney)
# kidney control general variables  ------------------------------------------------------
log2<- expression(Log[2])
fold.thres_tissue_kidney <- 0.5 # User-defined fold-change value for biological significance
p.thres_tissue_kidney <- 0.05 # User-defined p-val for statistical significance
p.thres2_tissue_kidney <- 1.31
neg.fold.thres_tissue_kidney <- -fold.thres_tissue_kidney
micos<- c("Apoo","Apool","Immt","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3")
paste(micos,collapse="|")

# kidney clean up data  -----------------------------------------------------------
my.data.clean_tissue_kidney <- my.data_tissue_kidney %>% 
  #left_join(mitocarta, match_fun = str_detect, by = c("Protein IDs" = "UniProt")) %>%
  left_join(mitocarta, by = c("Protein IDs" = "UniProt")) %>%
  #rename(`difcol_tissue_kidney`= contains("-Log2 Fold Change")) %>%
  rename(`difcol_tissue_kidney`= contains("Student's T-test Difference")) %>%
  # rename(`FDR` = contains("Student's T-test q-value")) %>%
  #rename(`Significant` = contains("Student's T-test Significant ")) %>%
  rename(`'-log10 p-value_tissue_kidney` = contains("-Log Student's T-test p-value")) %>%
  #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
  #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
  # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
  mutate(`sig_tissue_kidney` = ifelse(`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` > fold.thres_tissue_kidney) | (`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` < neg.fold.thres_tissue_kidney)), "S", "NS")) %>%
  mutate(`label_tissue_kidney` = case_when((`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` < neg.fold.thres_tissue_kidney) ~ "left"),
                                       (`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` > fold.thres_tissue_kidney) ~ "right"),
                                       (`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) ~ "other"))) %>%
  #Remove excess info from 'Gene names'
  mutate(`hakidneyICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
  mutate(`GenesNew_tissue_kidney` = `Gene names`) %>%
  #mutate(`GenesNew_tissue_kidney` = `Genes`) %>%
  mutate_at("GenesNew_tissue_kidney", list(~ gsub(";.*", " ", .))) %>%
  mutate_at("GenesNew_tissue_kidney", list(~ gsub("-.*", " ", .)))

# kidney run min/max & control visual variables  ----------------------------------------------------
# /* Define axes domains based on max/min values
round_up <- function(x, to = 10) {
  to * (x %/% to + as.logical(x %% to))
}
round_down <- function(x, to = 10) {
  to * (-x %/% to + as.logical(x %% to)) * -1
}

x_max_tissue_kidney <- round_up(max(my.data.clean_tissue_kidney $`difcol_tissue_kidney`, na.rm = TRUE), 1)
y_max_tissue_kidney <- round_up(max(my.data.clean_tissue_kidney$`'-log10 p-value_tissue_kidney`, na.rm = TRUE), 1)
x_min_tissue_kidney <- round_down(min(my.data.clean_tissue_kidney$`difcol_tissue_kidney`, na.rm = TRUE), 1)
y_min_tissue_kidney <- round_up(min(my.data.clean_tissue_kidney$`'-log10 p-value_tissue_kidney`, na.rm = TRUE), 1)
#y_min <- 0
y_max_tissue_kidney
y_min_tissue_kidney
x_min_tissue_kidney
x_max_tissue_kidney
# */

#adjust x limits
#x_limit_min <- #0 
#x_limit_min <- `x_min_tissue_kidney` 
x_limit_min <- 0
#x_limit_max <- `x_max_tissue_kidney`+2 #adjust after the plus
x_limit_max <- `x_max_tissue_kidney`
#adjust y limits
y_limit_min <- `y_min_tissue_kidney` 
#y_limit_max <- `y_max_tissue_kidney`+2 #adjust after the plus
y_limit_max <- `y_max_tissue_kidney`

my_theme_tissue_kidney <- theme_bw() +
  theme(plot.margin=unit(c(2,2,2,2), 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        # axis.text.x = element_text(size = rel(0.2)),
        # #face = "bold"),
        # axis.text.y = element_text(size = rel(0.2)),
        # #face = "bold"),
        # axis.line = element_line(colour = "black"),
        # axis.title.x = element_text(size = rel(1),
        #                             face = "bold"),
        # axis.title.y = element_text(size = rel(1),
        #                             face = "bold"),
        #              # axis.ticks.x = element_line(size = 4),
        #                 #                         angle = -90),
        # axis.ticks.length=unit(0.2, "cm"),
        # legend.justification = c(1, 0), 
        #legend.position=c(-0.8, -0.15),
        #legend.justification = c("right","bottom"),
        legend.text          = element_text(size = 15),
        # legend.box = "vertical",
        # legend.margin=margin(), 
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(colour = "black")
  )

my_layers_tissue_kidney = list(
  labs(title = "Kidney"),
  # x and y axes labels and scales
  xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
  #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
  ylab(expression(-log[10]~"p-value")),
  scale_x_continuous(
    expand = c(0,0),
    limits = c( `x_limit_min`, `x_limit_max`),
    breaks = seq(`x_min_tissue_kidney`, `x_max_tissue_kidney`, by = 1)), #spacing of axis 
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, `y_limit_max`)))

my.data.visual_tissue_kidney <- my.data.clean_tissue_kidney %>%
  mutate(`visual_tissue_kidney` = case_when(
    `Gene names` == "Immt" ~ "Bait",
    `sig_tissue_kidney` == "S" & hakidneyICOS == "MICOS" ~ "MICOS",
    # `sig_tissue_kidney` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
    # `sig_tissue_kidney` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
    `sig_tissue_kidney` == "S" & MitoCarta3.0_List == "MitoCarta3.0" ~ "MitoCarta3.0",
    `sig_tissue_kidney` == "S" & is.na(MitoCarta3.0_List) ~ "Other",
    `sig_tissue_kidney` == "S" ~ "Significant",
    `sig_tissue_kidney` == "NS" ~ "NS"
    
  ))
# 
# my.data.visual_tissue_kidney %>% 
#   filter(sig_tissue_kidney=="S")%>%
#     arrange(desc(difcol_tissue_kidney))
# head(my.data.visual_tissue_kidney)
#   

visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
names(visual.colors) <- levels(visual.levels)
# visual.colors <- c("gray", "black")
# visual.levels <- factor(levels=c("NS", "Significant"))
# names(visual.colors) <- levels(visual.levels)
# metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolikidney"))
# metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
# names(metals.colors) <- levels(metals.levels)

# imputed.levels <- factor(levels=c(0,1,2,3,4))
# imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
# names(imputed.colors) <- levels(imputed.levels)


# kidney run visual   --------------------------------------------------------
# p-value line label - adjust number for position changes 
line_horizontal_position_x <- `x_max_tissue_kidney` - 1
line_horizontal_position_y <- 1.31
line_horizontal_size <- 3 
# enrichment line label placement
line_vertical_positive_position_x <- `fold.thres_tissue_kidney`
line_vertical_positive_position_y <- `y_max_tissue_kidney` - 1 #negative line = same
line_vertical_negative_position_y <- `y_max_tissue_kidney` - 1 
line_vertical_negative_position_x <- `neg.fold.thres_tissue_kidney` 

# enrichment line label to actually label and ad the "x"
hlabel_format_tissue_kidney<- formattable(2^`fold.thres_tissue_kidney`, digits = 0, format = "f")
hlabel_tissue_kidney<- paste("x",`hlabel_format_tissue_kidney`, sep='')
line_vertical_size <- 3
my.plot_tissue_kidney <- my.data.visual_tissue_kidney %>%
  ggplot(mapping = aes(x = `difcol_tissue_kidney`, y= `'-log10 p-value_tissue_kidney`, label = `GenesNew_tissue_kidney`)) + 
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `visual_tissue_kidney` =="Significant"),
    aes(color = visual_tissue_kidney), #black
    size = 1.5) +
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `visual_tissue_kidney` =="NS"),
    aes(color = visual_tissue_kidney),#gray
    size = 1.5) +
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `visual_tissue_kidney` =="MICOS"),
    aes(color = visual_tissue_kidney), #blue
    size = 1.5) +
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `visual_tissue_kidney` =="Bait"),
    aes(color = visual_tissue_kidney), #red
    size = 1.5) +
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `sig_tissue_kidney` == "S" & MitoCarta3.0_List == "MitoCarta3.0" ),
    aes(color = visual_tissue_kidney),
    size = 1.5) +
  geom_point(
    data = subset(my.data.visual_tissue_kidney, `sig_tissue_kidney` == "S" & is.na(MitoCarta3.0_List)),
    aes(color = visual_tissue_kidney),
    size = 1.5) +
  # gal up 
  # geom_text_repel(  ###TO LABEL ALL DATA POINTS
  #   data= subset (my.data.visual_tissue_kidney, `label_tissue_kidney`== "right"),
  #   aes(color = visual_tissue_kidney),
  #   max.overlaps=1000,
  #   size=3,
  geom_text_repel(
    data= subset (my.data.visual_tissue_kidney, `hakidneyICOS`== "MICOS"),
    aes(color = visual_tissue_kidney),
    max.overlaps=1000,
    size=4.25, #qqSIZEshape
    # geom_text_repel(
    #   data= subset (my.data.visual_tissue_kidney,`difcol_tissue_kidney` >1.8),
    #   aes(color = visual_tissue_kidney),
    #   max.overlaps=1000,
    #   size=4,
    #box.padding = unit(2, "lines"),
    point.padding = 0,
    box.padding = 1, #white margin border as a box from corners of plot 
    #position = position_nudge_repel(x = 0, y = 0),
    nudge_x = 0.08,
    #box.padding = 10,
    force = 0.1,  #changes spacing around lines
    # force = 10,  #changes spacing around lines
    #nudge_y = 1,
    #segment.curvature = -0.1,
    #segment.ncp = 3,
    #segment.angle = 90,
    #segment.square= FALSE,
    min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
    seed = 42,
    # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
    # segment.ncp =3,
    # segment.angle = 20,
    segment.size = 0.4, #line thickness
    #point.size = 20
    # Repel just the labels and totally ignore the data points
    # p + geom_text_repel(point.size = NA)
    
    #segment.angle = -180 ## left vs right adjustment 
    #segment.angle = 180 ## left vs right adjustment 
    # Repel away from the left edge, not from the right.
    # xlim = c(NA, Inf),
    # # Do not repel from top or bottom edges.
    # ylim = c(-Inf, Inf)
  )+
  #alter the grey lines
  geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
  annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
  #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
  geom_vline(xintercept=`neg.fold.thres_tissue_kidney`, linetype='dashed', col = 'grey48') +
  annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_kidney`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
  geom_vline(xintercept=`fold.thres_tissue_kidney`, linetype='dashed', col = 'grey48') +
  annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_kidney`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
  scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))
# facet_zoom(xy = example_values_x > 0.8 & example_values_x < 1.8, zoom.data=zoom) +   # Note the zoom.data argument
# geom_label_repel(data = dftxt, aes(label = GenesNew_tissue_kidney))


p1 <- my.plot_tissue_kidney + my_theme_tissue_kidney + my_layers_tissue_kidney
p1

# heart upload   -------------------------------------------
my.data_tissue_heart <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/heart.txt")

colnames(my.data_tissue_heart)


# heart control general variables   ----------------------------------------------
log2<- expression(Log[2])
fold.thres_tissue_heart <- 0.5 # User-defined fold-change value for biological significance
p.thres_tissue_heart <- 0.05 # User-defined p-val for statistical significance
p.thres2_tissue_heart <- 1.31
neg.fold.thres_tissue_heart <- -fold.thres_tissue_heart
#micos<- c("Apoo","Apool","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3")
#paste(micos,collapse="|")

# clean up data heart  --------------------------------------------------
my.data.clean_tissue_heart <- my.data_tissue_heart %>%
 # left_join(mitocarta, match_fun = str_detect, by = c("Protein IDs" = "UniProt")) %>%
  #rename(`difcol_tissue_heart`= contains("-Log2 Fold Change")) %>%
  rename(`difcol_tissue_heart`= contains("Student's T-test Difference")) %>%
  # rename(`FDR` = contains("Student's T-test q-value")) %>%
  #rename(`Significant` = contains("Student's T-test Significant ")) %>%
  rename(`'-log10 p-value_tissue_heart` = contains("-Log Student's T-test p-value")) %>%
  #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
  #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
  # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
  mutate(`sig_tissue_heart` = ifelse(`'-log10 p-value_tissue_heart` > -log10(p.thres_tissue_heart) & (`difcol_tissue_heart` > fold.thres_tissue_heart) | (`'-log10 p-value_tissue_heart` > -log10(p.thres_tissue_heart) & (`difcol_tissue_heart` < neg.fold.thres_tissue_heart)), "S", "NS")) %>%
  mutate(`label_tissue_heart` = case_when((`'-log10 p-value_tissue_heart` > -log10(p.thres_tissue_heart) & (`difcol_tissue_heart` < neg.fold.thres_tissue_heart) ~ "left"),
                                          (`'-log10 p-value_tissue_heart` > -log10(p.thres_tissue_heart) & (`difcol_tissue_heart` > fold.thres_tissue_heart) ~ "right"),
                                          (`'-log10 p-value_tissue_heart` > -log10(p.thres_tissue_heart) ~ "other"))) %>%
  # Remove excess info from 'Gene names'
  mutate(`hakidneyICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
  mutate(`GenesNew_tissue_heart` = `Gene names`) %>%
  # mutate(`GenesNew_tissue_heart` = `Genes`) %>%
  mutate_at("GenesNew_tissue_heart", list(~ gsub(";.*", " ", .))) %>%
  mutate_at("GenesNew_tissue_heart", list(~ gsub("-.*", " ", .)))

# heart run min/max & control visual variables   --------------------------------


# /* Define axes domains based on max/min values
round_up <- function(x, to = 10) {
  to * (x %/% to + as.logical(x %% to))
}
round_down <- function(x, to = 10) {
  to * (-x %/% to + as.logical(x %% to)) * -1
}

# x_max_tissue_heart <- round_up(max(my.data.clean_tissue_heart $`difcol_tissue_heart`, na.rm = TRUE), 1)
x_max_tissue_heart <- 5.5
y_max_tissue_heart <- round_up(max(my.data.clean_tissue_heart$`'-log10 p-value_tissue_heart`, na.rm = TRUE), 1)
x_min_tissue_heart <- round_down(min(my.data.clean_tissue_heart$`difcol_tissue_heart`, na.rm = TRUE), 1)
y_min_tissue_heart <- round_up(min(my.data.clean_tissue_heart$`'-log10 p-value_tissue_heart`, na.rm = TRUE), 1)
#y_min <- 0
y_max_tissue_heart
y_min_tissue_heart
x_min_tissue_heart
x_max_tissue_heart
# */

#adjust x limits
#x_limit_min <- #0 
#x_limit_min <- `x_min_tissue_heart` 
x_limit_min <- 0
#x_limit_max <- `x_max_tissue_heart`+2 #adjust after the plus
x_limit_max <- `x_max_tissue_heart`
#adjust y limits
y_limit_min <- `y_min_tissue_heart` 
#y_limit_max <- `y_max_tissue_heart`+2 #adjust after the plus
y_limit_max <- `y_max_tissue_heart`

my_theme_tissue_heart <- theme_bw() +
  theme(plot.margin=unit(c(2,2,2,2), 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        # axis.text.x = element_text(size = rel(1.5)),
        # #face = "bold"),
        # axis.text.y = element_text(size = rel(1.5)),
        # #face = "bold"),
        # axis.line = element_line(colour = "black"),
        # axis.title.x = element_text(size = rel(1),
        #                             face = "bold"),
        # axis.title.y = element_text(size = rel(1),
        #                             face = "bold"),
        #              # axis.ticks.x = element_line(size = 4),
        #                 #                         angle = -90),
        # axis.ticks.length=unit(1.5, "cm"),
        # legend.justification = c(1, 0), 
        #legend.position=c(-0.8, -0.15),
        #legend.justification = c("right","bottom"),
        legend.text          = element_text(size = 15),
        # legend.box = "vertical",
        # legend.margin=margin(), 
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(colour = "black")
  )

my_layers_tissue_heart = list(
  labs(title = "Heart"),
  # x and y axes labels and scales
  xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
  #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
  ylab(expression(-log[10]~"p-value")),
  scale_x_continuous(
    expand = c(0,0),
    limits = c( `x_limit_min`, `x_limit_max`),
    breaks = seq(`x_min_tissue_heart`, `x_max_tissue_heart`, by = 1)), #spacing of axis 
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, `y_limit_max`)))

my.data.visual_tissue_heart <- my.data.clean_tissue_heart %>%
  mutate(`visual_tissue_heart` = case_when(
    `Gene names` == "Immt" ~ "Bait",
    `sig_tissue_heart` == "S" & hakidneyICOS == "MICOS" ~ "MICOS",
    # `sig_tissue_heart` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
    # `sig_tissue_heart` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
    `sig_tissue_heart` == "S" & MitoCarta3.0_List == "+" ~ "MitoCarta3.0",
    `sig_tissue_heart` == "S" & MitoCarta3.0_List != "+" ~ "Other",
    `sig_tissue_heart` == "S" ~ "Significant",
    `sig_tissue_heart` == "NS" ~ "NS"
    
  ))

visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
names(visual.colors) <- levels(visual.levels)
# visual.colors <- c("gray", "black")
# visual.levels <- factor(levels=c("NS", "Significant"))
# names(visual.colors) <- levels(visual.levels)
# metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolikidney"))
# metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
# names(metals.colors) <- levels(metals.levels)

# imputed.levels <- factor(levels=c(0,1,2,3,4))
# imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
# names(imputed.colors) <- levels(imputed.levels)

# heart run visual   -----------------------------------------------------

# p-value line label - adjust number for position changes 
line_horizontal_position_x <- `x_max_tissue_heart` - 1
line_horizontal_position_y <- 1.31
line_horizontal_size <- 3 
# enrichment line label placement
line_vertical_positive_position_x <- `fold.thres_tissue_heart`
line_vertical_positive_position_y <- `y_max_tissue_heart` - 1 #negative line = same
line_vertical_negative_position_y <- `y_max_tissue_heart` - 1
line_vertical_negative_position_x <- `neg.fold.thres_tissue_heart` 

# enrichment line label to actually label and ad the "x"
hlabel_format_tissue_heart<- formattable(2^`fold.thres_tissue_heart`, digits = 0, format = "f")
hlabel_tissue_heart<- paste("x",`hlabel_format_tissue_heart`, sep='')
line_vertical_size <- 3

 my.plot_tissue_heart<- my.data.visual_tissue_heart %>%
    ggplot(mapping = aes(x = `difcol_tissue_heart`, y= `'-log10 p-value_tissue_heart`, label = `GenesNew_tissue_heart`)) + 
    geom_point(
      data = subset(my.data.visual_tissue_heart, `visual_tissue_heart` =="Significant"),
      aes(color = visual_tissue_heart), #black
      size = 1.5) +
    geom_point(
      data = subset(my.data.visual_tissue_heart, `visual_tissue_heart` =="NS"),
      aes(color = visual_tissue_heart),#gray
      size = 1.5) +
    geom_point(
      data = subset(my.data.visual_tissue_heart, `visual_tissue_heart` =="MICOS"),
      aes(color = visual_tissue_heart), #blue
      size = 1.5) +
    geom_point(
      data = subset(my.data.visual_tissue_heart, `visual_tissue_heart` =="Bait"),
      aes(color = visual_tissue_heart), #red
      size = 1.5) +
    geom_point(
      data = subset(my.data.visual_tissue_heart, `sig_tissue_heart` == "S" & MitoCarta3.0_List == "+" ),
      aes(color = visual_tissue_heart),
      size = 1.5) +
    geom_point(
      data = subset(my.data.visual_tissue_heart, `sig_tissue_heart` == "S" & MitoCarta3.0_List != "+" ),
      aes(color = visual_tissue_heart),
      size = 1.5) +
    # geom_text_repel(  ###TO LABEL ALL DATA POINTS
    #   data= subset (my.data.visual_tissue_heart, `label_tissue_heart`== "right"),
    #   aes(color = visual_tissue_heart),
    #   max.overlaps=1000,
    #   size=3,
    geom_text_repel(
       data= subset (my.data.visual_tissue_heart, `hakidneyICOS`== "MICOS"),
     # data= subset (my.data.visual_tissue_heart, `difcol_tissue_heart`> 1.8 & `sig_tissue_heart`=="S"),
      aes(color = visual_tissue_heart),
      max.overlaps=1000,
      size=4.25, #qqSIZEshape
      #box.padding = unit(2, "lines"),
      point.padding = 0,
      box.padding = 1, #white margin border as a box from corners of plot 
      #position = position_nudge_repel(x = 0, y = 0),
      nudge_x = 0.08,
      #box.padding = 10,
      force = 0.1,  #changes spacing around lines
      # force = 10,  #changes spacing around lines
      #nudge_y = 1,
      #segment.curvature = -0.1,
      #segment.ncp = 3,
      #segment.angle = 90,
      #segment.square= FALSE,
      min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
      seed = 42,
      # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
      # segment.ncp =3,
      # segment.angle = 20,
      segment.size = 0.4, #line thickness
      #point.size = 20
      # Repel just the labels and totally ignore the data points
      # p + geom_text_repel(point.size = NA)
      
      #segment.angle = -180 ## left vs right adjustment 
      #segment.angle = 180 ## left vs right adjustment 
      # # Repel away from the left edge, not from the right.
      # xlim = c(NA, Inf),
      # # Do not repel from top or bottom edges.
      # ylim = c(-Inf, Inf)
      
    )+
    #alter the grey lines
    geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
    annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
    #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
    geom_vline(xintercept=`neg.fold.thres_tissue_heart`, linetype='dashed', col = 'grey48') +
    annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_heart`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
    geom_vline(xintercept=`fold.thres_tissue_heart`, linetype='dashed', col = 'grey48') +
    annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_heart`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
    scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))

 # BAT upload  --------------------------------------------------
 
 my.data_tissue_BAT <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/BAT.txt")
 
 colnames(my.data_tissue_BAT)
 # BAT control general variables   ------------------------------------------------------
 log2<- expression(Log[2])
 fold.thres_tissue_BAT <- 0.5 # User-defined fold-change value for biological significance
 p.thres_tissue_BAT <- 0.05 # User-defined p-val for statistical significance
 p.thres2_tissue_BAT <- 1.31
 neg.fold.thres_tissue_BAT <- -fold.thres_tissue_BAT
 # micos<- c("Apoo","Apool","Immt","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3","Immt")
 # paste(micos,collapse="|")

 # BAT clean up data  -----------------------------------------------------------
 my.data.clean_tissue_BAT <- my.data_tissue_BAT %>%
   # left_join(mitocarta, by = c("Protein IDs" = "UniProt")) %>%
   #rename(`difcol_tissue_BAT`= contains("-Log2 Fold Change")) %>%
   rename(`difcol_tissue_BAT`= contains("Student's T-test Difference")) %>%
   # rename(`FDR` = contains("Student's T-test q-value")) %>%
   #rename(`Significant` = contains("Student's T-test Significant ")) %>%
   rename(`'-log10 p-value_tissue_BAT` = contains("-Log Student's T-test p-value")) %>%
   #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
   #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
   # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
   mutate(`sig_tissue_BAT` = ifelse(`'-log10 p-value_tissue_BAT` > -log10(p.thres_tissue_BAT) & (`difcol_tissue_BAT` > fold.thres_tissue_BAT) | (`'-log10 p-value_tissue_BAT` > -log10(p.thres_tissue_BAT) & (`difcol_tissue_BAT` < neg.fold.thres_tissue_BAT)), "S", "NS")) %>%
   mutate(`label_tissue_BAT` = case_when((`'-log10 p-value_tissue_BAT` > -log10(p.thres_tissue_BAT) & (`difcol_tissue_BAT` < neg.fold.thres_tissue_BAT) ~ "left"),
                                         (`'-log10 p-value_tissue_BAT` > -log10(p.thres_tissue_BAT) & (`difcol_tissue_BAT` > fold.thres_tissue_BAT) ~ "right"),
                                         (`'-log10 p-value_tissue_BAT` > -log10(p.thres_tissue_BAT) ~ "other"))) %>%
   #Remove excess info from 'Gene names'
   mutate(`hasMICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
   mutate(`GenesNew_tissue_BAT` = `Gene names`) %>%
   #mutate(`GenesNew_tissue_BAT` = `Genes`) %>%
   mutate_at("GenesNew_tissue_BAT", list(~ gsub(";.*", " ", .))) %>%
   mutate_at("GenesNew_tissue_BAT", list(~ gsub("-.*", " ", .)))
 
 
 
 # BAT run min/max & control visual variables  ----------------------------------------------------
 
 
 # /* Define axes domains based on max/min values
 round_up <- function(x, to = 10) {
   to * (x %/% to + as.logical(x %% to))
 }
 round_down <- function(x, to = 10) {
   to * (-x %/% to + as.logical(x %% to)) * -1
 }
 
 # x_max_tissue_BAT <- round_up(max(my.data.clean_tissue_BAT $`difcol_tissue_BAT`, na.rm = TRUE), 1)
 x_max_tissue_BAT <- 4.5
 y_max_tissue_BAT <- round_up(max(my.data.clean_tissue_BAT$`'-log10 p-value_tissue_BAT`, na.rm = TRUE), 1)
 x_min_tissue_BAT <- round_down(min(my.data.clean_tissue_BAT$`difcol_tissue_BAT`, na.rm = TRUE), 1)
 y_min_tissue_BAT <- round_up(min(my.data.clean_tissue_BAT$`'-log10 p-value_tissue_BAT`, na.rm = TRUE), 1)
 #y_min <- 0
 y_max_tissue_BAT
 y_min_tissue_BAT
 x_min_tissue_BAT
 x_max_tissue_BAT
 # */
 
 
 
 
 #adjust x limits
 #x_limit_min <- #0 
 #x_limit_min <- `x_min_tissue_BAT` 
 x_limit_min <- 0
 #x_limit_max <- `x_max_tissue_BAT`+2 #adjust after the plus
 x_limit_max <- `x_max_tissue_BAT`
 #adjust y limits
 y_limit_min <- `y_min_tissue_BAT` 
 #y_limit_max <- `y_max_tissue_BAT`+2 #adjust after the plus
 y_limit_max <- `y_max_tissue_BAT`
 
 my_theme_tissue_BAT <- theme_bw() +
   theme(plot.margin=unit(c(2,2,2,2), 'cm'),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.title = element_text(hjust = 0.5),
         panel.border = element_blank(),
         # axis.text.x = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.text.y = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.line = element_line(colour = "black"),
         # axis.title.x = element_text(size = rel(1),
         #                             face = "bold"),
         # axis.title.y = element_text(size = rel(1),
         #                             face = "bold"),
         #              # axis.ticks.x = element_line(size = 4),
         #                 #                         angle = -90),
         # axis.ticks.length=unit(1.5, "cm"),
         # legend.justification = c(1, 0), 
         #legend.position=c(-0.8, -0.15),
         #legend.justification = c("right","bottom"),
         legend.text          = element_text(size = 15),
         # legend.box = "vertical",
         # legend.margin=margin(), 
         #legend.box.just = "right",
         #legend.margin = margin(6, 6, 6, 6),
         legend.box.background = element_rect(colour = "black")
   )
 
 my_layers_tissue_BAT = list(
   labs(title = "Brown Adipose Tissue"),
   # x and y axes labels and scales
   xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
   #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
   ylab(expression(-log[10]~"p-value")),
   scale_x_continuous(
     expand = c(0,0),
     limits = c( `x_limit_min`, `x_limit_max`),
     breaks = seq(`x_min_tissue_BAT`, `x_max_tissue_BAT`, by = 1)), #spacing of axis 
   scale_y_continuous(
     expand = c(0,0),
     limits = c(0, `y_limit_max`)))
 
 my.data.visual_tissue_BAT <- my.data.clean_tissue_BAT %>%
   mutate(`visual_tissue_BAT` = case_when(
     `Gene names` == "Immt" ~ "Bait",
     `sig_tissue_BAT` == "S" & hasMICOS == "MICOS" ~ "MICOS",
     # `sig_tissue_BAT` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
     # `sig_tissue_BAT` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
     `sig_tissue_BAT` == "S" & MitoCarta3.0_List == "+" ~ "MitoCarta3.0",
     `sig_tissue_BAT` == "S" & MitoCarta3.0_List !="+" ~ "Other",
     `sig_tissue_BAT` == "S" ~ "Significant",
     `sig_tissue_BAT` == "NS" ~ "NS"
     
   ))
 
 visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
 visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
 names(visual.colors) <- levels(visual.levels)
 # visual.colors <- c("gray", "black")
 # visual.levels <- factor(levels=c("NS", "Significant"))
 # names(visual.colors) <- levels(visual.levels)
 # metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolism"))
 # metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
 # names(metals.colors) <- levels(metals.levels)
 
 # imputed.levels <- factor(levels=c(0,1,2,3,4))
 # imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
 # names(imputed.colors) <- levels(imputed.levels)
 
 # BAT run visual   --------------------------------------------------------
 # p-value line label - adjust number for position changes 
 line_horizontal_position_x <- `x_max_tissue_BAT` -1 
 line_horizontal_position_y <- 1.31
 line_horizontal_size <- 3 
 # enrichment line label placement
 line_vertical_positive_position_x <- `fold.thres_tissue_BAT`
 line_vertical_positive_position_y <- `y_max_tissue_BAT` - 1#negative line = same
 line_vertical_negative_position_y <- `y_max_tissue_BAT`- 1
 line_vertical_negative_position_x <- `neg.fold.thres_tissue_BAT` 
 
 # enrichment line label to actually label and ad the "x"
 hlabel_format_tissue_BAT<- formattable(2^`fold.thres_tissue_BAT`, digits = 0, format = "f")
 hlabel_tissue_BAT<- paste("x",`hlabel_format_tissue_BAT`, sep='')
 line_vertical_size <- 3
 
 my.plot_tissue_BAT <- my.data.visual_tissue_BAT %>%
   ggplot(mapping = aes(x = `difcol_tissue_BAT`, y= `'-log10 p-value_tissue_BAT`, label = `GenesNew_tissue_BAT`)) + 
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `visual_tissue_BAT` =="Significant"),
     aes(color = visual_tissue_BAT), #black
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `visual_tissue_BAT` =="NS"),
     aes(color = visual_tissue_BAT),#gray
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `visual_tissue_BAT` =="MICOS"),
     aes(color = visual_tissue_BAT), #blue
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `visual_tissue_BAT` =="Bait"),
     aes(color = visual_tissue_BAT), #red
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `sig_tissue_BAT` == "S" & MitoCarta3.0_List == "+" ),
     aes(color = visual_tissue_BAT),
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_BAT, `sig_tissue_BAT` == "S" & MitoCarta3.0_List != "+" ),
     aes(color = visual_tissue_BAT),
     size = 1.5) +
   # gal up 
   # geom_text_repel(  ###TO LABEL ALL DATA POINTS
   #   data= subset (my.data.visual_tissue_BAT, `label_tissue_BAT`== "right"),
   #   aes(color = visual_tissue_BAT),
   #   max.overlaps=1000,
   #   size=3,
   # geom_text_repel(
   #   data= subset (my.data.visual_tissue_BAT, `hasMICOS`== "MICOS"),
   #   aes(color = visual_tissue_BAT),
   #   max.overlaps=1000,
   #   size=3,
 geom_text_repel(
   data= subset (my.data.visual_tissue_BAT,`hasMICOS` =="MICOS"),
   aes(color = visual_tissue_BAT),
   max.overlaps=1000,
   size=4.25, #qqSIZEshape
   #box.padding = unit(2, "lines"),
   point.padding = 0,
   box.padding = 1, #white margin border as a box from corners of plot 
   #position = position_nudge_repel(x = 0, y = 0),
   nudge_x = 0.08,
   #box.padding = 10,
   force = 0.1,  #changes spacing around lines
   # force = 10,  #changes spacing around lines
   #nudge_y = 1,
   #segment.curvature = -0.1,
   #segment.ncp = 3,
   #segment.angle = 90,
   #segment.square= FALSE,
   min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
   seed = 42,
   # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
   # segment.ncp =3,
   # segment.angle = 20,
   segment.size = 0.4, #line thickness
   #point.size = 20
   # Repel just the labels and totally ignore the data points
   # p + geom_text_repel(point.size = NA)
   
   #segment.angle = -180 ## left vs right adjustment 
   #segment.angle = 180 ## left vs right adjustment 
   # Repel away from the left edge, not from the right.
   # xlim = c(NA, Inf),
   # # Do not repel from top or bottom edges.
   # ylim = c(-Inf, Inf)
 )+
   #alter the grey lines
   geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
   #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
   geom_vline(xintercept=`neg.fold.thres_tissue_BAT`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_BAT`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   geom_vline(xintercept=`fold.thres_tissue_BAT`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_BAT`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))
 # facet_zoom(xy = example_values_x > 0.8 & example_values_x < 1.8, zoom.data=zoom) +   # Note the zoom.data argument
 # geom_label_repel(data = dftxt, aes(label = GenesNew_tissue_BAT))
 
 
 # WAT upload   -------------------------------------------
 my.data_tissue_WAT <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/WAT.txt")
 
 colnames(my.data_tissue_WAT)
 
 
 # WAT control general variables  ----------------------------------------------
 log2<- expression(Log[2])
 fold.thres_tissue_WAT <- 0.5 # User-defined fold-change value for biological significance
 p.thres_tissue_WAT <- 0.05 # User-defined p-val for statistical significance
 p.thres2_tissue_WAT <- 1.31
 neg.fold.thres_tissue_WAT <- -fold.thres_tissue_WAT
 micos<- c("Apoo","Apool","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3", "Immt")
 paste(micos,collapse="|")
 
 # clean up data WAT --------------------------------------------------
 my.data.clean_tissue_WAT <- my.data_tissue_WAT %>%
   # left_join(mitocarta, match_fun = str_detect, by = c("Protein IDs" = "UniProt")) %>%
   #rename(`difcol_tissue_WAT`= contains("-Log2 Fold Change")) %>%
   rename(`difcol_tissue_WAT`= contains("Student's T-test Difference")) %>%
   # rename(`FDR` = contains("Student's T-test q-value")) %>%
   #rename(`Significant` = contains("Student's T-test Significant ")) %>%
   rename(`'-log10 p-value_tissue_WAT` = contains("-Log Student's T-test p-value")) %>%
   #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
   #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
   # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
   mutate(`sig_tissue_WAT` = ifelse(`'-log10 p-value_tissue_WAT` > -log10(p.thres_tissue_WAT) & (`difcol_tissue_WAT` > fold.thres_tissue_WAT) | (`'-log10 p-value_tissue_WAT` > -log10(p.thres_tissue_WAT) & (`difcol_tissue_WAT` < neg.fold.thres_tissue_WAT)), "S", "NS")) %>%
   mutate(`label_tissue_WAT` = case_when((`'-log10 p-value_tissue_WAT` > -log10(p.thres_tissue_WAT) & (`difcol_tissue_WAT` < neg.fold.thres_tissue_WAT) ~ "left"),
                                         (`'-log10 p-value_tissue_WAT` > -log10(p.thres_tissue_WAT) & (`difcol_tissue_WAT` > fold.thres_tissue_WAT) ~ "right"),
                                         (`'-log10 p-value_tissue_WAT` > -log10(p.thres_tissue_WAT) ~ "other"))) %>%
   # Remove excess info from 'Gene names'
   mutate(`hasMICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
   mutate(`GenesNew_tissue_WAT` = `Gene names`) %>%
   # mutate(`GenesNew_tissue_WAT` = `Genes`) %>%
   mutate_at("GenesNew_tissue_WAT", list(~ gsub(";.*", " ", .))) %>%
   mutate_at("GenesNew_tissue_WAT", list(~ gsub("-.*", " ", .)))
 
 
 # WAT run min/max & control visual variables   --------------------------------
 
 
 # /* Define axes domains based on max/min values
 round_up <- function(x, to = 10) {
   to * (x %/% to + as.logical(x %% to))
 }
 round_down <- function(x, to = 10) {
   to * (-x %/% to + as.logical(x %% to)) * -1
 }
 
 # x_max_tissue_WAT <- round_up(max(my.data.clean_tissue_WAT $`difcol_tissue_WAT`, na.rm = TRUE), 1)
 x_max_tissue_WAT <- 3.5
 y_max_tissue_WAT <- round_up(max(my.data.clean_tissue_WAT$`'-log10 p-value_tissue_WAT`, na.rm = TRUE), 1)
 x_min_tissue_WAT <- round_down(min(my.data.clean_tissue_WAT$`difcol_tissue_WAT`, na.rm = TRUE), 1)
 y_min_tissue_WAT <- round_up(min(my.data.clean_tissue_WAT$`'-log10 p-value_tissue_WAT`, na.rm = TRUE), 1)
 #y_min <- 0
 y_max_tissue_WAT
 y_min_tissue_WAT
 x_min_tissue_WAT
 x_max_tissue_WAT
 # */
 
 
 #adjust x limits
 #x_limit_min <- #0 
 #x_limit_min <- `x_min_tissue_WAT` 
 x_limit_min <- 0 
 #x_limit_max <- `x_max_tissue_WAT`+2 #adjust after the plus
 x_limit_max <- `x_max_tissue_WAT`
 #adjust y limits
 y_limit_min <- `y_min_tissue_WAT` 
 #y_limit_max <- `y_max_tissue_WAT`+2 #adjust after the plus
 y_limit_max <- `y_max_tissue_WAT`
 
 my_theme_tissue_WAT <- theme_bw() +
   theme(plot.margin=unit(c(2,2,2,2), 'cm'),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.title = element_text(hjust = 0.5),
         panel.border = element_blank(),
         # axis.text.x = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.text.y = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.line = element_line(colour = "black"),
         # axis.title.x = element_text(size = rel(1),
         #                             face = "bold"),
         # axis.title.y = element_text(size = rel(1),
         #                             face = "bold"),
         #              # axis.ticks.x = element_line(size = 4),
         #                 #                         angle = -90),
         # axis.ticks.length=unit(1.5, "cm"),
         # legend.justification = c(1, 0), 
         #legend.position=c(-0.8, -0.15),
         #legend.justification = c("right","bottom"),
         legend.text          = element_text(size = 15),
         # legend.box = "vertical",
         # legend.margin=margin(), 
         #legend.box.just = "right",
         #legend.margin = margin(6, 6, 6, 6),
         legend.box.background = element_rect(colour = "black")
   )
 
 my_layers_tissue_WAT = list(
   labs(title = "White Adipose Tissue"),
   # x and y axes labels and scales
   xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
   #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
   ylab(expression(-log[10]~"p-value")),
   scale_x_continuous(
     expand = c(0,0),
     limits = c( `x_limit_min`, `x_limit_max`),
     breaks = seq(`x_min_tissue_WAT`, `x_max_tissue_WAT`, by = 1)), #spacing of axis 
   scale_y_continuous(
     expand = c(0,0),
     limits = c(0, `y_limit_max`)))
 
 my.data.visual_tissue_WAT <- my.data.clean_tissue_WAT %>%
   mutate(`visual_tissue_WAT` = case_when(
     `Gene names` == "Immt" ~ "Bait",
     `sig_tissue_WAT` == "S" & hasMICOS == "MICOS" ~ "MICOS",
     # `sig_tissue_WAT` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
     # `sig_tissue_WAT` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
     `sig_tissue_WAT` == "S" & MitoCarta3.0_List == "+" ~ "MitoCarta3.0",
     `sig_tissue_WAT` == "S" & MitoCarta3.0_List != "+" ~ "Other",
     `sig_tissue_WAT` == "S" ~ "Significant",
     `sig_tissue_WAT` == "NS" ~ "NS"
     
   ))
 
 visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
 visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
 names(visual.colors) <- levels(visual.levels)
 # visual.colors <- c("gray", "black")
 # visual.levels <- factor(levels=c("NS", "Significant"))
 # names(visual.colors) <- levels(visual.levels)
 # metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolism"))
 # metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
 # names(metals.colors) <- levels(metals.levels)
 
 # imputed.levels <- factor(levels=c(0,1,2,3,4))
 # imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
 # names(imputed.colors) <- levels(imputed.levels)
 
 # WAT run visual  -----------------------------------------------------
 
 # p-value line label - adjust number for position changes 
 line_horizontal_position_x <- `x_max_tissue_WAT` - 1
 line_horizontal_position_y <- 1.31
 line_horizontal_size <- 3 
 # enrichment line label placement
 line_vertical_positive_position_x <- `fold.thres_tissue_WAT`
 line_vertical_positive_position_y <- `y_max_tissue_WAT` - 1 #negative line = same
 line_vertical_negative_position_y <- `y_max_tissue_WAT` - 1
 line_vertical_negative_position_x <- `neg.fold.thres_tissue_WAT` 
 
 # enrichment line label to actually label and ad the "x"
 hlabel_format_tissue_WAT<- formattable(2^`fold.thres_tissue_WAT`, digits = 0, format = "f")
 hlabel_tissue_WAT<- paste("x",`hlabel_format_tissue_WAT`, sep='')
 line_vertical_size <- 3
 
 
 my.plot_tissue_WAT<- my.data.visual_tissue_WAT %>%
   ggplot(mapping = aes(x = `difcol_tissue_WAT`, y= `'-log10 p-value_tissue_WAT`, label = `GenesNew_tissue_WAT`)) + 
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `visual_tissue_WAT` =="Significant"),
     aes(color = visual_tissue_WAT), #black
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `visual_tissue_WAT` =="NS"),
     aes(color = visual_tissue_WAT),#gray
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `visual_tissue_WAT` =="MICOS"),
     aes(color = visual_tissue_WAT), #blue
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `visual_tissue_WAT` =="Bait"),
     aes(color = visual_tissue_WAT), #red
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `sig_tissue_WAT` == "S" & MitoCarta3.0_List == "+" ),
     aes(color = visual_tissue_WAT),
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_WAT, `sig_tissue_WAT` == "S" & MitoCarta3.0_List != "+" ),
     aes(color = visual_tissue_WAT),
     size = 1.5) +
   # geom_text_repel(  ###TO LABEL ALL DATA POINTS
   #   data= subset (my.data.visual_tissue_WAT, `label_tissue_WAT`== "right"),
   #   aes(color = visual_tissue_WAT),
   #   max.overlaps=1000,
   #   size=3,
   geom_text_repel(
     # data= subset (my.data.visual_tissue_WAT, `hasMICOS`== "MICOS"),
     data= subset (my.data.visual_tissue_WAT, `hasMICOS`=="MICOS"),
     aes(color = visual_tissue_WAT),
     max.overlaps=1000,
     size=4.25, #qqSIZEshape
     #box.padding = unit(2, "lines"),
     point.padding = 0,
     box.padding = 1, #white margin border as a box from corners of plot 
     #position = position_nudge_repel(x = 0, y = 0),
     nudge_x = 0.08,
     #box.padding = 10,
     force = 0.1,  #changes spacing around lines
     # force = 10,  #changes spacing around lines
     #nudge_y = 1,
     #segment.curvature = -0.1,
     #segment.ncp = 3,
     #segment.angle = 90,
     #segment.square= FALSE,
     min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
     seed = 42,
     # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
     # segment.ncp =3,
     # segment.angle = 20,
     segment.size = 0.4, #line thickness
     #point.size = 20
     # Repel just the labels and totally ignore the data points
     # p + geom_text_repel(point.size = NA)
     
     #segment.angle = -180 ## left vs right adjustment 
     #segment.angle = 180 ## left vs right adjustment 
     # # Repel away from the left edge, not from the right.
     # xlim = c(NA, Inf),
     # # Do not repel from top or bottom edges.
     # ylim = c(-Inf, Inf)
     
   )+
   #alter the grey lines
   geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
   #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
   geom_vline(xintercept=`neg.fold.thres_tissue_WAT`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_WAT`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   geom_vline(xintercept=`fold.thres_tissue_WAT`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_WAT`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))
 

 # Skeletal muscle upload  --------------------------------------------------
 
 my.data_tissue_SM <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/sm.txt")
 
 colnames(my.data_tissue_SM)
 
 # Skeletal muscle control general variables  ------------------------------------------------------
 log2<- expression(Log[2])
 fold.thres_tissue_SM <- 0.5 # User-defined fold-change value for biological significance
 p.thres_tissue_SM <- 0.05 # User-defined p-val for statistical significance
 p.thres2_tissue_SM <- 1.31
 neg.fold.thres_tissue_SM <- -fold.thres_tissue_SM
 micos<- c("Apoo","Apool","Immt","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3")
 paste(micos,collapse="|")
 
 # Skeletal muscle clean up  -----------------------------------------------------------
 my.data.clean_tissue_SM <- my.data_tissue_SM %>% 
   #left_join(mitocarta, match_fun = str_detect, by = c("Protein IDs" = "UniProt")) %>%
   left_join(mitocarta, by = c("Protein IDs" = "UniProt")) %>%
   #rename(`difcol_tissue_SM`= contains("-Log2 Fold Change")) %>%
   rename(`difcol_tissue_SM`= contains("Student's T-test Difference")) %>%
   # rename(`FDR` = contains("Student's T-test q-value")) %>%
   #rename(`Significant` = contains("Student's T-test Significant ")) %>%
   rename(`'-log10 p-value_tissue_SM` = contains("-Log Student's T-test p-value")) %>%
   #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
   #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
   # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
   mutate(`sig_tissue_SM` = ifelse(`'-log10 p-value_tissue_SM` > -log10(p.thres_tissue_SM) & (`difcol_tissue_SM` > fold.thres_tissue_SM) | (`'-log10 p-value_tissue_SM` > -log10(p.thres_tissue_SM) & (`difcol_tissue_SM` < neg.fold.thres_tissue_SM)), "S", "NS")) %>%
   mutate(`label_tissue_SM` = case_when((`'-log10 p-value_tissue_SM` > -log10(p.thres_tissue_SM) & (`difcol_tissue_SM` < neg.fold.thres_tissue_SM) ~ "left"),
                                        (`'-log10 p-value_tissue_SM` > -log10(p.thres_tissue_SM) & (`difcol_tissue_SM` > fold.thres_tissue_SM) ~ "right"),
                                        (`'-log10 p-value_tissue_SM` > -log10(p.thres_tissue_SM) ~ "other"))) %>%
   #Remove excess info from 'Gene names'
   mutate(`hasMICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
   mutate(`GenesNew_tissue_SM` = `Gene names`) %>%
   #mutate(`GenesNew_tissue_SM` = `Genes`) %>%
   mutate_at("GenesNew_tissue_SM", list(~ gsub(";.*", " ", .))) %>%
   mutate_at("GenesNew_tissue_SM", list(~ gsub("-.*", " ", .)))
 
 #  Skeletal muscle  run min/max & control visual variables----------------------------------------------------
 # /* Define axes domains based on max/min values
 round_up <- function(x, to = 10) {
   to * (x %/% to + as.logical(x %% to))
 }
 round_down <- function(x, to = 10) {
   to * (-x %/% to + as.logical(x %% to)) * -1
 }
 
 # x_max_tissue_SM <- round_up(max(my.data.clean_tissue_SM $`difcol_tissue_SM`, na.rm = TRUE), 1)
 x_max_tissue_SM <- 3.5
 y_max_tissue_SM <- round_up(max(my.data.clean_tissue_SM$`'-log10 p-value_tissue_SM`, na.rm = TRUE), 1)
 x_min_tissue_SM <- round_down(min(my.data.clean_tissue_SM$`difcol_tissue_SM`, na.rm = TRUE), 1)
 y_min_tissue_SM <- round_up(min(my.data.clean_tissue_SM$`'-log10 p-value_tissue_SM`, na.rm = TRUE), 1)
 #y_min <- 0
 y_max_tissue_SM
 y_min_tissue_SM
 x_min_tissue_SM
 x_max_tissue_SM
 # */
 
 #adjust x limits
 #x_limit_min <- #0 
 #x_limit_min <- `x_min_tissue_SM` 
 x_limit_min <- 0
 #x_limit_max <- `x_max_tissue_SM`+2 #adjust after the plus
 x_limit_max <- `x_max_tissue_SM`
 #adjust y limits
 y_limit_min <- `y_min_tissue_SM` 
 #y_limit_max <- `y_max_tissue_SM`+2 #adjust after the plus
 y_limit_max <- `y_max_tissue_SM`
 
 my_theme_tissue_SM <- theme_bw() +
   theme(plot.margin=unit(c(2,2,2,2), 'cm'),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.title = element_text(hjust = 0.5),
         panel.border = element_blank(),
         # axis.text.x = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.text.y = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.line = element_line(colour = "black"),
         # axis.title.x = element_text(size = rel(1),
         #                             face = "bold"),
         # axis.title.y = element_text(size = rel(1),
         #                             face = "bold"),
         #              # axis.ticks.x = element_line(size = 4),
         #                 #                         angle = -90),
         # axis.ticks.length=unit(1.5, "cm"),
         # legend.justification = c(1, 0), 
         #legend.position=c(-0.8, -0.15),
         #legend.justification = c("right","bottom"),
         legend.text          = element_text(size = 15),
         # legend.box = "vertical",
         # legend.margin=margin(), 
         #legend.box.just = "right",
         #legend.margin = margin(6, 6, 6, 6),
         legend.box.background = element_rect(colour = "black")
   )
 
 my_layers_tissue_SM = list(
   labs(title = "Skeletal muscle"),
   # x and y axes labels and scales
   xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
   #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
   ylab(expression(-log[10]~"p-value")),
   scale_x_continuous(
     expand = c(0,0),
     limits = c( `x_limit_min`, `x_limit_max`),
     breaks = seq(`x_min_tissue_SM`, `x_max_tissue_SM`, by = 1)), #spacing of axis 
   scale_y_continuous(
     expand = c(0,0),
     limits = c(0, `y_limit_max`)))
 
 my.data.visual_tissue_SM <- my.data.clean_tissue_SM %>%
   mutate(`visual_tissue_SM` = case_when(
     `Gene names` == "Immt" ~ "Bait",
     `sig_tissue_SM` == "S" & hasMICOS == "MICOS" ~ "MICOS",
     # `sig_tissue_SM` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
     # `sig_tissue_SM` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
     `sig_tissue_SM` == "S" & MitoCarta3.0_List == "MitoCarta3.0" ~ "MitoCarta3.0",
     `sig_tissue_SM` == "S" & is.na(MitoCarta3.0_List) ~ "Other",
     `sig_tissue_SM` == "S" ~ "Significant",
     `sig_tissue_SM` == "NS" ~ "NS"
     
   ))
 
 visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
 visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
 names(visual.colors) <- levels(visual.levels)
 # visual.colors <- c("gray", "black")
 # visual.levels <- factor(levels=c("NS", "Significant"))
 # names(visual.colors) <- levels(visual.levels)
 # metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolism"))
 # metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
 # names(metals.colors) <- levels(metals.levels)
 
 # imputed.levels <- factor(levels=c(0,1,2,3,4))
 # imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
 # names(imputed.colors) <- levels(imputed.levels)
 
 
 #Skeletal muscle run visual  --------------------------------------------------------
 # p-value line label - adjust number for position changes 
 line_horizontal_position_x <- `x_max_tissue_SM` - 1
 line_horizontal_position_y <- 1.31
 line_horizontal_size <- 3 
 # enrichment line label placement
 line_vertical_positive_position_x <- `fold.thres_tissue_SM`
 line_vertical_positive_position_y <- `y_max_tissue_SM` - 1 #negative line = same
 line_vertical_negative_position_y <- `y_max_tissue_SM` - 1 
 line_vertical_negative_position_x <- `neg.fold.thres_tissue_SM` 
 
 # enrichment line label to actually label and ad the "x"
 hlabel_format_tissue_SM<- formattable(2^`fold.thres_tissue_SM`, digits = 0, format = "f")
 hlabel_tissue_SM<- paste("x",`hlabel_format_tissue_SM`, sep='')
 line_vertical_size <- 3
 my.plot_tissue_SM <- my.data.visual_tissue_SM %>%
   ggplot(mapping = aes(x = `difcol_tissue_SM`, y= `'-log10 p-value_tissue_SM`, label = `GenesNew_tissue_SM`)) + 
   geom_point(
     data = subset(my.data.visual_tissue_SM, `visual_tissue_SM` =="Significant"),
     aes(color = visual_tissue_SM), #black
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_SM, `visual_tissue_SM` =="NS"),
     aes(color = visual_tissue_SM),#gray
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_SM, `visual_tissue_SM` =="MICOS"),
     aes(color = visual_tissue_SM), #blue
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_SM, `visual_tissue_SM` =="Bait"),
     aes(color = visual_tissue_SM), #red
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_SM, `sig_tissue_SM` == "S" & MitoCarta3.0_List == "MitoCarta3.0" ),
     aes(color = visual_tissue_SM),
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_SM, `sig_tissue_SM` == "S" & is.na(MitoCarta3.0_List)),
     aes(color = visual_tissue_SM),
     size = 1.5) +
   # gal up 
   # geom_text_repel(  ###TO LABEL ALL DATA POINTS
   #   data= subset (my.data.visual_tissue_SM, `label_tissue_SM`== "right"),
   #   aes(color = visual_tissue_SM),
   #   max.overlaps=1000,
   #   size=3,
   geom_text_repel(
     data= subset (my.data.visual_tissue_SM, `hasMICOS`== "MICOS"),
     aes(color = visual_tissue_SM),
     max.overlaps=1000,
     size=4.25, #qqSIZEshape
     # geom_text_repel(
     #   data= subset (my.data.visual_tissue_SM,`difcol_tissue_SM` >1.8),
     #   aes(color = visual_tissue_SM),
     #   max.overlaps=1000,
     #   size=4,
     #box.padding = unit(2, "lines"),
     point.padding = 0,
     box.padding = 1, #white margin border as a box from corners of plot 
     #position = position_nudge_repel(x = 0, y = 0),
     nudge_x = 0.08,
     #box.padding = 10,
     force = 0.1,  #changes spacing around lines
     # force = 10,  #changes spacing around lines
     #nudge_y = 1,
     #segment.curvature = -0.1,
     #segment.ncp = 3,
     #segment.angle = 90,
     #segment.square= FALSE,
     min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
     seed = 42,
     # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
     # segment.ncp =3,
     # segment.angle = 20,
     segment.size = 0.4, #line thickness
     #point.size = 20
     # Repel just the labels and totally ignore the data points
     # p + geom_text_repel(point.size = NA)
     
     #segment.angle = -180 ## left vs right adjustment 
     #segment.angle = 180 ## left vs right adjustment 
     # Repel away from the left edge, not from the right.
     # xlim = c(NA, Inf),
     # # Do not repel from top or bottom edges.
     # ylim = c(-Inf, Inf)
   )+
   #alter the grey lines
   geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
   #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
   geom_vline(xintercept=`neg.fold.thres_tissue_SM`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_SM`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   geom_vline(xintercept=`fold.thres_tissue_SM`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_SM`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))
 # facet_zoom(xy = example_values_x > 0.8 & example_values_x < 1.8, zoom.data=zoom) +   # Note the zoom.data argument
 # geom_label_repel(data = dftxt, aes(label = GenesNew_tissue_SM))
 
 
 # testes upload    -------------------------------------------
 my.data_tissue_testes <- read_tsv("/Users/vlis2/OneDrive - Monash University/analysis/exp/[mic60_tissueIP]/testes.txt")
 
 colnames(my.data_tissue_testes)
 
 
 # testes control general variables  ----------------------------------------------
 log2<- expression(Log[2])
 fold.thres_tissue_testes <- 0.5 # User-defined fold-change value for biological significance
 p.thres_tissue_testes <- 0.05 # User-defined p-val for statistical significance
 p.thres2_tissue_testes <- 1.31
 neg.fold.thres_tissue_testes <- -fold.thres_tissue_testes
 #micos<- c("Apoo","Apool","Qil1","Minos1","Chchd3","Chchd6","Samm50","Mtx1","Mtx2","Mtx3")
 #paste(micos,collapse="|")
 
 # testes clean up data   --------------------------------------------------
 my.data.clean_tissue_testes <- my.data_tissue_testes %>%
   # left_join(mitocarta, match_fun = str_detect, by = c("Protein IDs" = "UniProt")) %>%
   #rename(`difcol_tissue_testes`= contains("-Log2 Fold Change")) %>%
   rename(`difcol_tissue_testes`= contains("Student's T-test Difference")) %>%
   # rename(`FDR` = contains("Student's T-test q-value")) %>%
   #rename(`Significant` = contains("Student's T-test Significant ")) %>%
   rename(`'-log10 p-value_tissue_testes` = contains("-Log Student's T-test p-value")) %>%
   #rename(`Test statistic` = contains("Student's T-test Test statistic ")) %>%
   #mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` < -log2(fold.thres) | `difcol` > log2(fold.thres)), "S", "NS"))
   # mutate(`sig` = ifelse(`'-log10 p-value` > -log10(p.thres) & (`difcol` > log2(fold.thres)), "S", "NS"))
   mutate(`sig_tissue_testes` = ifelse(`'-log10 p-value_tissue_testes` > -log10(p.thres_tissue_testes) & (`difcol_tissue_testes` > fold.thres_tissue_testes) | (`'-log10 p-value_tissue_testes` > -log10(p.thres_tissue_testes) & (`difcol_tissue_testes` < neg.fold.thres_tissue_testes)), "S", "NS")) %>%
   mutate(`label_tissue_testes` = case_when((`'-log10 p-value_tissue_testes` > -log10(p.thres_tissue_testes) & (`difcol_tissue_testes` < neg.fold.thres_tissue_testes) ~ "left"),
                                            (`'-log10 p-value_tissue_testes` > -log10(p.thres_tissue_testes) & (`difcol_tissue_testes` > fold.thres_tissue_testes) ~ "right"),
                                            (`'-log10 p-value_tissue_testes` > -log10(p.thres_tissue_testes) ~ "other"))) %>%
   # Remove excess info from 'Gene names'
   mutate(`hasMICOS` = case_when(str_detect(`Gene names`, paste(micos,collapse="|")) ~ "MICOS"))  %>%
   mutate(`GenesNew_tissue_testes` = `Gene names`) %>%
   # mutate(`GenesNew_tissue_testes` = `Genes`) %>%
   mutate_at("GenesNew_tissue_testes", list(~ gsub(";.*", " ", .))) %>%
   mutate_at("GenesNew_tissue_testes", list(~ gsub("-.*", " ", .)))
 
 # testes run min/max & control visual variables  --------------------------------
 
 
 # /* Define axes domains based on max/min values
 round_up <- function(x, to = 10) {
   to * (x %/% to + as.logical(x %% to))
 }
 round_down <- function(x, to = 10) {
   to * (-x %/% to + as.logical(x %% to)) * -1
 }
 
 x_max_tissue_testes <- round_up(max(my.data.clean_tissue_testes $`difcol_tissue_testes`, na.rm = TRUE), 1)
 y_max_tissue_testes <- round_up(max(my.data.clean_tissue_testes$`'-log10 p-value_tissue_testes`, na.rm = TRUE), 1)
 x_min_tissue_testes <- round_down(min(my.data.clean_tissue_testes$`difcol_tissue_testes`, na.rm = TRUE), 1)
 y_min_tissue_testes <- round_up(min(my.data.clean_tissue_testes$`'-log10 p-value_tissue_testes`, na.rm = TRUE), 1)
 #y_min <- 0
 y_max_tissue_testes
 y_min_tissue_testes
 x_min_tissue_testes
 x_max_tissue_testes
 # */
 
 #adjust x limits
 #x_limit_min <- #0 
 #x_limit_min <- `x_min_tissue_testes` 
 x_limit_min <- 0
 #x_limit_max <- `x_max_tissue_testes`+2 #adjust after the plus
 x_limit_max <- `x_max_tissue_testes`
 #adjust y limits
 y_limit_min <- `y_min_tissue_testes` 
 #y_limit_max <- `y_max_tissue_testes`+2 #adjust after the plus
 y_limit_max <- `y_max_tissue_testes`
 
 my_theme_tissue_testes <- theme_bw() +
   theme(plot.margin=unit(c(2,2,2,2), 'cm'),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.title = element_text(hjust = 0.5),
         panel.border = element_blank(),
         # axis.text.x = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.text.y = element_text(size = rel(1.5)),
         # #face = "bold"),
         # axis.line = element_line(colour = "black"),
         # axis.title.x = element_text(size = rel(1),
         #                             face = "bold"),
         # axis.title.y = element_text(size = rel(1),
         #                             face = "bold"),
         #              # axis.ticks.x = element_line(size = 4),
         #                 #                         angle = -90),
         # axis.ticks.length=unit(1.5, "cm"),
         # legend.justification = c(1, 0), 
         #legend.position=c(-0.8, -0.15),
         #legend.justification = c("right","bottom"),
         legend.text          = element_text(size = 15),
         # legend.box = "vertical",
         # legend.margin=margin(), 
         #legend.box.just = "right",
         #legend.margin = margin(6, 6, 6, 6),
         legend.box.background = element_rect(colour = "black")
   )
 
 my_layers_tissue_testes = list(
   labs(title = "Testes"),
   # x and y axes labels and scales
   xlab(expression(log[2]~enrichment~Mic60^FLAG/Control)), #add names!!!!
   #xlab(expression(log[2]~mean~ratio~CCDC127^KO/Control)),
   ylab(expression(-log[10]~"p-value")),
   scale_x_continuous(
     expand = c(0,0),
     limits = c( `x_limit_min`, `x_limit_max`),
     breaks = seq(`x_min_tissue_testes`, `x_max_tissue_testes`, by = 1)), #spacing of axis 
   scale_y_continuous(
     expand = c(0,0),
     limits = c(0, `y_limit_max`)))
 
 my.data.visual_tissue_testes <- my.data.clean_tissue_testes %>%
   mutate(`visual_tissue_testes` = case_when(
     `Gene names` == "Immt" ~ "Bait",
     `sig_tissue_testes` == "S" & hasMICOS == "MICOS" ~ "MICOS",
     # `sig_tissue_testes` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*MICOS complex.*") ~"MICOS",
     # `sig_tissue_testes` == "S" & str_detect(MitoCarta3.0_MitoPathways.x,".*SAM.*") ~ "MICOS",
     `sig_tissue_testes` == "S" & MitoCarta3.0_List == "+" ~ "MitoCarta3.0",
     `sig_tissue_testes` == "S" & MitoCarta3.0_List != "+" ~ "Other",
     `sig_tissue_testes` == "S" ~ "Significant",
     `sig_tissue_testes` == "NS" ~ "NS"
     
   ))
 
 visual.colors <- c("blue", "gray", "black", "red", "darkgreen", "darkorange4")
 visual.levels <- factor(levels=c("MICOS", "NS", "Significant", "Bait", "MitoCarta3.0", "Other"))
 names(visual.colors) <- levels(visual.levels)
 # visual.colors <- c("gray", "black")
 # visual.levels <- factor(levels=c("NS", "Significant"))
 # names(visual.colors) <- levels(visual.levels)
 # metals.levels <- factor(levels=c("S", "NS", "Fe-S", "Heme synthesis", "Heme containing", "Copper metabolism"))
 # metals.colors <- c("black", "gray", "red", "orange", "purple", "brown")
 # names(metals.colors) <- levels(metals.levels)
 
 # imputed.levels <- factor(levels=c(0,1,2,3,4))
 # imputed.colors <- c("gray", "yellow", "orange", "red", "purple")
 # names(imputed.colors) <- levels(imputed.levels)
 
 # testes run visual  -----------------------------------------------------
 
 # p-value line label - adjust number for position changes 
 line_horizontal_position_x <- `x_max_tissue_testes` - 1
 line_horizontal_position_y <- 1.31
 line_horizontal_size <- 3 
 # enrichment line label placement
 line_vertical_positive_position_x <- `fold.thres_tissue_testes`
 line_vertical_positive_position_y <- `y_max_tissue_testes` - 1 #negative line = same
 line_vertical_negative_position_y <- `y_max_tissue_testes` - 1
 line_vertical_negative_position_x <- `neg.fold.thres_tissue_testes` 
 
 # enrichment line label to actually label and ad the "x"
 hlabel_format_tissue_testes<- formattable(2^`fold.thres_tissue_testes`, digits = 0, format = "f")
 hlabel_tissue_testes<- paste("x",`hlabel_format_tissue_testes`, sep='')
 line_vertical_size <- 3
 
 my.plot_tissue_testes<- my.data.visual_tissue_testes %>%
   ggplot(mapping = aes(x = `difcol_tissue_testes`, y= `'-log10 p-value_tissue_testes`, label = `GenesNew_tissue_testes`)) + 
   geom_point(
     data = subset(my.data.visual_tissue_testes, `visual_tissue_testes` =="Significant"),
     aes(color = visual_tissue_testes), #black
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_testes, `visual_tissue_testes` =="NS"),
     aes(color = visual_tissue_testes),#gray
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_testes, `visual_tissue_testes` =="MICOS"),
     aes(color = visual_tissue_testes), #blue
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_testes, `visual_tissue_testes` =="Bait"),
     aes(color = visual_tissue_testes), #red
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_testes, `sig_tissue_testes` == "S" & MitoCarta3.0_List == "+" ),
     aes(color = visual_tissue_testes),
     size = 1.5) +
   geom_point(
     data = subset(my.data.visual_tissue_testes, `sig_tissue_testes` == "S" & MitoCarta3.0_List != "+" ),
     aes(color = visual_tissue_testes),
     size = 1.5) +
   # geom_text_repel(  ###TO LABEL ALL DATA POINTS
   #   data= subset (my.data.visual_tissue_testes, `label_tissue_testes`== "right"),
   #   aes(color = visual_tissue_testes),
   #   max.overlaps=1000,
   #   size=3,
   geom_text_repel(
     data= subset (my.data.visual_tissue_testes, `hasMICOS`== "MICOS"),
     # data= subset (my.data.visual_tissue_testes, `difcol_tissue_testes`> 1.8 & `sig_tissue_testes`=="S"),
     aes(color = visual_tissue_testes),
     max.overlaps=1000,
     size=4.25,#qqSIZEshape
     #box.padding = unit(2, "lines"),
     point.padding = 0,
     box.padding = 1, #white margin border as a box from corners of plot 
     #position = position_nudge_repel(x = 0, y = 0),
     nudge_x = 0.08,
     #box.padding = 10,
     force = 0.1,  #changes spacing around lines
     # force = 10,  #changes spacing around lines
     #nudge_y = 1,
     #segment.curvature = -0.1,
     #segment.ncp = 3,
     #segment.angle = 90,
     #segment.square= FALSE,
     min.segment.length = 0,# Draw all line segments OR INF TO NOT DRAW SEGS
     seed = 42,
     # segment.curvature = -0.1, #1 increases right-hand curvature, negative values would increase left-hand curvature, 0 makes straight lines
     # segment.ncp =3,
     # segment.angle = 20,
     segment.size = 0.4, #line thickness
     #point.size = 20
     # Repel just the labels and totally ignore the data points
     # p + geom_text_repel(point.size = NA)
     
     #segment.angle = -180 ## left vs right adjustment 
     #segment.angle = 180 ## left vs right adjustment 
     # # Repel away from the left edge, not from the right.
     # xlim = c(NA, Inf),
     # # Do not repel from top or bottom edges.
     # ylim = c(-Inf, Inf)
     
   )+
   #alter the grey lines
   geom_hline(yintercept=1.31, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_horizontal_position_x` , y =`line_horizontal_position_y`, label = "0.05", vjust = -0.5, size = `line_horizontal_size`, col = 'grey48') +
   #annotate("text", x =-4, y =1.31, label = "0.05", vjust = -0.5, size = 5, col = 'grey48') +
   geom_vline(xintercept=`neg.fold.thres_tissue_testes`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_negative_position_x`, y =`line_vertical_negative_position_y`, label = `hlabel_tissue_testes`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   geom_vline(xintercept=`fold.thres_tissue_testes`, linetype='dashed', col = 'grey48') +
   annotate("text", x =`line_vertical_positive_position_x`, y =`line_vertical_positive_position_y`, label = `hlabel_tissue_testes`, vjust = -0.5, size = `line_vertical_size`, col = 'grey48', angle=90)+
   scale_color_manual("", values = visual.colors)+ guides(color = guide_legend(nrow = 3, byrow = FALSE))
 

# combine -----------------------------------------------------------------
 
 my_theme_combine <- theme(
   # X and Y numbers
   axis.text.x = element_text(size = rel(1.4)),
   #face = "bold"),
   axis.text.y = element_text(size = rel(1.4)),
   #face = "bold"),
   # X and Y writing 
   axis.title.x = element_text(size = rel(1.1), 
                               face = "bold"),
   axis.title.y = element_text(size = rel(1.1), 
                               face = "bold"),
   # ticks 
   # axis.ticks.x = element_line(size = 4),
   #                         angle = -90),
   axis.ticks.length=unit(0.2, "cm"),
   axis.line = element_line(colour = "black")
   # legend.position = 'bottom'
   
   
 )
p1 <- my.plot_tissue_kidney + my_theme_tissue_kidney + my_layers_tissue_kidney + my_theme_combine 
p2 <-  my.plot_tissue_heart + my_theme_tissue_heart + my_layers_tissue_heart + my_theme_combine 
p3 <- my.plot_tissue_BAT + my_theme_tissue_BAT + my_layers_tissue_BAT + my_theme_combine 
p4 <-  my.plot_tissue_WAT + my_theme_tissue_WAT + my_layers_tissue_WAT + my_theme_combine 
p5 <- my.plot_tissue_SM + my_theme_tissue_SM + my_layers_tissue_SM + my_theme_combine 
p6<-  my.plot_tissue_testes + my_theme_tissue_testes + my_layers_tissue_testes + my_theme_combine 

design <- "
123
456
"

firstSet <- p1 + p2 + p3 + plot_layout(guides='collect', nrow =3, ncol=1) &
  theme(legend.position = 'bottom') &
  theme(legend.position = "none")
  coord_cartesian(clip = "off")
     #   (plot.margin = unit(c(3, 3, 3, 3), "cm")))
  firstSet 
  ggsave(paste0("plot_",naming,".pdf"),width = 210, height = 297, units= "mm", dpi = 300)
  
  firstSet <- p4 + p5 + p6 + plot_layout(guides='collect', nrow =3, ncol=1) &
    theme(legend.position = 'bottom') &
    theme(legend.position = "none")
  coord_cartesian(clip = "off")
  #   (plot.margin = unit(c(3, 3, 3, 3), "cm")))
 secondSet 
  ggsave(paste0("plot_2_",naming,".pdf"),width = 210, height = 297, units= "mm", dpi = 300)

# plot_final<- p1 + p2 + p3 + p4 + p5 + p6 +  plot_spacer() + plot_layout(design=design) &
#   theme(legend.position = 'bottom') &
#   theme(legend.position = "none") 
# 
# coord_cartesian(clip = "off")
# #   (plot.margin = unit(c(3, 3, 3, 3), "cm")))
# 
# # plot_final<- p1 + p2 + p3 + p4 + p5 + p6 +  plot_spacer() + plot_layout(guides='collect', nrow =3, ncol=3) &
# #   theme(legend.position = 'bottom') &
# #   theme(legend.position = "none")   
# #   coord_cartesian(clip = "off")
# #      #   (plot.margin = unit(c(3, 3, 3, 3), "cm")))
# plot_final
# ggsave(paste0("plot_",naming,".pdf"),width = 210, height = 297, units= "mm", dpi = 300)





# # test_zoom ---------------------------------------------------------------
# my.data.visual_tissue_kidney <- my.data.visual_tissue_kidney %>%
# mutate(`zoomed1` = case_when((`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` > 1.8 )~ "S")))
#                                          # (`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) & (`difcol_tissue_kidney` > fold.thres_tissue_kidney) ~ "right"),
#                                          # (`'-log10 p-value_tissue_kidney` > -log10(p.thres_tissue_kidney) ~ "other")))
# 
# dftxt <- dplyr::filter(my.dakjjta.visual_tissue_kidney, difcol_tissue_kidney > 0.5 ) %>%
#   dplyr::mutate( zoom = TRUE )      ## All entries to appear in the zoom panel only