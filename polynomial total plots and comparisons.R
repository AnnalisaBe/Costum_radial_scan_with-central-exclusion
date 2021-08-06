


#Author: Annalisa Bellandi, Faulkner group (John Innes centre)

#Title: Polynimial total plots and comparisons between genotypes/treatments

#Content: Given one or multiple series of polynomials (ie various replicates for one or multiple conditions/genotypes) described by 6 coefficients, 
#this script allows to build an average polynomial for each series (condition/genotype). 
#Further, allows to compare the average polynomials with a t-test at each time point.


##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ packages needed  --------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

library('drc')
library('dplyr')
library('ggplot2')
library('reshape2')
library('tidyr')
library('xlsx')
library('ggpubr')
library('purrr')
library('readxl')
library('Rmisc') #important because one key function used below is the summarySE()
library('XLConnect')
library('readxl')
library('tidyverse')
library('R2wd')
library('openxlsx')
library("RColorBrewer")
library("svglite")

##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ set my colour palette --------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

my_pal_div <- RColorBrewer::brewer.pal(11, "BrBG")[2:11]
my_pal_quant_1 <- RColorBrewer::brewer.pal(9, "Oranges")
my_pal_quant_2 <- RColorBrewer::brewer.pal(9, "Blues")
my_pal_gray <- RColorBrewer::brewer.pal(9, "Greys")
okabe_pal <- c("#E69F00","#56B4E9","#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7")

n <- max(length(my_pal_div), length(my_pal_quant_1), length(my_pal_quant_2), length(my_pal_gray), length(okabe_pal))

length(my_pal_div) <- n
length(my_pal_quant_1) <- n
length(my_pal_quant_2) <- n
length(my_pal_gray) <- n
length(okabe_pal) <- n

my_pal_gray_d <- data.frame(my_pal_gray)
my_pal_quant_1_d <- data.frame(my_pal_quant_1)
my_pal_quant_2_d <- data.frame(my_pal_quant_2)
my_pal_div_d <- data.frame(my_pal_div)
okabe_pal_d <- data.frame(okabe_pal)

my_col_set <- (0)
my_col_set <- cbind(my_pal_gray_d, my_pal_quant_1_d)
my_col_set <- cbind(my_col_set, my_pal_quant_2_d)
my_col_set <- cbind(my_col_set, okabe_pal_d)
my_col_set <- cbind(my_col_set, my_pal_div_d)

my_col_set_df <- data.frame(my_col_set)

order <- c(1:10)
my_col_set_df1 <- cbind(my_col_set_df, order)
my_col_set_df1

long_color <- melt(my_col_set_df1,
                   id.vars = "order",
                   variable.name = "palette",
                   value.name = "color")

my_colors_plot <- ggplot(long_color, aes(x = palette, y = order, fill = color)) +
  geom_tile(aes(width=0.93, height=0.95)) +
  scale_fill_identity() +
  scale_y_continuous(breaks=c(1:n)) +
  theme_light()+
  geom_label(aes(label=color),colour = "black", fill= "white", fontface = "bold", size = 4)

my_colors_plot


##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------


##=======================================================================================================================
##=======================================  to manually set before start  ====================================================

#You will have an excel workbook that contains one sheet for each condition/genotype. Each sheet will contain 6 coefficients that describe the polynomial
#corresponding to each replicate


folder_exp1 <- "..." #insert working dyrectory path
exp1 <- "..." #name your experiment

setwd(folder_exp1)
list.files() #lists the files in the wd, check it is correct!

poly_file_max <- "..." #copy paste from the console below to here the name of the excel file that contains the coefficients that describe the polymomials
                       #The Excel workbook will have one sheet for each genotype/condition


#==========================================================================================================================
#================================ calculation of the average polynomial, describing various replicates ====================
#==========================================================================================================================


filename <- poly_file_max

df_total <- c()
df_D_total <- c()

#retrieves numbers of conditions/genotypes
num <- length(getSheetNames(filename))
sheets <- 1:num

wb <- createWorkbook("Polynomial fit_result.xlsx")

#loops over the sheets to calculate an average polynomial and derivative for each sheet
for (k in sheets) {
  
  df <- data.frame(read.xlsx(filename, k))
  GEN <- df[1,1]
  
  xm <- 0:190
  df_reps <- c()
  steps <- rownames(df)
  
  
  #----------------- Polynomial
  
  #loops over replicates of the same sheet to retrieve coefficients and build polynomials
  for (i in steps) {
    y <- df[i,4]*xm + df[i,5]*(xm^2) + df[i,6]*(xm^3) + df[i,7]*(xm^4) + df[i,8]*(xm^5) + df[i,9]*(xm^6)
    replicate <- cbind(xm, y)
    df_reps <- rbind(df_reps, replicate)
  }
  
  
  df_reps <- data.frame(df_reps)
  
  df_reps <- cbind(df_reps, GEN)
  
  head(df_reps,10)
  
  tail(df_reps,10)
  
  #summarySE requires Rmisc
  #creates summary of the postion of the signal at each time point (time=xm)
  summ_df <- summarySE(data=df_reps, measurevar="y", groupvars="xm")
  
  head(summ_df,10)
  
  summ_df <- cbind(summ_df, GEN)
  
  head(summ_df,10)
  
  
  
  #----------------- Derivative
  df_D_reps <- c()
  
  #loops over replicates of the same sheet to retrieve coefficients and build derivatives
  for (i in steps) {
    D <- df[i,4] + 2*(df[i,5]*(xm)) + 3*(df[i,6]*(xm^2)) + 4*(df[i,7]*(xm^3)) + 5*(df[i,8]*(xm^4)) + 6*(df[i,9]*(xm^5))
    replicateD <- cbind(xm, D)
    df_D_reps <- rbind(df_D_reps, replicateD)
  }
  
  df_D_reps <- data.frame(df_D_reps)
  
  df_D_reps <- cbind(df_D_reps, GEN)
  
  head(df_D_reps,10)
  
  #summarySE requires Rmisc
  #creates summary of the velocity of the signal at each time point (time=xm)
  summ_df_D <- summarySE(data=df_D_reps, measurevar="D", groupvars="xm")
  
  head(summ_df_D,10)
  
  summ_df_D <- cbind(summ_df_D, GEN)
  
  head(summ_df_D,10)
  
  
  #-------------------total space vs time: average polynomial
  
  df_total <- rbind(df_total, summ_df)
  
  head(df_total,10)
  
  tail(df_total,10)
  
  
  
  #-------------------total velocity vs time: average derivative
  
  df_D_total <- rbind(df_D_total, summ_df_D)
  
  head(df_D_total,10)
  
  tail(df_D_total,10)
  
  
  poly_reps_sheet <- paste("poly_reps", GEN)      
  
  addWorksheet(wb, poly_reps_sheet)
  writeData(wb, poly_reps_sheet, df_reps)      
  
  poly_D_reps_sheet <- paste("poly_D_reps", GEN)
  
  addWorksheet(wb, poly_D_reps_sheet)
  writeData(wb, poly_D_reps_sheet, df_D_reps)       
  
} #closes for k in num sheet


addWorksheet(wb, 'polynomial_tot')
writeData(wb,'polynomial_tot',df_total)
saveWorkbook(wb, file="Polynomial fit_result.xlsx", overwrite = TRUE)

addWorksheet(wb, 'derivative_tot')
writeData(wb,'derivative_tot',df_D_total)
saveWorkbook(wb, file="Polynomial fit_result.xlsx", overwrite = TRUE)



####--------------plots polynomial

#SE

plot_se <- ggplot(data=df_total, aes(x=xm, y=y, colour = GEN)) +
  geom_errorbar(aes(ymin=y-se, ymax=y+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,500)) + 
  xlab("Time [s]") +
  ylab(paste("Distance from the vein", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  #ggtitle (paste("Avrg dist of the max from the vein over time", "\n", "error bars=SE"))+#, '\n', "N=14")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,3,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_se


#saves PNG version for easy acces
ggsave(paste(EXP, "total distance front_se", ".png"), plot_se, width = 8, height = 6, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "total distance front_se", ".svg"),width=8,height=6)
figure <- plot_se
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total distance front_se", ".pdf"), plot_se, width = 8, height = 6, dpi = 450)



#-----------------plots derivative

#SE

plotD_se <- ggplot(data=df_D_total, aes(x=xm, y=D, colour = GEN)) +
  geom_errorbar(aes(ymin=D-se, ymax=D+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,10)) + 
  xlab("Time [s]") +
  ylab(paste("Velocity [\U003BCm/s]"))+
  theme_bw(base_size=18)+
  #ggtitle (paste("Average velocity over time", "\n", "error bars=SE"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,3,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plotD_se


#saves PNG version for easy acces
ggsave(paste(EXP, "total velocity_se", ".png"), plotD_se, width = 8, height = 6, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "total velocity_se", ".svg"),width=8,height=6)
figure <- plotD_se
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total velocity_se", ".pdf"), plotD_se, width = 8, height = 6, dpi = 450)




#=================================================== Statistical comparisons across genotypes/conditions

#------ for polynomials (space vs time)

wb_test <- createWorkbook("t_test.xlsx")

df_1 <- data.frame(read.xlsx("Polynomial fit_result.xlsx", 1)) #select the sheet where the data for the first item of the comaprison are
df_2 <- data.frame(read.xlsx("Polynomial fit_result.xlsx", 3)) #select the sheet where the data for the second item of the comparison are

head(df_1,10)
tail(df_1,10)

head(df_2,10)
tail(df_2,10)

xm <- unique(df_1[,1])
pvalue <-c()
ttest <- data.frame(cbind(xm, pvalue))
head(ttest,10)

#loops along each time point to compare the values of item 1 and item 2 at each time point
for (i in xm) {
  pvalue <- rbind(pvalue, (t.test(df_1[which(df_1[1]==i), 2], df_2[which(df_2[1]==i), 2])$p.value))
}
pvalue
ttest <- cbind(xm, pvalue)

ttest_poly <- data.frame(ttest)

addWorksheet(wb_test, 'poly_ttest')
writeData(wb_test,'poly_ttest',ttest_poly)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)


#------ for derivative (velocity vs time)
df_1 <- data.frame(read.xlsx("Polynomial fit_result.xlsx", 1)) #select the sheet where the data for the first item of the comaprison are
df_2 <- data.frame(read.xlsx("Polynomial fit_result.xlsx", 3)) #select the sheet where the data for the second item of the comaprison are
head(df_1,10)
tail(df_1,10)

head(df_2,10)
tail(df_2,10)

xm <- unique(df_1[,1])
pvalue <-c()
ttest <- data.frame(cbind(xm, pvalue))
head(ttest,10)

#loops along each time point to compare the values of item 1 and item 2 at each time point
for (i in xm) {
  pvalue <- rbind(pvalue, (t.test(df_1[which(df_1[1]==i), 2], df_2[which(df_2[1]==i), 2])$p.value))
}
pvalue
ttest <- cbind(xm, pvalue)
ttest
ttest_deriv <- data.frame(ttest)

addWorksheet(wb_test, 'deriv_ttest')
writeData(wb_test,'deriv_ttest',ttest_deriv)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)


#=============================================================== Plots

#-------------------------------------------------polynomial

#the basic plot

plot_se_basic <- ggplot()+
  geom_errorbar(data=df_total, mapping=aes(x=xm, ymin=y-se, ymax=y+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_total, mapping=aes(x=xm, y=y, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,300)) + 
  xlab("Time [s]") +
  ylab(paste("Distance from the vein", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Avrg dist of the max from the vein over time", "\n", "error bars=SE"))+#, '\n', "N=14")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_se_basic



#make sure the tttest_poly is a dataframe  
ttest_poly <- data.frame(ttest_poly)


#-----------------------evalaute when the p values calculated above, fall below chosen tresholds

#--------------treshold 0.05 = P05

rl <- rle(ttest_poly[,2]<0.05)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p05 <- data.frame(
  xmin = ttest_poly$xm[starts],
  xmax = ttest_poly$xm[ends],
  type = rl$values
)

df_p05

# Filter on type
df_p05 <- df_p05[df_p05$type == TRUE, ] # Satisfied threshold criterium



##----------------treshold 0.01 = P01

rl <- rle(ttest_poly[,2]<0.01)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p01 <- data.frame(
  xmin = ttest_poly$xm[starts],
  xmax = ttest_poly$xm[ends],
  type = rl$values
)

df_p01

# Filter on type
df_p01 <- df_p01[df_p01$type == TRUE, ] # Satisfied threshold criterium



#--------------------treshold 0.001 = P001

rl <- rle(ttest_poly[,2]<0.001)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p001 <- data.frame(
  xmin = ttest_poly$xm[starts],
  xmax = ttest_poly$xm[ends],
  type = rl$values
)

df_p001

# Filter on type
df_p001 <- df_p001[df_p001$type == TRUE, ] # Satisfied threshold criterium


#----- save these df
addWorksheet(wb_test, 'poly p<0.05')
writeData(wb_test,'poly p<0.05',df_p05)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)

addWorksheet(wb_test, 'poly p<0.01')
writeData(wb_test,'poly p<0.01',df_p01)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)

addWorksheet(wb_test, 'poly p<0.001')
writeData(wb_test,'poly p<0.001',df_p001)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)


#-------------Add shades on the plots that indicate p values falling below the tresholds

plot_p05 <- ggplot()+
  geom_rect(df_p05,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_total, mapping=aes(x=xm, ymin=y-se, ymax=y+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_total, mapping=aes(x=xm, y=y, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,300)) + 
  xlab("Time [s]") +
  ylab(paste("Distance from the vein", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Avrg dist of the max from the vein over time", "\n", "error bars=SE"))+#, '\n', "N=14")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p05

#saves PNG version for easy acces
ggsave(paste(EXP, "total distance max_se_sigplot_p05", ".png"), plot_p05, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
plot_p05
svglite(file=paste(EXP, "total distance max_se_sigplot_p05", ".svg"),width=8,height=8)
figure <- plot_p05
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total distance max_se_sigplot_p05", ".pdf"), plot_p05, width = 8, height = 8, dpi = 450)





plot_p01 <- ggplot()+
  geom_rect(df_p01,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_total, mapping=aes(x=xm, ymin=y-se, ymax=y+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_total, mapping=aes(x=xm, y=y, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,300)) + 
  xlab("Time [s]") +
  ylab(paste("Distance from the vein", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Avrg dist of the max from the vein over time", "\n", "error bars=SE"))+#, '\n', "N=14")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p01

#saves PNG version for easy acces
ggsave(paste(EXP, "total distance max_se_sigplot_p01", ".png"), plot_p01, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
plot_p01
svglite(file=paste(EXP, "total distance max_se_sigplot_p01", ".svg"),width=8,height=8)
figure <- plot_p01
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total distance max_se_sigplot_p01", ".pdf"), plot_p01, width = 8, height = 8, dpi = 450)




plot_p001 <- ggplot()+
  geom_rect(df_p001,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_total, mapping=aes(x=xm, ymin=y-se, ymax=y+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_total, mapping=aes(x=xm, y=y, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,520)) + 
  xlab("Time [s]") +
  ylab(paste("Distance from the vein", "\n" , "[\U003BCm]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Avrg dist of the max from the vein over time", "\n", "error bars=SE"))+#, '\n', "N=14")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p001

#saves PNG version for easy acces
ggsave(paste(EXP, "total distance max_se_sigplot_p001", ".png"), plot_p001, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
plot_p001
svglite(file=paste(EXP, "total distance max_se_sigplot_p001", ".svg"),width=8,height=8)
figure <- plot_p001
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total distance max_se_sigplot_p001", ".pdf"), plot_p001, width = 8, height = 8, dpi = 450)




#----------------------------------------------derivative

ttest_deriv <- data.frame(ttest_deriv)
#if you need to start again after you stopped
#ttest_deriv <- data.frame(read.delim('clipboard', sep = '\t', header = T))

#--------------P05

rl <- rle(ttest_deriv[,2]<0.05)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p05 <- data.frame(
  xmin = ttest_deriv$xm[starts],
  xmax = ttest_deriv$xm[ends],
  type = rl$values
)

df_p05

# Filter on type
df_p05 <- df_p05[df_p05$type == TRUE, ] # Satisfied threshold criterium



##----------------P01

rl <- rle(ttest_poly[,2]<0.01)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p01 <- data.frame(
  xmin = ttest_deriv$xm[starts],
  xmax = ttest_deriv$xm[ends],
  type = rl$values
)

df_p01

# Filter on type
df_p01 <- df_p01[df_p01$type == TRUE, ] # Satisfied threshold criterium



#--------------------P001

rl <- rle(ttest_deriv[,2]<0.001)

ends <- cumsum(rl$lengths)
starts <- ends - rl$lengths+1

df_p001 <- data.frame(
  xmin =ttest_deriv$xm[starts],
  xmax = ttest_deriv$xm[ends],
  type = rl$values
)

df_p001

# Filter on type
df_p001 <- df_p001[df_p001$type == TRUE, ] # Satisfied threshold criterium


##save these df

addWorksheet(wb_test, 'deriv p<0.05')
writeData(wb_test,'deriv p<0.05',df_p05)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)

addWorksheet(wb_test, 'deriv p<0.01')
writeData(wb_test,'deriv p<0.01',df_p01)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)

addWorksheet(wb_test, 'deriv p<0.001')
writeData(wb_test,'deriv p<0.001',df_p001)
saveWorkbook(wb_test, file="t_test.xlsx", overwrite = TRUE)




#-------------Add shades on the plots that indicate p values falling below the tresholds

#basic plot
plot_p05 <- ggplot()+
  geom_rect(df_p05,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_D_total,mapping=aes(x=xm, ymin=D-se, ymax=D+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_D_total, mapping=aes(x=xm, y=D, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,15)) + 
  xlab("Time [s]") +
  ylab(paste("Velocity [\U003BCm/s]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Average velocity over time", "\n", "error bars=SE"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p05

#saves PNG version for easy acces
ggsave(paste(EXP, "total velocity_se_sigplot_p05", ".png"), plot_p05, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "total velocity_se_sigplot_p05", ".svg"),width=8,height=8)
figure <- plot_p05
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total velocity_se_sigplot_p05", ".pdf"), plot_p05, width = 8, height = 8, dpi = 450)







plot_p01 <- ggplot()+
  geom_rect(df_p01,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_D_total,mapping=aes(x=xm, ymin=D-se, ymax=D+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_D_total, mapping=aes(x=xm, y=D, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,15)) + 
  xlab("Time [s]") +
  ylab(paste("Velocity [\U003BCm/s]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Average velocity over time", "\n", "error bars=SE"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p01

#saves PNG version for easy acces
ggsave(paste(EXP, "total velocity_se_sigplot_p01", ".png"), plot_p01, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "total velocity_se_sigplot_p01", ".svg"),width=8,height=8)
figure <- plot_p01
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total velocity_se_sigplot_p01", ".pdf"), plot_p01, width = 8, height = 8, dpi = 450)



plot_p001 <- ggplot()+
  geom_rect(df_p001,
            mapping=aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf), fill=my_pal_gray[2], alpha=0.5)+
  geom_errorbar(data=df_D_total,mapping=aes(x=xm, ymin=D-se, ymax=D+se), width=1, color=my_pal_gray[5], alpha=0.5) +  
  geom_line(data=df_D_total, mapping=aes(x=xm, y=D, colour = GEN), size=1.4, linetype=1) +
  coord_cartesian(xlim=c(0,180), ylim=c(0,15)) + 
  xlab("Time [s]") +
  ylab(paste("Velocity [\U003BCm/s]"))+
  theme_bw(base_size=18)+
  ggtitle (paste("Average velocity over time", "\n", "error bars=SE"))+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank())+
  scale_color_manual(values=okabe_pal[c(1,2,4)]) +
  theme(text =element_text(size=20)) +
  #theme(legend.key.size = unit(1.5, "cm"),legend.key.width = unit(0.5,"cm"), legend.spacing.y = unit(1, "cm"))+
  theme(legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(1.5, "cm"), legend.position = 'bottom')

plot_p001

#saves PNG version for easy acces
ggsave(paste(EXP, "total velocity_se_sigplot_p001", ".png"), plot_p001, width = 8, height = 8, dpi = 450, device = 'png')

#saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
#strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
svglite(file=paste(EXP, "total velocity_se_sigplot_p001", ".svg"),width=8,height=8)
figure <- plot_p001
print(figure)
dev.off()

#saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
#illustrator, however this renders correctly in inkscape
ggsave(paste(EXP,"total velocity_se_sigplot_p001", ".pdf"), plot_p001, width = 8, height = 8, dpi = 450)




