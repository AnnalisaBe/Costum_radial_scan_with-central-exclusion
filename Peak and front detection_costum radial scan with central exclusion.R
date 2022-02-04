
#Author: Annalisa Bellandi, Faulkner group (John Innes centre)

#Title: Peak and front detection after costum radial scan with central exclusion

#Content: given data output of the costum radila scan with central exclusion fiji macro
#retrieves number of radii and area of the central exclusion
#calculates average profile of the sginal at each time point across all radii
#calculates and accounts for delay between the application of teh stimulus and the detection of the response
#detects peak and front of the wave
#outputs polynomials to describe the distance reached by the signal over time and its velcoity over time for each replicate of each genotype



##-------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ packages needed  --------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

#packages needed
library('rlang')
library("purrr")
library("tidyr")
library("Rmisc")
library("reshape2")
library("dplyr")
library("openxlsx")
library("ggplot2")
library("zoo")
library("RColorBrewer")
library("svglite")
library('ggsignif')
library('ggpubr')


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




##=======================================================================================================================
##=======================================  to manually set before starting  ====================================================

EXP <- '...' #experiment name


experiment_folder <- "..." #path to the folder that contains all genotypes folder


genotype1 <- "..." #path to the folder that contains data for teh genotype 1

GEN1 <-  "..." #name of the genotype 1 as it is in displayed the files



genotype2 <- "..." #path to the folder that contains data for teh genotype 2

GEN2 <-  "..." #name of the genotype 2 as it is in displayed the files




CON <- "..." #condition of the experiment (ie cotyleodns or detached leaves... etc)



#---------Things to adjust only if necessary: decide the optimal spans, window width and treshold
w=4
span1=0.2 
span2=0.5
treshold <- 1.96

#note there is no pre-defined number of radii, number of radii are stored in the fiji output and retrieved for each sample

##=============================================================================================================================
##-----------------------------------------------------------------------------------------------------------------------------
##------------------------------------------ defining FUNCTIONS ---------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------

#----------------------------------------- FUNCTION 1

##this function fits a loess curve on the data, then finds the max and then saves their position along the x.
#All of this, on one time point


##this function fits a loess curve on the data, then finds the max and then saves their position along the x.
#All of this, on one time point

d_max <- function(x, y, w, span) {
  require(zoo)
  n <- length(y)
  
  #the span influences how smooth the spline will be, this will influence also the detection of the max, 
  #the bigger the value the smoother the loess, ie little bumps will be ignored
  # https://rafalab.github.io/dsbook/smoothing.html best explanation of loess fitting
  y_smooth <- loess(y ~ x, span=span)$fitted
  
  #zoo lists the values of y_smooth assigning to them an increasing index, then we divide the y_smooth in small windows 
  #of  2*w+1 width, within those windows we find the max value of y_smooth, we assign this value to the index of the 
  #centre of the window. As the window rolls down the y_smooth data we form a y_max list of values consituted by the 
  #max value of each window with the indexes corresponding to the centre of the window
  #how rollapply works: https://stackoverflow.com/questions/50014268/shifting-of-running-mean-output-from-the-rollapply-function-in-r
  y_max <- rollapply(zoo(y_smooth), 2*w+1, max, align="center")
  
  #now we want to compare the y_max values (the max of each window) and the y_smooth values so that to determine the exact position of the
  #maxima, but to do so we need to remember that y_max indexes starts 
  #form the centre of the first window, so first we need to centre y_max on y_smooth and trim y_smooth to exclude all the 
  #values before the centre of the first window and after the centre of the last window.
  y_smooth_trimmed <- y_smooth[-c(1:w, n+1-1:w)]
  #y_smooth_trimmed now starts with the value in position w+1, which is the centre of the first window used in the rollapply command. 
  #now we can compare the two dataset - a max listed in the y_max it is a real max only if in that pont the loess function is actually equal
  #or bigger than the value of the max.This comparison allow us to get out 
  #only the real local maxima (see graphic rapresentation to understand it)
  delta <- y_max - y_smooth_trimmed
  
  
  #the indexes of our max values now refer to the indexes of the y_max, but we want to bring them back to the indexes of our original 
  #dataset (y_smooth) so that we can fish out the correct distance from the centre that corresponds to each max, so we add w. 
  index_max_all <- which(delta <= 0) + w
  
  #subset the list of all the detected peaks  considering only the peaks in forward, this 
  #should guarantee that we follow only the major peak and we don't get distracted with peaks from the
  #blind side
  index_max_all_maj <- subset(index_max_all, index_max_all>=(peaks$i-2))
  
  #sometimes at later time point, the detected peak jumps forward without a real reason, even if the real peak is sunk into
  #the noise. This can mess up my polynomial fitting, so I here I fix it by taking advantage of the fact that at later
  #time points it is unusual so the max to jump many steps forward.
  if (j>20) {
    #not too far away
    index_max_all_maj <- subset(index_max_all, index_max_all<=(peaks$i+65))
    #not too far back
    index_max_all_maj <- subset(index_max_all_maj, index_max_all_maj>=(peaks$i-2))
  }
  
  if (j>115) {
    #not too far away
    index_max_all_maj <- subset(index_max_all, index_max_all<=(peaks$i+15))
    #not too far back
    index_max_all_maj <- subset(index_max_all_maj, index_max_all_maj>=(peaks$i-2))
  }
  
  #this ensures that if there are no peaks in front, the detected wavefront remains in the previous position
  if (length(index_max_all_maj)==0) {index_max_all_maj <- peaks$i}
  
  #creates a df with the values of the fitted loess in one column and the indexing in the other
  intensity <- y_smooth_trimmed
  index_num <- index(y_max)
  df_y_smooth_trimmed <- data.frame(cbind(index_num, intensity))
  
  #in some datasets, for unknown reasons, the lenght of the rays it is not always the same at each time point
  #in this situation, if the previous max is further away than the legth of the current rays, script would crash beacuse searching
  #for max at a distance higher than the available at the moment and so it reamins without index of the max.
  #This line assigns the highet possible distance to the index in such case
  if (min(index_max_all_maj) > max(index_num)) {index_max_all_maj <- max(index_num)}
  
  #subsets the df to extract the rows containing the maxima
  max_values_all <- subset(df_y_smooth_trimmed, index_num %in% index_max_all_maj)
  
  #In normal conditions:
  #find highest max, to discard smaller scondary peaks
  max_value_main <- max(max_values_all$intensity)
  
  #subset dataframe to have a row that contains both the value of the chosen max and the index of it
  max_value <- df_y_smooth_trimmed[which(intensity==max_value_main),]
  
  
  #Sometimes at the very beginning of the time series the peak it is not correctly detected because:
  #-a) falls in the blind spot at the beginning of the data (so the index is < w+1)
  #-b) it is out of the blind zone at 
  #     the beginning of the time series, but it is not a well enough formed peak to be recognised as a max 
  #-c) there is a clear peak, but it is lower than the central peak
  
  if (j<30) {
    #case c - if more than one peak is detected in the first 1/3 of the distance, most likely the one closer to the centre
    #is the central peak and the second furthest away is the real peak. To individuate it:
    
    #subsets the max value all and selects only the peaks in the first 1/3 of the distance
    max_values_all_3 <- subset(max_values_all, index_num <= max(index(y_max))/3)
    
    #if the detected peaks in that section are more than one the real peak could be ignored if it is less intense
    #than the wounding site, so if the peaks are more than one....
    
    if ( nrow(max_values_all_3)>1) { 
      
      #normally if the real peak and the central peak are both there, they are pretty close in intensity, the following
      #'if' makes sure that we don't consider situations where we have a nice clear peak and some smaller peaks further
      #away. So we consider two peaks only if their intensities are less than 15% of the higher peak away from each other
      
      int_peak1 <- max_values_all_3[1,2]
      
      int_peak2 <- max_values_all_3[2,2]
      
      if (abs(int_peak1-int_peak2) <= (int_peak1*0.15)) {
        
        n <- length(max_values_all_3$intensity)
        #sorts the peaks according to the dispance from the centre and choose the second one
        second_peak_ind <- sort(max_values_all_3$index_num)[2]
        #retreive the intensity of that peak, we don't care if it is less intense
        second_peak_int <- max_values_all[which(max_values_all$index_num==second_peak_ind), "intensity"]
        #the df_y_smooth_trimmed is not subsetted based on the main peak, but the second peak
        max_value <- df_y_smooth_trimmed[which(intensity==second_peak_int),]
        
      }
    }
    
    #If the detected max is further than half of the distance, the real peak went missing eaither for reason a or b
    
    
    if (max_value$index_num > max(index_num)/2) {
      
      index_max_all_maj <- peaks$i
      
      max_values_all <- subset(df_y_smooth_trimmed, index_num %in% index_max_all_maj)
      
      max_value_main <- max(max_values_all$intensity)
      
      #subset dataframe to have a row that contains both the value of the chosen max and the index of it
      max_value <- df_y_smooth_trimmed[which(intensity==max_value_main),]
    }
    
    #if the detected maxima is closer than the half of the distance but is not high, the real peak went missing for reason a or b
    
    if (max_value$intensity < (max(y_smooth))/2) {
      
      index_max_all_maj <- peaks$i
      
      max_values_all <- subset(df_y_smooth_trimmed, index_num %in% index_max_all_maj)
      
      max_value_main <- max(max_values_all$intensity)
      
      #subset dataframe to have a row that contains both the value of the chosen max and the index of it
      max_value <- df_y_smooth_trimmed[which(intensity==max_value_main),]
    }
    
  }
  
  
  
  ##---------------------added bit for the front detection
  
  
  
  #newline
  #insert here part to find the front of the signal, from the ysmooth I think find where the values are higher that the baseline+sd. The further value at the end of a continuus line of good values of maybe 10 values?
  surface <- baseline[1]+base_sd[1]
  float <- y_smooth - surface
  
  rl <- rle(float>0)
  
  rl
  
  ends <- cumsum(rl$lengths)
  starts <- ends - rl$lengths+1
  ends
  
  starts
  
  df <- data.frame(xmin=starts, xmax=ends, type=rl$values)
  df
  
  df_T<- df[df$type == TRUE, ] 
  df_T
  
  if (dim(df_T)[1] == 0) { #if there is nothing floating on top of the baseline, just wait where you are (default wait before start is the edge of the window)
    
    front_index <- peaks$front_ind
    
  } else
    
  {
    
    wide <- df_T$xmax-df_T$xmin
    
    df_w <- data.frame(cbind(df_T, wide))
    
    long_int <- which(df_w$wide>10) #which are the long intervals
    
    
    
    if (length(long_int)>1) { #if there are more than one long lintervals floating on the surface, then choose the bigger one
      
      long_int <- max(long_int)
      front_index <- df_w[long_int, "xmax"]
      
      #but if the bigger one is behind, then you are ion the centre of the doughnut, so choose the second long intervall             
      
      if (front_index<(peaks$front_ind-20)) {
        long_int <- 2 
        front_index <- df_w[long_int, "xmax"]
        
        
        #but if the smaller one is very far away in front compared to the previous front, maybe this is not a legitimate jump
        #i.e. a random noisy peak may have created a front further away
        #in this case i would stick with the first intervall
        
        if (front_index>(peaks$front_ind+20)) {
          long_int <- 1 
          front_index <- df_w[long_int, "xmax"]
        }
        
        #of course this can't be valid at the very start, because the start front_ind is w+1 and the front is often quite far from there 
        #even at the startbut at early timem, there shouldn't be two intervals, so it should be fine
        
      }
    } else {
      
      
      if ((length(long_int))==1) { #if there is a good long intervall, just choose that
        
        long_int <- 1 
        front_index <- df_w[long_int, "xmax"]
      }
      
      
      if ((length(long_int))==0) { #if there is no long interval floating on the surface, then it means that maybe it is early time, so just wait at the edge of the first
        front_index <- peaks$front_ind
      }
      
    }
    
  }
  
  #when the wave starts fading there are problems because the front would receed, so at later times I will force it to stay in place as to mark how far the signal went
  #what i am saying is'try to find the peak, but if is late times when the wave is likely to be fading and what you find it is too far behind, then stay in place
  
  #if (j>150) { this was initially a time-dependent idea, but then I realize the front should never go back
  
  if (front_index<(peaks$front_ind)) {
    
    front_index <- peaks$front_ind
    
  }
  
  # }
  
  
  #in general, apart the first 15 sec when we may have a front already far away from the base front_index set at the start
  #, the front shouldn't make big jumps in front either and if it does is because some randon noise is surpassing the threshold
  # so i constrict the front to not jump too far, ideally as I do with the peak
  
  if(j>15) {
    if (front_index>(peaks$front_ind+35)) {
      front_index <- peaks$front_ind
    }
  }
  
  
  
  #in any case the peak shouldn't be further away than the front. 
  #So if the detected peak is further than the front, then just stay where the previous peak was. This will save me from any other situations where the peak may bounce forward for some reason
  
  if (max_value$index_num>front_index) {
    
    index_max_all_maj <- peaks$i
    
    max_values_all <- subset(df_y_smooth_trimmed, index_num %in% index_max_all_maj)
    
    max_value_main <- max(max_values_all$intensity)
    
    #subset dataframe to have a row that contains both the value of the chosen max and the index of it
    max_value <- df_y_smooth_trimmed[which(intensity==max_value_main),]
    
  }
  
  
  #this lists:
  #index_max: the value of the x (profile_t_davg$d) that correspond to the index of the max (so my distance of the max 
  #from the centre!)
  #i: the index of those maxima in the y_max list
  #y.hat: the values of intesity of y_smooth
  #y.values: the value of intesity of y_smooth that correspond to the maxima (where do you want your red ball, what 
  #is the intesnity of the max)
  list(x_max=x[max_value$index_num], i=max_value$index_num, y.hat=y_smooth, y_values=max_value$intensity, baseline=baseline[1], base_sd=base_sd[1], x_front=x[front_index], y_front=y_smooth[front_index], front_ind=front_index)
}



#----------------------------------------- FUNCTION 2

#this is the function that calculates the noise for one replicate

estimate_noise <- function(treshold, avg_profile) { 
  
  #treshold is how many times the sd you want to add to the baseline to reach the real peak level
  
  ## to calculate the noise we need an average sd and a baseline
  sd <- (avg_profile$se) * (sqrt(avg_profile[1,"n_rad"])) 
  
  #average standard deviation across all data 
  avg_sd <- sum(sd)/length(sd)
  
  #baseline, calculated as the mean of the fluorescence in the first 30 seconds at a distance 
  #from the centre higher than 300 um
  avg_profile_base <- avg_profile[avg_profile$d>350,]
  avg_profile_base <- avg_profile_base[avg_profile_base$t<20,]
  
  baseline <- mean(avg_profile_base$avrg)
  
  real_peak <- baseline+avg_sd*treshold
  
  s_n <- avg_profile$avrg/real_peak
  avg_profile_sn <- data.frame(cbind(avg_profile,s_n))
  
  #newline
  avg_profile_sn <- data.frame(cbind(avg_profile_sn,baseline))
  avg_profile_sn <- data.frame(cbind(avg_profile_sn,avg_sd))
  
  addWorksheet(workbook_profiles, NAME)
  writeData(workbook_profiles, NAME , avg_profile_sn)
  saveWorkbook(workbook_profiles, file=paste(GEN, "profiles", ".xlsx"), overwrite = TRUE)
  
  return(avg_profile_sn)
  
}

#----------------------------------------- FUNCTION 3
#this is the function that: runs the d_max function, makes plots and saves them

fit_loess_peaks <- function(w, span, j, SE, x, profile_t) {
  
  #this saves the list tha is the output of the function d_max under the name "peaks"
  peaks <- d_max(x, y, w=w, span=span)
  
  max_j <- data.frame(cbind(t=j, distance_max=peaks$x_max, intensity_max=peaks$y_values))
  front_j <-data.frame(cbind(t=j, distance_front=peaks$x_front, intensity_front=peaks$y_front, baseline=peaks$baseline, base_sd=peaks$base_sd))
  
  model_j <- data.frame(cbind(t=j, space=x, profile_smooth=peaks$y.hat))
  
  ribbon_min <- (front_j$baseline)-(front_j$base_sd)
  
  ribbon_max <- (front_j$baseline)+(front_j$base_sd)
  
  y_axes_min <- ribbon_min-(front_j$base_sd)
    
  profile_t1<- profile_t %>% mutate(Group = cut(s_n, breaks = c(0,1.06,1.5,2, Inf), include.lowest = TRUE))
  
    ggplot() +
    geom_errorbar(profile_t1, 
                  mapping=aes(x=d, ymin=avrg-se, ymax=avrg+se), width=0.1, color="gray") +
    geom_point(profile_t1, 
               mapping=aes(x=d, y=avrg, fill=Group),
               shape=21, col=my_pal_gray[6], stroke=0.01, size=2.5, alpha=1) +
    xlab(paste("Distance from the edge of the droplet", "\n" , "[\U003BCm]")) +
    ylab("Fluorescence intensity") + #use unicode in ggplot for greek letters
    theme_bw(base_size=18)+
    ggtitle(paste("t = ", j, ", w = ", w, ", span = ", span, sep=""))+
    labs(fill="signal/noise")+
    
    scale_fill_manual(breaks = levels(profile_t1$Group), drop = FALSE, values = my_pal_quant_2[c(1,3,5,8)],
                      labels=c(expression("s/n" <= "1"), "1 < s/n < 1.5", expression("1.5"<="s/n < 2"), expression("s/n" >= "2")))+
    
    geom_line(data=model_j,mapping=aes(x=model_j$space, y=model_j$profile_smooth), color=my_pal_gray[9], lty=1, lwd=0.7, alpha=0.4)+
    
    #baseline line
    geom_hline(data=front_j,mapping=aes(yintercept=front_j$baseline), color=my_pal_gray[9], lty=1, lwd=0.7, alpha=0.4)+
    
    #ribbon for the baseline
    annotate('ribbon', x = c(-Inf, Inf), ymin = ribbon_min[1], ymax =ribbon_max[1], alpha = 0.2, fill = okabe_pal[1])+
   
    #segment for the max
    geom_segment(data=max_j, mapping=aes(x = distance_max, y = y_axes_min , xend = distance_max, yend = intensity_max), colour = my_pal_gray[9],lty=2, lwd=0.7)+
    
    #segment for the front
    geom_segment(data=front_j, mapping=aes(x = distance_front, y = y_axes_min , xend = distance_front, yend = intensity_front), colour = my_pal_gray[9],lty=2, lwd=0.7)+
    
    guides(aesthetics = "fill", fill = guide_legend(reverse = TRUE, override.aes = list(shape = 21, size= 10))) +
    
    #point for the max
    geom_point(data=max_j, mapping=aes(x=distance_max, y=intensity_max, fill="detected /n peak"), shape=23, size=3.5, fill=my_pal_quant_1[6])+
  
    #point for the front
    geom_point(data=front_j, mapping=aes(x=distance_front, y=intensity_front, fill="detected /n front"), shape=24, size=3.5, fill=okabe_pal[3])+
    
    scale_y_continuous(limits = c(y_axes_min, NA))
    
  ggsave(file = file.path(plotdir, paste("Max_time_point", j, ".png", sep = "")),width = 8, height = 5, dpi = 450, device = 'png')
  
  #dev.off()
  
  
  
}


#-------------------------------FUNCTION 4
#from the previous script version, creates the average_profile dataset according to the wound frame and stable frmate 
#number, starting from the fiji output


wave_profile <- function(PLANT, GEN) {
  
  NAME <- paste(GEN," ", "series", " ", PLANT,".csv", sep="")
  
  df <- read.csv(NAME, header=T)
  tail(df,5)
  head(df,10)
  
  #cut away anything that is before the wound event
  dftrim <- subset(df, frame>=woundframe)
  dftrim <- dftrim[-1]
  head(dftrim, 10)
  
  #sets the wonding time as zero
  stime <- min(dftrim[,"t"])
  dftrimt0 <- transform(dftrim, t = t-stime)
  
  #remove times higher than 350 sec, as I know that sometimes I recorded longer but not much happens further than that time and analysing it would
  #just slow down the analysis
  dftrimt0 <- subset(dftrimt0, t<=350) 
  
  #removing the data before the stability saves from detecting the wrong max beacuse the leaf is not stable and then remaining anchored in the wrong location
  dftrimt0 <- subset(dftrimt0, frame>=stableframe) 
  head(dftrimt0, 10)
  tail(dftrimt0, 10)
  
  # ------------------ profile of the wave
  
  dflong<-melt(dftrimt0, id=c("t","frame","d", "r.drop", "n.radii"))
  head(dflong, 10)
  
  n_rad=floor(dflong[1,"n.radii"])
  n_rad
  r_drop=dflong[1,"r.drop"]
  r_drop
  
  dflong_edge <- subset(dflong, d>=r_drop)
  head(dflong_edge, 10)
  dflong_edge0 <- transform(dflong_edge, d = d-r_drop)
  head(dflong_edge0, 10)
  
  avg_profile <-  dflong_edge0 %>% dplyr::group_by(t, d) %>% summarise(avrg=mean(value), se=sd(value)/sqrt(n_rad))
  head(avg_profile,10)
  
  avg_profile <- data.frame(avg_profile)
  avg_profile <- cbind(avg_profile,n_rad)
  avg_profile <- cbind(avg_profile, r_drop)
  head(avg_profile,10)
  
  return(avg_profile)
  
}



#-------------------------------FUNCTION 5
#from the previous script version
#fits 6th grade polynomial on the average distance of the max from the centre over time
#calculates derivatives
#plots all

Fit_poly_6_max <- function (df_maxima, my_pal_gray, base_plot) {
  
  model_poly6 <- lm(df_maxima$distance_max ~ poly(df_maxima$t,6, raw = T) -1) # the -1 forces the poly to pass from 0,0
  cf0 = coef(model_poly6)
  cf0
  
  xm <- 0:max(df_maxima$t)
  y0 <- cf0[1]*xm + cf0[2]*(xm^2) + cf0[3]*(xm^3) + cf0[4]*(xm^4) + cf0[5]*(xm^5) + cf0[6]*(xm^6)
  model_poly6_df <- data.frame(cbind(xm,y0))
  
  
  coefficient_poly <- data.frame(genotype=GEN, sample_name=PLANT, condition=CON)
  
  coefficient_poly$coeff1 <- cf0[1]
  coefficient_poly$coeff2 <- cf0[2]
  coefficient_poly$coeff3 <- cf0[3]
  coefficient_poly$coeff4 <- cf0[4]
  coefficient_poly$coeff5 <- cf0[5]
  coefficient_poly$coeff6 <- cf0[6]
  
  coefficient_poly
  
  coefficient_poly6 <- coefficient_poly
  
  #now plot the model line on top and then save
  
  poly6_plot <- base_plot +
    geom_line(data=model_poly6_df, mapping=aes(x=model_poly6_df$xm, y=model_poly6_df$y0),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  poly6_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 6th", ".png"), poly6_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "polynomial constrained fit 6th", ".svg"),width=8,height=5)
  #figure <- poly6_plot
  #print(figure)
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 6th", ".pdf"), poly6_plot, width = 8, height = 5, dpi = 450)
  
  
  
  #---------------find the derivative of my fitted polynomial
  
  #first derivative of my model - this gives us the velocity
  xmD1 <- 0:max(df_maxima$t)
  
  D1 = cf0[1] + 2*cf0[2]*xmD1 + 3*cf0[3]*xmD1^2 + 4*cf0[4]*xmD1^3 + 5*cf0[5]*xmD1^4 + 6*cf0[6]*xmD1^5
  
  model6tD1_df <- data.frame(cbind(xmD1,D1))
  
  poly6_D1 <- ggplot()+
    geom_line(data=model6tD1_df, mapping=aes(x=xmD1, y=D1),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
    coord_cartesian(xlim=c(0,300), ylim=c(0,20)) + 
    xlab("Time [s]") + 
    ylab("Velocity [\U003BCm/s]") + 
    ggtitle (paste("Velocity of the wave over time", '\n', GEN, "series", PLANT, '\n', CON)) + 
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none")
  
  
  poly6_D1
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "derivative poly 6th", ".png"), poly6_D1, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "derivative poly 6th", ".svg"),width=8,height=5)
  #figure <- poly6_plot
  #print(figure)
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "derivative poly 6th", ".pdf"), poly6_D1, width = 8, height = 5, dpi = 450)
  
  
  #this is important beacuse R normally does not remember what happens inside a function, but I can chose one thing to remember as output of a function
  return(coefficient_poly6)
}



#from the previous script version
#fits 6th grade polynomial on the average distance of the front from the centre over time
#calculates derivatives
#plots all

Fit_poly_6_front <- function (df_front, my_pal_gray, base_plot) {
  
  model_poly6 <- lm(df_front$distance_front ~ poly(df_front$t,6, raw = T) -1) # the -1 forces the poly to pass from 0,0
  cf0 = coef(model_poly6)
  cf0
  
  xm <- 0:max(df_maxima$t)
  y0 <- cf0[1]*xm + cf0[2]*(xm^2) + cf0[3]*(xm^3) + cf0[4]*(xm^4) + cf0[5]*(xm^5) + cf0[6]*(xm^6)
  model_poly6_df <- data.frame(cbind(xm,y0))
  
  
  coefficient_poly <- data.frame(genotype=GEN, sample_name=PLANT, condition=CON)
  
  coefficient_poly$coeff1 <- cf0[1]
  coefficient_poly$coeff2 <- cf0[2]
  coefficient_poly$coeff3 <- cf0[3]
  coefficient_poly$coeff4 <- cf0[4]
  coefficient_poly$coeff5 <- cf0[5]
  coefficient_poly$coeff6 <- cf0[6]
  
  coefficient_poly
  
  coefficient_poly6 <- coefficient_poly
  
  #now plot the model line on top and then save
  
  poly6_plot_front <- base_plot +
    geom_line(data=model_poly6_df, mapping=aes(x=model_poly6_df$xm, y=model_poly6_df$y0),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  poly6_plot_front
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 6th_front", ".png"), poly6_plot_front, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "polynomial constrained fit 6th", ".svg"),width=8,height=5)
  #figure <- poly6_plot
  #print(figure)
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 6th", ".pdf"), poly6_plot, width = 8, height = 5, dpi = 450)
  
  
  
  #---------------find the derivative of my fitted polynomial
  
  #first derivative of my model - this gives us the velocity
  xmD1 <- 0:max(df_front$t)
  
  D1 = cf0[1] + 2*cf0[2]*xmD1 + 3*cf0[3]*xmD1^2 + 4*cf0[4]*xmD1^3 + 5*cf0[5]*xmD1^4 + 6*cf0[6]*xmD1^5
  
  model6tD1_df <- data.frame(cbind(xmD1,D1))
  
  poly6_D1_front <- ggplot()+
    geom_line(data=model6tD1_df, mapping=aes(x=xmD1, y=D1),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
    coord_cartesian(xlim=c(0,300), ylim=c(0,20)) + 
    xlab("Time [s]") + 
    ylab("Velocity [\U003BCm/s]") + 
    ggtitle (paste("Velocity of the front over time", '\n', GEN, "series", PLANT, '\n', CON)) + 
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none")
  
  
  poly6_D1_front
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "derivative poly 6th_front", ".png"), poly6_D1_front, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "derivative poly 6th", ".svg"),width=8,height=5)
  #figure <- poly6_plot
  #print(figure)
  #
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "derivative poly 6th", ".pdf"), poly6_D1, width = 8, height = 5, dpi = 450)
  
  
  #this is important beacuse R normally does not remember what happens inside a function, but I can chose one thing to remember as output of a function
  return(coefficient_poly6)
}


#-------------------------------FUNCTION 6
#Fits 4th grade polynomial on the average distance of the max from the centre over time
#calculates derivatives
#plots all

Fit_poly_4 <- function (df_maxima, my_pal_gray, base_plot) {
  
  model_poly4 <- lm(df_maxima$distance_max ~ poly(df_maxima$t,4, raw = T) -1) # the -1 forces the poly to pass from 0,0
  cf0 = coef(model_poly4)
  cf0
  
  xm <- 0:max(df_maxima$t)
  y0 <- cf0[1]*xm + cf0[2]*(xm^2) + cf0[3]*(xm^3) + cf0[4]*(xm^4)
  model_poly4_df <- data.frame(cbind(xm,y0))
  
  
  coefficient_poly <- data.frame(genotype=GEN, sample_name=PLANT, condition=CON)
  
  coefficient_poly$coeff1 <- cf0[1]
  coefficient_poly$coeff2 <- cf0[2]
  coefficient_poly$coeff3 <- cf0[3]
  coefficient_poly$coeff4 <- cf0[4]
  
  coefficient_poly4 <- coefficient_poly
  
  coefficient_poly4
  
  #now plot the model line on top and then save
  
  poly4_plot <- base_plot +
    geom_line(data=model_poly4_df, mapping=aes(x=model_poly4_df$xm, y=model_poly4_df$y0),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  poly4_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 4th", ".png"), poly4_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  svglite(file=paste(GEN, "series", PLANT, "polynomial constrained fit 4th", ".svg"),width=8,height=5)
  figure <- poly4_plot
  print(figure)
  dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  ggsave(paste(GEN, "series", PLANT, "polynomial constrained fit 4th", ".pdf"), poly4_plot, width = 8, height = 5, dpi = 450)
  
  
  
  #---------------find the derivative of my fitted polynomial
  
  #first derivative of my model - this gives us the velocity
  xmD1 <- 0:max(df_maxima$t)
  
  D1 = cf0[1] + 2*cf0[2]*xmD1 + 3*cf0[3]*xmD1^2 + 4*cf0[4]*xmD1^3
  
  model4tD1_df <- data.frame(cbind(xmD1,D1))
  
  poly4_D1 <- ggplot()+
    geom_line(data=model4tD1_df, mapping=aes(x=xmD1, y=D1),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
    coord_cartesian(xlim=c(0,300), ylim=c(0,20)) + 
    xlab("Time [s]") + 
    ylab("Velocity [\U003BCm/s]") + 
    ggtitle (paste("Velocity of the wave over time", '\n', GEN, "series", PLANT, '\n', CON)) + 
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none")
  
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "derivative poly 4th", ".png"), poly4_D1, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  svglite(file=paste(GEN, "series", PLANT, "derivative poly 4th", ".svg"),width=8,height=5)
  figure <- poly4_plot
  print(figure)
  dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  ggsave(paste(GEN, "series", PLANT, "derivative poly 4th", ".pdf"), poly4_D1, width = 8, height = 5, dpi = 450)
  
  
  #this is important beacuse R normally does not remember what happens inside a function, but I can chose one thing to remember as output of a function
  return(coefficient_poly4)
}




#-------------------------------FUNCTION 6
#Fits negative exponential on the average distance of the max from the centre over time
#calculates derivatives
#plots all

Fit_negative_exponential <- function (df_maxima, my_pal_gray, base_plot) {
  
  y <- df_maxima$distance_max
  t <- df_maxima$t
  
  plateau_start <- mean(tail(y,50))
  
  halft_start <- df_maxima[which.min(abs(y - plateau_start/2)), "t"]
  
  model_negexp<-NULL
  
  try(model_negexp <-nls(y ~ alfa * (1 - exp(-beta * t)), 
                         df_maxima, 
                         start=list(alfa = plateau_start , beta = log(2)/halft_start))
  )# does not stop in the case of error
  
  if(!is.null(model_negexp))break
  
  coeff_df <- data.frame(coefficients(model_negexp))
  
  coefficient_negexp <- data.frame(genotype=GEN, sample_name=PLANT, condition=CON)
  coefficient_negexp$alfa <- coeff_df[1,1]
  coefficient_negexp$beta <- coeff_df[2,1]
  
  half_time <- log(2)/coefficient_negexp$beta
  
  coefficient_negexp$half_time <- half_time
  
  coefficient_negexp
  
  x_mod <- c(0:max(t))
  y_mod <-predict(model_negexp, newdata=x_mod)
  
  y_mod <- (coefficient_negexp$alfa * (1 - exp(-coefficient_negexp$beta * x_mod)))
  
  model_negexp_df <- data.frame(cbind(x_mod, y_mod))
  
  
  negexp_plot <- maxd_plot +
    geom_line(data=model_negexp_df, mapping=aes(x=model_negexp_df$x_mod, y=model_negexp_df$y_mod), 
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3)
  
  negexp_plot
  
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "Negative exponential fit", ".png"), negexp_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "Negative exponential fit", ".svg"),width=8,height=5)
  #figure <- negexp_plot
  #print(figure)
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "Negative exponential fit", ".pdf"), negexp_plot, width = 8, height = 5, dpi = 450)
  
  
  #---------------find the derivative of neg exponential 
  #https://mathinsight.org/exploring_derivative_exponential_function
  #https://www.wolframalpha.com/input/?i=4-4e%5E%28-2x%29
  
  #first derivative of my model - this gives us the velocity
  xmD1 <- 0:max(df_maxima$t)
  
  #D1 = alfa*beta*exp(-beta*x)
  D1 = (coefficient_negexp$alfa)*(coefficient_negexp$beta)*exp(-(coefficient_negexp$beta)*xmD1)
  
  modelnegexpD1_df <- data.frame(cbind(xmD1,D1))
  
  negexp_D1 <- ggplot()+
    geom_line(data=modelnegexpD1_df, mapping=aes(x=xmD1, y=D1),
              color=my_pal_gray[9], lty=1, lwd=1.5, alpha=0.3) +
    coord_cartesian(xlim=c(0,300), ylim=c(0,20)) + 
    xlab("Time [s]") + 
    ylab("Velocity [\U003BCm/s]") + 
    ggtitle (paste("Velocity of the wave over time", '\n', GEN, "series", PLANT, '\n', CON)) + 
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none")
  
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "derivative negative exponential", ".png"), negexp_D1, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  #svglite(file=paste(GEN, "series", PLANT, "derivative poly 4th", ".svg"),width=8,height=5)
  #figure <- negexp_D1
  #print(figure)
  #dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  #ggsave(paste(GEN, "series", PLANT, "derivative negative exponential", ".pdf"), negexp_D1, width = 8, height = 5, dpi = 450)
  
  
  #this is important beacuse R normally does not remember what happens inside a function, but I can chose one thing to remember as output of a function
  return(coefficient_negexp)
}



##----------------------------------------------------------------------------------------------------------
##-------------------------------------------- end of functions ---------------------------------------------------
##-------------------------------------------------------------------------------------------------------



#---------loops that run on all folders on all plants on all data points

#set wd to experiment folder
setwd(experiment_folder)

#create workbooks that you will need for the coefficients

#workbook coefficients polinomial 6th grade for peak
workbook_coeff_poly6 <- createWorkbook(paste(EXP, "coeff_6", ".xlsx"))

#workbook coefficients polinomial 6th grade for front
workbook_coeff_poly6_front <- createWorkbook(paste(EXP, "coeff_6_front", ".xlsx"))

#workbook coefficients polinomial 4th grade
workbook_coeff_poly4 <- createWorkbook(paste(EXP, "coeff_4", ".xlsx"))

#workbook coefficients polinomial negative exponential
workbook_coeff_negexp <- createWorkbook(paste(EXP, "coeff_negexp", ".xlsx"))

#workbook for information about the droplet size and delay
drop_info_book <- createWorkbook(paste(EXP, "drop_info", ".xlsx"))

peaks <- list(x_max=numeric(), i=numeric(), y.hat=numeric(), y_values=numeric(), baseline=numeric(), base_sd=numeric(), x_front=numeric(), y_front=numeric(), front_ind=numeric())


folders <- c(genotype1, genotype2)


for (f in folders) {
  
  setwd(f)
  
  list.files()
  
  if (f == genotype1) {
    GEN <- GEN1}
  if (f == genotype2) {
    GEN <- GEN2}

  
  #NAME <- paste(GEN," ", "series"," ", PLANT, sep="")
  #check name matches the file name in the directory
  #NAME
  
  #creates a workbook where to save the position of the maxima over time
  workbook_loess_max <- createWorkbook(paste(GEN, "maxd_loess", ".xlsx"))
  
  #creates a workbook where to save the position of the front over time
  workbook_loess_front <- createWorkbook(paste(GEN, "frontd_loess", ".xlsx"))
  
  #start_data is a excel file where for each plant of this genotype I store information on which is the wounding frame (wf) and which is the frame where the leaf is in focus and stabe again (sf)
  feed_fun <- data.frame(read.xlsx("start_data.xlsx", 1))
  
  #this creates an excel file for the genotype where I will save the data of the average profile fo the wave over time
  workbook_profiles <- createWorkbook(paste(GEN, "profiles", ".xlsx"))
  
  #profile value according to loess fitted curve
  workbook_profiles_loess <- createWorkbook(paste(GEN, "profiles_loess", ".xlsx"))
  
  # this creates an empty df that I will progressively populate with data of each replicate
  coeff_poly6_total=NULL
  coeff_poly6_total <- data.frame(genotype=character(), sample_name=character(), condition=character(), coeff1=numeric(), coeff2=numeric(), coeff3=numeric(), coeff4=numeric(), coeff5=numeric(), coeff6=numeric())
  
  
  coeff_poly6_total_front=NULL
  coeff_poly6_total_front <- data.frame(genotype=character(), sample_name=character(), condition=character(), coeff1=numeric(), coeff2=numeric(), coeff3=numeric(), coeff4=numeric(), coeff5=numeric(), coeff6=numeric())
  
  
  Drop_total=NULL
  Drop_total <- data.frame(delay=numeric(), r_drop=numeric(), vol_ul=numeric())
  
  
  #x11()
  
  
  #loop along all the files in the folder
  plants <- feed_fun[,"pl"]
  
  for (k in plants) { #with this small loop I feed my function
    
    #reset this to the edge of my window, this guarantees that:
    #- the detection of the peak in the following plant is not going to be biased by the position of the last max in the previous plant
    #- if at the start of the data I can't detect any peak, default goes to the minimum distance possible
    peaks$i <- (w+1)
    peaks$front_ind <- (w+1)
    
    df_maxima=NULL
    max_j=NULL
    df_front=NULL
    front_j=NULL
    df_maxima <- data.frame(t=numeric(), distance_max=numeric(), intensity_max=numeric(), s_n_max=numeric())
    df_front <- data.frame(t=numeric(), distance_front=numeric(), intensity_front=numeric(), baseline=numeric(), base_sd=numeric())
    
    model_j=NULL
    model_profile=NULL
    model_profile <- data.frame(t=numeric(), space=numeric(), profile_smooth=numeric())
    
    PLANT <- feed_fun[which(feed_fun$pl==k), "pl"]
    woundframe <- feed_fun[which(feed_fun$pl==k), "wf"]
    stableframe <- feed_fun[which(feed_fun$pl==k), "sf"]
    
    NAME <- paste(GEN," ", "series", " ", PLANT,".csv", sep="")
    
    #calculates the profile of the wave for the whole dataset
    avg_profile <- wave_profile(PLANT, GEN)
    avg_profile <- data.frame(avg_profile)
    
    #calculates noise for the whole dataset
    avg_profile_sn <- estimate_noise(treshold, avg_profile)

    
    # create a directory for the plots that you will generate
    plotdir <- file.path(getwd(), paste(GEN, " ", 'series', " ", PLANT, "profile_loess6"))
    dir.create(plotdir)
    
    #times<-as.numeric(levels(as.factor(avg_profile_sn$t)))
    times <- unique(avg_profile_sn$t)
    
    #loop that applies the two functions on all time points of your replicate
    for (j in times) {
      
      profile_t <- avg_profile_sn[avg_profile_sn$t==j,]
      y <- profile_t$avrg
      
      if (var(y)==0) {
        print(paste(NAME, "skipped t=", j,"for var=0" ))
        next
      }
      x <- profile_t$d
      SE <- profile_t$se
      z <- profile_t$s_n
          
      
      #newline
      baseline <- profile_t$baseline
      base_sd <- profile_t$avg_sd
    
    
      if (j <= 150.00) {
        fit_loess_peaks (w, span1, j, SE, x, profile_t)
        
        peaks <- d_max(x, y, w=w, span=span1)
      }
      else {
        fit_loess_peaks (w, span2, j, SE, x, profile_t)
        
        peaks <- d_max(x, y, w=w, span=span2)
      }
      
      max_j <- cbind(t=j, distance_max=peaks$x_max, intesity_max=peaks$y_values, s_n_max=profile_t[which(profile_t$d==peaks$x_max), "s_n"])
      front_j <- cbind(t=j, distance_front=peaks$x_front, intensity_front=peaks$y_front, baseline=peaks$baseline, base_sd=peaks$base_sd)
     
       df_maxima <- rbind(df_maxima, max_j)
       df_front <- rbind(df_front, front_j)
      
      model_j <- data.frame(cbind(t=j, space=x, profile_smooth=peaks$y.hat))
      model_profile <- rbind(model_profile, model_j)
      
    }
    
    
    #sometimes I get NA in the front, I am not sure why at the moment but they prevent analyis to proceed, so here I clean them up
    df_front <- df_front[complete.cases(df_front),]
    df_maxima <- df_maxima[complete.cases(df_front),]
    
    
    addWorksheet(workbook_loess_max, NAME)
    writeData(workbook_loess_max, NAME , df_maxima)
    saveWorkbook(workbook_loess_max, file=paste(GEN, "maxd_loess", ".xlsx"), overwrite = TRUE)
    
    addWorksheet(workbook_loess_front, NAME)
    writeData(workbook_loess_front, NAME , df_front)
    saveWorkbook(workbook_loess_front, file=paste(GEN, "frontd_loess", ".xlsx"), overwrite = TRUE)
    
    addWorksheet(workbook_profiles_loess, NAME)
    writeData(workbook_profiles_loess, NAME , model_profile)
    saveWorkbook(workbook_profiles_loess, file=paste(GEN, "profiles_loess", ".xlsx"), overwrite = TRUE)    
    
    ##---------------------------------------------------------------------------------------------------------    
    #----------------basic plot of position of the max from the centre over time
    
    df_maxima1<- df_maxima %>% mutate(Group = cut(s_n_max, breaks = c(0,1.06,1.5,2, Inf), include.lowest = TRUE))
    
    maxd_plot <- ggplot()+
      geom_point(data=df_maxima1, aes(x=t, y=distance_max, fill=Group), shape=21, col=my_pal_gray[5], stroke=0.01, size=3.5, alpha=1)+
      xlab("Time [s]") +
      ylab(paste("Distance from the centre", "\n" , "[\U003BCm]"))+
      theme_bw(base_size=18)+
      theme(plot.title = element_text(hjust = 0.5))+
      ggtitle(paste("Position of max over time", '\n', GEN, "series", PLANT, '\n', CON, sep=""))+
      labs(fill="signal/noise") +
      scale_fill_manual(breaks = levels(df_maxima1$Group), drop = FALSE, values = my_pal_quant_2[c(1,3,5,8)],
                        labels=c(expression("s/n" <= "1"), "1 < s/n < 1.5", expression("1.5"<="s/n < 2"), expression("s/n" >= "2")))+
      guides(aesthetics = "fill", fill = guide_legend(reverse = TRUE, override.aes = list(shape = 21, size= 10)))
    
    maxd_plot
    
    #saves PNG version for easy acces
    ggsave(paste(GEN, "series", PLANT, "distance max", ".png"), maxd_plot, width = 8, height = 5, dpi = 450, device = 'png')
    
    #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
    #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
    svglite(file=paste(GEN, "series", PLANT, "distance max", ".svg"),width=8,height=5)
    figure <- maxd_plot
    print(figure)
    dev.off()
    
    #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
    #illustrator, however this renders correctly in inkscape
    ggsave(paste(GEN, "series", PLANT, "distance max", ".pdf"), maxd_plot, width = 8, height = 5, dpi = 450)
    
    
    #----- 6th grade polynomial fitting and delay decision
    
    #unlike the wounding datasets, the elicitors dataset presents a delay between the
    #deposition of the droplet and the elicitation/start of a wave, also sometimes there is no signal at 
    #all being elicited from the droplet.
    
    #first what happens if there is no wave and the max stays where it is - if there is no progress at all, I don't even go searching for the delay
    
    if (var(df_maxima$distance_max, na.rm = TRUE)==0) {
      print(paste(NAME, "no wave, dist max has var=0" ))
      
      delay <- 0
      
      result <- Fit_poly_6_max (df_maxima, my_pal_gray, base_plot=maxd_plot)
      
      coeff_poly6_total <- rbind(coeff_poly6_total, result)
    }
    
    # Then if the var is not 0, we procced to fit annotate the delay and fit a curve. I think I is more precise 
    # to fit a curve that starts when the wave starts, so that the delay in the elicitation 
    # does not affect the calciulation of the speed of the wave.
    # first thing to do is to check when the wavefront starts moving:
    
    else {
      
      start <- Position(function(x) x > min(df_maxima$distance_max), df_maxima$distance_max)
      
      start
      
      #a start of 2 means that the wave did not sit still for even one single frame, so it means no
      #delay. However a start higher than 2 means the wave did not kickstart immediately and so it has a
      #delay.
      
      #if delay is > 160 sec, probably it is not a delay that preceeds a real wave, so I may as well proceed with no trimming...
      
      if (start >= 80) {
        
        delay <- 0
        
        result <- Fit_poly_6_max (df_maxima, my_pal_gray, base_plot=maxd_plot)
        
        coeff_poly6_total <- rbind(coeff_poly6_total, result)
        
      } 
      
      
      #while, if the start is earlier that 160 sec...
      
      if (start < 80) {
        
        #...but there is a small delay, ie not immediate start, then I want to proceed with trimming
        
        if (start > 2) {
          
          delay <- df_maxima$t[start]
          
          #cut away anything that is before the wave kicks off
          dat_trim <- subset(df_maxima, t>=delay)
          head(dat_trim, 10)
          
          #sets the wonding time as zero
          stime <- min(dat_trim[,"t"])
          dat_trimt0 <- transform(dat_trim, t = t-stime)
          
          df_maxima2 <- data.frame(dat_trimt0)
          
          df_maxima2<- df_maxima2 %>% mutate(Group = cut(s_n_max, breaks = c(0,1.06,1.5,2, Inf), include.lowest = TRUE))
          
          maxd_plot_trim <- ggplot()+
            geom_point(data=df_maxima2, aes(x=t, y=distance_max, fill=Group), shape=21, col=my_pal_gray[5], stroke=0.01, size=3.5, alpha=1)+
            xlab("Time [s]") +
            ylab(paste("Distance from the centre", "\n" , "[\U003BCm]"))+
            theme_bw(base_size=18)+
            theme(plot.title = element_text(hjust = 0.5))+
            ggtitle(paste("Position of max over time", '\n', GEN, "series", PLANT, '\n', CON, sep=""))+
            labs(fill="signal/noise") +
            scale_fill_manual(breaks = levels(df_maxima2$Group), drop = FALSE, values = my_pal_quant_2[c(1,3,5,8)],
                              labels=c(expression("s/n" <= "1"), "1 < s/n < 1.5", expression("1.5"<="s/n < 2"), expression("s/n" >= "2")))+
            guides(aesthetics = "fill", fill = guide_legend(reverse = TRUE, override.aes = list(shape = 21, size= 10)))
          
          maxd_plot_trim
          
          #saves PNG version for easy acces
          ggsave(paste(GEN, "series", PLANT, "distance max_trim", ".png"), maxd_plot_trim, width = 8, height = 5, dpi = 450, device = 'png')
          
          #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
          #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
          svglite(file=paste(GEN, "series", PLANT, "distance max_trim", ".svg"),width=8,height=5)
          figure <- maxd_plot_trim
          print(figure)
          dev.off()
          
          #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
          #illustrator, however this renders correctly in inkscape
          ggsave(paste(GEN, "series", PLANT, "distance max_trim", ".pdf"), maxd_plot_trim, width = 8, height = 5, dpi = 450)
          
          
          #fit curves on the trimmed data so that the start of the curve is at the start of the wave
          result <- Fit_poly_6_max (df_maxima2, my_pal_gray, base_plot=maxd_plot_trim)
          
          coeff_poly6_total <- rbind(coeff_poly6_total, result)
          
          
          #----- 4th grade polynomial    
          
          #result <- Fit_poly_4 (df_maxima, my_pal_gray, base_plot=maxd_plot) 
          
          #coeff_poly4_total <- rbind(coeff_poly4_total, result)
          
          
          #----- negative exponential
          #result <- Fit_negative_exponential (df_maxima, my_pal_gray, base_plot=maxd_plot) 
          
          #coeff_negexp_total <- rbind(coeff_negexp_total, result)
          
          
          
          #if a trim has been necessary, I want to save the trimmed df max
          NICK <- gsub('.{4}$', '', NAME)
          
          if (nchar(NICK)>26){
            NICK <- gsub("[[:blank:]]", "", NICK)
          }
          
          addWorksheet(workbook_loess_max, paste('trim',NICK))
          writeData(workbook_loess_max, paste('trim',NICK) , df_maxima2)
          saveWorkbook(workbook_loess_max, file=paste(GEN, "maxd_loess", ".xlsx"), overwrite = TRUE)
          
          
        }
        
        #...whereas if  start is <50 but not >2, we say that this wave had no delay and the original untrimmed data are used to fit the curves on
        
        else {delay <- 0
        
        result <- Fit_poly_6_max (df_maxima, my_pal_gray, base_plot=maxd_plot)
        
        coeff_poly6_total <- rbind(coeff_poly6_total, result)
        
        
        
        #----- 4th grade polynomial    
        
        #result <- Fit_poly_4 (df_maxima, my_pal_gray, base_plot=maxd_plot) 
        
        #coeff_poly4_total <- rbind(coeff_poly4_total, result)
        
        
        #----- negative exponential
        #result <- Fit_negative_exponential (df_maxima, my_pal_gray, base_plot=maxd_plot) 
        
        #coeff_negexp_total <- rbind(coeff_negexp_total, result)
        
        }
        
      }
      
    } #this closes the delay or not delay decision, then, before proceeding to next plant, I want to save the delay and the size of the droplet
    
    r_drop <- avg_profile_sn[2,'r_drop']
    
    Drop_k = NULL
    Drop_k <- data.frame(cbind(delay, r_drop))
    
    #approximate the volume to the volume of half sphere
    vol_cub_um <- 2*(pi*((r_drop)^3))/3
    
    vol_ul <- vol_cub_um/(10^9)
    
    Drop_k <- data.frame(cbind(Drop_k, vol_ul))
    
    
  ##---------------------------------------------------------------------------------------------------------    
  #----------------basic plot of position of the front from the centre over time
  

  
  frontd_plot <- ggplot()+
    geom_point(data=df_front, aes(x=t, y=distance_front), shape=21, col=my_pal_gray[5], stroke=0.01, size=3.5, fill=okabe_pal[3], alpha=0.4)+
    xlab("Time [s]") +
    ylab(paste("Distance from the centre", "\n" , "[\U003BCm]"))+
    theme_bw(base_size=18)+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle(paste("Position of front over time", '\n', GEN, "series", PLANT, '\n', CON, sep=""))

  frontd_plot
  
  #saves PNG version for easy acces
  ggsave(paste(GEN, "series", PLANT, "distance front", ".png"), frontd_plot, width = 8, height = 5, dpi = 450, device = 'png')
  
  #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
  #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
  svglite(file=paste(GEN, "series", PLANT, "distance front", ".svg"),width=8,height=5)
  figure <- frontd_plot
  print(figure)
  dev.off()
  
  #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
  #illustrator, however this renders correctly in inkscape
  ggsave(paste(GEN, "series", PLANT, "distance front", ".pdf"), frontd_plot, width = 8, height = 5, dpi = 450)
  
  
  
  #----- 6th grade polynomial fitting and delay decision
  
  #unlike the wounding datasets, the elicitors dataset presents a delay between the
  #deposition of the droplet and the elicitation/start of a wave, also sometimes there is no signal at 
  #all being elicited from the droplet.
  
  #first what happens if there is no wave and the max stays where it is - if there is no progress at all, I don't even go searching for the delay
  
  if (var(df_front$distance_front, na.rm = TRUE)==0) {
    print(paste(NAME, "no wave, dist front has var=0" ))
    
    delay <- 0
    
    result <- Fit_poly_6_front (df_front, my_pal_gray, base_plot=frontd_plot)
    
    coeff_poly6_total_front <- rbind(coeff_poly6_total_front, result)
  }
  
  # Then if the var is not 0, we procced to fit annotate the delay and fit a curve. I think I is more precise 
  # to fit a curve that starts when the wave starts, so that the delay in the elicitation 
  # does not affect the calciulation of the speed of the wave.
  # first thing to do is to check when the wavefront starts moving:
  
  else {
    
    start <- Position(function(x) x > min(df_front$distance_front), df_front$distance_front)
    
    start
    
    #a start of 2 means that the wave did not sit still for even one single frame, so it means no
    #delay. However a start higher than 2 means the wave did not kickstart immediately and so it has a
    #delay.
    
    #if delay is > 100 sec, probably it is not a delay that preceeds a real wave, so I may as well proceed with no trimming...
    
    if (start >= 50) {
      
      delay <- 0
      
      result <- Fit_poly_6_front (df_front, my_pal_gray, base_plot=frontd_plot)
      
      coeff_poly6_total_front <- rbind(coeff_poly6_total_front, result)
      
    } 
    
    
    #while, if the start is earlier that 100 sec...
    
    if (start < 50) {
      
      #...but there is a small delay, ie not immediate start, then I want to proceed with trimming
      
      if (start > 2) {
        
        delay <- df_front$t[start]
        #cut away anything that is before the wave kicks off
        dat_trim <- subset(df_front, t>=delay)
        head(dat_trim, 10)
        
        #sets the wonding time as zero
        stime <- min(dat_trim[,"t"])
        dat_trimt0 <- transform(dat_trim, t = t-stime)
        
        df_front <- data.frame(dat_trimt0)
        
        
       frontd_plot_trim <- ggplot()+
          geom_point(data=df_front, aes(x=t, y=distance_front), shape=21, col=my_pal_gray[5], fill=okabe_pal[3], stroke=0.01, size=3.5, alpha=0.4)+
          xlab("Time [s]") +
          ylab(paste("Distance from the centre", "\n" , "[\U003BCm]"))+
          theme_bw(base_size=18)+
          theme(plot.title = element_text(hjust = 0.5))+
          ggtitle(paste("Position of front over time", '\n', GEN, "series", PLANT, '\n', CON, sep=""))
          
        frontd_plot_trim
        
        #saves PNG version for easy acces
        ggsave(paste(GEN, "series", PLANT, "distance front_trim", ".png"), frontd_plot_trim, width = 8, height = 5, dpi = 450, device = 'png')
        
        #saves svg version for easy editing, this is correctly visualized only in illustrator while in inkscape looses the 
        #strokes and looses size of the legend. To insert in thesis, open in illustrator, save and save as pdf
        svglite(file=paste(GEN, "series", PLANT, "distance front_trim", ".svg"),width=8,height=5)
        figure <- frontd_plot_trim
        print(figure)
        dev.off()
        
        #saves pdf version for easy insertion in thesis. Problem with this pdf is that the data point are rendered as squares in 
        #illustrator, however this renders correctly in inkscape
        ggsave(paste(GEN, "series", PLANT, "distance front_trim", ".pdf"), frontd_plot_trim, width = 8, height = 5, dpi = 450)
        
        
        #fit curves on the trimmed data so that the start of the curve is at the start of the wave
        result <- Fit_poly_6_front (df_front, my_pal_gray, base_plot=frontd_plot_trim)
        
        coeff_poly6_total_front <- rbind(coeff_poly6_total_front, result)
        
        
        #if a trim has been necessary, I want to save the trimmed df max
        NICK <- gsub('.{4}$', '', NAME)
        
        if (nchar(NICK)>26){
          NICK <- gsub("[[:blank:]]", "", NICK)
        }
        
        addWorksheet(workbook_loess_front, paste('trim',NICK))
        writeData(workbook_loess_front, paste('trim',NICK) , df_front)
        saveWorkbook(workbook_loess_front, file=paste(GEN, "frontd_loess", ".xlsx"), overwrite = TRUE)
        
        
      }
      
      #...whereas if  start is <50 but not >2, we say that this wave had no delay and the original untrimmed data are used to fit the curves on
      
      else {delay <- 0
      
      result <- Fit_poly_6_front (df_front, my_pal_gray, base_plot=frontd_plot)
      
      coeff_poly6_total_front <- rbind(coeff_poly6_total_front, result)
      
      
      
      #----- 4th grade polynomial    
      
      #result <- Fit_poly_4 (df_maxima, my_pal_gray, base_plot=maxd_plot) 
      
      #coeff_poly4_total <- rbind(coeff_poly4_total, result)
      
      
      #----- negative exponential
      #result <- Fit_negative_exponential (df_maxima, my_pal_gray, base_plot=maxd_plot) 
      
      #coeff_negexp_total <- rbind(coeff_negexp_total, result)
      
      }
      
    }
    
  } #this closes the delay or not delay decision on the front 
  
  
  delay_front <- delay
  
  Drop_k <- data.frame(cbind(Drop_k, delay_front))
  
  Drop_total <- rbind(Drop_total, Drop_k)
  
    
  } #this closes for k in plants
  
  
  #before proceeding to the next genotype, I need to save the polynomials and the drop info
  
  setwd(experiment_folder)
  
  addWorksheet(drop_info_book, GEN)
  writeData(drop_info_book, GEN , Drop_total)
  saveWorkbook(drop_info_book, file=(paste(EXP, "drop_info", ".xlsx")), overwrite = TRUE)
  
  addWorksheet(workbook_coeff_poly6, GEN)
  writeData(workbook_coeff_poly6, GEN , coeff_poly6_total)
  saveWorkbook(workbook_coeff_poly6, file=(paste(EXP, "coeff_6", ".xlsx")), overwrite = TRUE)
  
  addWorksheet(workbook_coeff_poly6_front, GEN)
  writeData(workbook_coeff_poly6_front, GEN , coeff_poly6_total_front)
  saveWorkbook(workbook_coeff_poly6_front, file=(paste(EXP, "coeff_6_front", ".xlsx")), overwrite = TRUE)
  
  
} #this closes for f in folders

