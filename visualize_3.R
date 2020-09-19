  library(tidyverse)
  library(ggplot2)
  
  
  #Bekannte TADS einlesen
  TAD_BED_bal<-read.delim("log2_O_E/bonev_mESC_MAPQ30_KR.balanced.chr1.bed",header = F)
  TAD_BED_imbal<-read.delim("log2_O_E/bonev_mESC_MAPQ30_KR.imbalanced.chr1.bed",header = F)
  
  
  dateipfad <- ("../../Aufgabe 1/chr1_10kb_mESC_mapq30_KR_10Mb_diag.txt")
  
  for (i in (1:30)){
    TAD_One_bal<-TAD_BED_bal[,2]
    TAD_Two_bal <- TAD_BED_bal[,3]
    T1<-TAD_One_bal[i]
    T2<- TAD_Two_bal[i]
    dif <- abs(T2-T1)
    T1 <- T1-dif
    T2 <- T2 + dif
    dummy = read.delim(dateipfad, header= FALSE, col.names = c("start_one","start_two","cont")) %>%
      filter(start_one >= T1 & start_one <= T2 & start_two >= T1 & start_two <= T2)
    NoNa <- dummy %>% drop_na()
    maxCont <-max(NoNa[,3])
    midCont <- median(NoNa[,3])
    qu95<-as.numeric(quantile(NoNa[,3],probs=0.95))
    bin_size <- 10000
    
    
    NoNa %>%
      ggplot(aes(start_one,start_two)) +
      geom_bin2d(bins = bin_size, binwidth = c(bin_size/2, bin_size/2), aes(fill = after_stat(NoNa[,3]))) +
      #scale_fill_gradientn(colours = heat.colors(n = 2, alpha = 1, rev = TRUE)) +
      scale_fill_gradient2(low = "white", mid = "yellow", high="darkred", limit = c(0,maxCont), space = "Lab", na.value = "black") +
      expand_limits(x = range(NoNa[,1]), y = range(NoNa[,2])) +
      theme_bw() +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey")) +
      coord_fixed()
    ggsave(filename = paste(T1,"_", T2,"HeatMapRAW.png"), limitsize = FALSE)   
    
    
    data_shifted = NoNa %>% 
      mutate(cont = case_when(start_one==start_two ~ NaN, start_one != start_two ~ cont)) %>%
      mutate(cont = case_when(cont >= qu95 ~ qu95, cont < qu95 ~ cont))
    data_shifted %>%    
      ggplot(aes(start_one,start_two)) +
      geom_bin2d(bins = bin_size, binwidth = c(bin_size/2, bin_size/2), aes(fill = after_stat(data_shifted[,3]))) +
      scale_fill_gradientn(colours = heat.colors(n = 2, alpha = 1, rev = TRUE)) +
      #scale_fill_gradient(low = "white", high="darkred", limit = c(0,181), space = "Lab", na.value = "black") +
      expand_limits(x = range(data_shifted[,1]), y = range(data_shifted[,2])) +
      theme_bw() +
      theme(
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "grey")) +
      coord_fixed()
    ggsave(filename = paste(T1,"_", T2,"HeatMapNew.png"), limitsize = FALSE)  
  }  
  
  
  # dtt %>% 
  #   ggplot(aes(V1,V2)) +
  #   geom_bin2d(bins = 100, binwidth = c(0.5, 0.5), aes(fill = after_stat(dtt[,3]))) +
  #   scale_fill_gradientn(colours = heat.colors(n = 2, alpha = 1, rev = TRUE)) +
  #   #scale_fill_gradient(low = "white", high="darkred", limit = c(0,181), space = "Lab", na.value = "black") +
  #   expand_limits(x = range(alles[,1]), y = range(alles[,1])) +
  #   theme_bw() +
  #   theme(
  #     plot.background = element_rect(fill = "white"),
  #     panel.background = element_rect(fill = "white"),
  #     axis.line.x = element_line(color = "grey")) +
  #   coord_fixed()