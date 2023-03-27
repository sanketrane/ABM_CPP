## clearing the environment
rm(list = ls())  
gc()    

library(tidyverse)
s
####################################################################################
myTheme <- theme(text = element_text(size = 12), axis.text = element_text(size = 12), axis.title =  element_text(size = 12, face = "bold"),
                 plot.title = element_text(size=12, face = 'bold',  hjust = 0.5), legend.text = element_text(size=12),
                 legend.title = element_text(size = 12, face = 'bold'), legend.background = element_blank())

# setting ggplot theme for rest fo the plots
theme_set(theme_bw())

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

log10minorbreaks=as.numeric(1:10 %o% 10^(4:8))


####################################################################################
## loading required datasets for plotting
counts_file <- file.path("data", "Counts_naiTreg.csv")
counts_data <- read.csv(counts_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4')))) %>%
  gather(c(Thymus, Periphery), key='location', value = 'total_counts')

Nfd_file <- file.path("data", "Nfd_naiTreg.csv")
Nfd_data <- read.csv(Nfd_file) %>% 
  arrange(age.at.S1K)%>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))))%>%
  gather(c(Thymus, Periphery), key='location', value = 'Nfd')

hostki_file <- file.path("data", "hostKi67_naiTreg.csv")
hostki_data <- read.csv(hostki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Host') %>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

donorki_file <- file.path("data", "donorKi67_naiTreg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         subcomp='Donor')%>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

ki_data <- rbind(donorki_data, hostki_data)

####################################################################################
# importing ABM output'
mod_out1 <- read.csv('mytest_output.csv')
#mod_out2 <- read.csv('myparallel_output.csv')

legn_labels <- c('6-8', '8-10', '10-12', '12-25')

ggplot(mod_out1) +
  geom_line(aes(x=time, y=Normalized_fd), col='navy', size=1.05)+
  geom_point(data = Nfd_data, aes(x = age.at.S1K, y = Nfd, color = ageBMT_bin), size=2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(1, 350), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  guides(fill='none')+ myTheme


ggplot(mod_out1) +
  geom_line(aes(x=time, y=physiol_counts), col="navy", size=1.05) +
  geom_point(data = counts_data, aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(60, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(5e3, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  guides(fill = 'none') + myTheme 

ki_df <- mod_out1 %>%
  select(time, contains("Ki67")) %>%
  na_if(0)
  
# thymic Ki67 fractions

fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')

ggplot(ki_df) +
  geom_line(aes(x=time, y=Donor_Ki67_pos*100), col=2) + 
  geom_line(aes(x=time, y=Host_Ki67_pos*100), col=4) + 
  geom_point(data = ki_data, aes(x = age.at.S1K, y = prop_ki*100, color = subcomp), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in thymic naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 50), breaks = c(0, 10, 20, 30, 40, 50))+ 
  guides(fill='none') + myTheme + theme(legend.title = element_blank())




