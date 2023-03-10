## clearing the environment
rm(list = ls())  
gc()    

library(tidyverse)
####################################################################################
# importing ABM output'
mod_out1 <- read.csv('mytest_output.csv')
mod_out2 <- read.csv('myparallel_output.csv')

ggplot(mod_out1) +
  labs(x="Host age (days)", y=NULL, title = "Normalized donor fraction") +
  ylim(0, 1) +
  geom_line(aes(x=time, y=Normalized_fd), col='navy', size=1.05)


ggplot(mod_out1) +
  geom_line(aes(x=time, y=physiol_counts), col="navy", size=1.05) +
  labs(x="Host age (days)", y=NULL, title = "Cell counts") +
  scale_y_log10(limits=c(1e4, 1e7))



ki_df <- mod_out1 %>%
  select(time, contains("Ki67")) %>%
  na_if(0)
  
ggplot(ki_df) +
  geom_line(aes(x=time, y=Donor_Ki67_pos), col=2) + 
  geom_line(aes(x=time, y=Host_Ki67_pos), col=4) + 
  labs(x="Host age (days)", y=NULL, title = "Fraction of Ki67 positive") +
  ylim(0, 1)
  
  





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