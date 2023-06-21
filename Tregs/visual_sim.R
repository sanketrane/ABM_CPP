## clearing the environment
rm(list = ls())  
gc()    

library(tidyverse)

### Model name
modelName <- 'Neutral'

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
         popln='Host') %>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

donorki_file <- file.path("data", "donorKi67_naiTreg.csv")
donorki_data <- read.csv(donorki_file) %>% 
  arrange(age.at.S1K) %>%
  mutate(ageBMT_bin = ifelse(age.at.BMT <= 56, 'agebin1',
                             ifelse(age.at.BMT <= 70, 'agebin2',
                                    ifelse(age.at.BMT <= 84, 'agebin3', 'agebin4'))),
         popln='Donor')%>%
  gather(c(Thymus, Periphery), key='location', value = 'prop_ki')

ki_data <- rbind(donorki_data, hostki_data)

####################################################################################
#import and merge all three CSV files into one data frame
df <- list.files(path= file.path("output_csv", modelName, "finished_runs"), full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows %>%
  rename(host_age = time.int, 
         thymic_tregs = physiol_counts,
         fp3neg_SP= sp.numbers)

df_single <- read_csv(file.path("output_csv", modelName, "outfile_m1.csv")) %>%
  rename(host_age = time.int, 
         thymic_tregs = physiol_counts,
         fp3neg_SP= sp.numbers)


counts_df <- df_single %>%
  select("host_age", "thymic_tregs") %>%
  gather(-host_age, key='popln', value='counts') %>% na.omit() %>%
  group_by(host_age) %>%
  summarize(lb = quantile(counts, probs = 0.05),
            median = quantile(counts, probs = 0.5),
            ub = quantile(counts, probs = 0.95))

nfd_df <- df_single %>%
  select("host_age", "Normalized_fd") %>%
  gather(-host_age, key='popln', value='Normalized_fd') %>% na.omit() %>%
  group_by(host_age) %>%
  summarize(lb = quantile(Normalized_fd, probs = 0.05),
            median = quantile(Normalized_fd, probs = 0.5),
            ub = quantile(Normalized_fd, probs = 0.95))

ki_df <- df_single %>%
  select("host_age", "Donor_Ki67_pos", "Host_Ki67_pos") %>% na.omit() %>%
  rename(Donor = Donor_Ki67_pos, Host = Host_Ki67_pos) %>%
  gather(-host_age, key='popln', value='counts') %>%
  group_by(host_age, popln) %>%
  summarize(lb = quantile(counts, probs = 0.05),
            median = quantile(counts, probs = 0.5),
            ub = quantile(counts, probs = 0.95))
  


legn_labels <- c('6-8', '8-10', '10-12', '12-25')

ggplot(nfd_df) +
  geom_line(aes(x=host_age, y=median))+
  geom_ribbon(data = counts_df, aes(x = host_age, ymin = lb, ymax = ub), alpha = 0.15)+
  geom_point(data = filter(Nfd_data, location=="Thymus"), aes(x = age.at.S1K, y = Nfd, color = ageBMT_bin), size=2) +
  labs(x = "Host age (days)", y = NULL, title = "Normalised Chimerism in naive Tregs") +
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(30, 450), breaks = c(0,100,200,300, 400, 500))+
  scale_y_continuous(limits =c(0, 1.02), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  #facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
  guides(fill='none')+ myTheme


ggplot(counts_df) +
  geom_line(aes(x=host_age, y=median)) +
  geom_ribbon(data = counts_df, aes(x = host_age, ymin = lb, ymax = ub), alpha = 0.15)+
  geom_point(data = filter(counts_data, location=="Thymus"), aes(x = age.at.S1K, y = total_counts, color = ageBMT_bin), size=2) +
  labs(title=paste('Total counts of naive Tregs'),  y=NULL, x= "Host age (days)") + 
  scale_color_discrete(name="Host age at \n BMT (Wks)", labels=legn_labels)+
  scale_x_continuous(limits = c(40, 450) , trans="log10", breaks=c(10, 30, 100, 300))+
  scale_y_continuous(limits = c(5e3, 5e6), trans="log10", breaks=c(1e4, 1e5, 1e6, 1e7, 1e8), minor_breaks = log10minorbreaks, labels =fancy_scientific) +
  facet_wrap(~ factor(location, levels = c('Thymus', "Periphery")))+
  guides(fill = 'none') + myTheme 


fac_labels <- c(`agebin1`= '6-8 weeks', `agebin2`= '8-10 weeks', `agebin3`= '10-12 weeks', `agebin4`= '12-25 weeks')

ggplot(ki_df) +
  geom_line(aes(x=host_age, y=median*100, col=popln)) + 
  geom_ribbon(aes(x = host_age, ymin = lb, ymax = ub), alpha = 0.15)+
  geom_point(data = ki_data, aes(x = age.at.S1K, y = prop_ki*100, color = popln), size=1.5) +
  labs(x = "Host age (days)", y = NULL, title = "% Ki67hi in thymic naive Tregs") +
  scale_x_continuous(limits = c(60, 450), breaks = c(0,100,200,300, 400, 500)) +
  scale_y_continuous(limits =c(0, 100), breaks = c(0, 60, 20, 80, 40, 100)) + 
  guides(fill='none') + myTheme + theme(legend.title = element_blank())




