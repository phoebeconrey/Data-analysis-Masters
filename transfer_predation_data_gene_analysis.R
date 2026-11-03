library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)
library(RefSeeker)

ct_data <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/ACT_GAP_LAC_transfer_predation_samples_4_nov_2026.txt",
                      header = TRUE, sep = "\t") 

colnames(ct_data) <- c("Well", "Target", "Content", "Sample", "Cq")

#data_long <- ct_data %>% pivot_longer(cols = -c(Well, Target, Content),
                                      #names_to = "Sample", value_to = "Cq")
ct_clean <- ct_data %>%
  filter(Sample != "NTC", Sample != "NRT") %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Sample, Target) %>%
  summarise(Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  mutate(Cq = ifelse(is.nan(Cq), NA, Cq))

ct_samples_rows <- ct_clean %>%
  pivot_wider(names_from = Target, values_from = Cq) %>%
  select(-Sample)

print(ct_samples_rows)

rs_results <- rs_reffinder(ct_samples_rows)

print(rs_results)
