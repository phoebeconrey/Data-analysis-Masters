library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)
library(RefSeeker)

ct_data <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/ACT_GAP_LAC_transfer_predation_samples_4_nov_2026.txt",
                      header = TRUE, sep = "\t") 

colnames(ct_data) <- c("Group", "Well", "Target", "Content", "Sample", "Cq")

# Check stability, collapse into only genes, then run rs ----
ct_clean <- ct_data %>%
  filter(Sample != "NTC", Sample != "NRT") %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Sample, Target) %>%
  summarise(Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  mutate(Cq = ifelse(is.nan(Cq), NA, Cq))

ct_samples_rows <- ct_clean %>%
  pivot_wider(names_from = Target, values_from = Cq) %>%
  dplyr::select(-Sample)

print(ct_samples_rows)

rs_results <- rs_reffinder(ct_samples_rows)

print(rs_results)

#Delta Ct Analysis ----
delta_data <- read_Ct_long(path = "C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/ACT_GAP_LAC_transfer_predation_samples_4_nov_2026.txt",
                           sep = "\t",
                           dec = ".",
                           skip = 0,
                           add.column.Flag = FALSE,
                           column.Sample = 5,
                           column.Gene = 3,
                           column.Ct = 6,
                           column.Group = 1)
delta_data_clean <- delta_data %>%
  filter(Sample != "NTC", Sample != "NRT")

print(delta_data_clean)

data.CtF.ready <- make_Ct_ready(data = delta_data_clean,
                                imput.by.mean.within.groups = TRUE)

data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         normalise = TRUE,
                         ref = c("Actin", "GAPDH"),
                         transform = TRUE)

results.pn.dCt <- RQ_dCt(data = data.dCt.exp,
                         do.tests = TRUE,
                         group.study = "Predated",
                         group.ref = "Naive")

results.vn.dCt <- RQ_dCt(data = data.dCt.exp,
                         do.tests = TRUE,
                         group.study = "Temp",
                         group.ref = "Naive")

results.zero.dCt <- RQ_dCt(data = data.dCt.exp,
                         do.tests = FALSE,
                         group.study = "Naive",
                         group.ref = "Naive")

print(results.pn.dCt)
print(results.vn.dCt)

naive_comparison <- results.pn.dCt %>%
  filter(Gene %in% c("Lac 1")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Predated vs Naive")

vnaive_comparison <- results.vn.dCt %>%
  filter(Gene %in% c("Lac 1")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "37 °C vs Naive")

zero_comparison <- results.zero.dCt %>%
  filter(Gene %in% c("Lac 1")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Naive vs Naive")

plot_data <- bind_rows(naive_comparison, vnaive_comparison, zero_comparison)

plot_data <- plot_data %>%
  mutate(Comparison = fct_recode(Comparison,
                                 "Predated" = "Predated vs Naive",
                                 "37 °C" = "37 °C vs Naive",
                                 "Naive" = "Naive vs Naive"))
plot_data <- plot_data %>%
  mutate(Comparison = fct_relevel(Comparison,
                                  "Naive",
                                  "Predated",
                                  "37 °C"))

log10FCh.plot <- ggplot(plot_data, aes(x = Comparison, y = FCh, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  ylab("Relative Gene Expression")+
  xlab("")+
  scale_fill_manual(values = c("Lac 1" = "#5e3c99"))+
  theme_classic()+
  theme(legend.text = element_text(face = "italic", size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))

print(log10FCh.plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/relative_gene_expre.png", 
       plot = log10FCh.plot, 
       width = 10, 
       height = 8, 
       dpi = 600)
