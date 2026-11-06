library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)
library(RefSeeker)

ct_data_1 <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/ACT_GAP_LAC_transfer_predation_samples_4_nov_2026.txt",
                      header = TRUE, sep = "\t") 

colnames(ct_data_1) <- c("Group", "Well", "Target", "Content", "Sample", "Cq")

ct_data_2 <- read.table("C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/URE_CAP_COX_transfer_predation_ct_values.txt",
                        header = TRUE, sep = "\t")

colnames(ct_data_2) <- c("Well", "Target", "Content", "Sample", "Group", "Cq")

# Check stability, collapse into only genes, then run rs ----
ct_clean <- ct_data_1 %>%
  filter(Sample != "NTC", Sample != "NRT", Target != "Lac 1") %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Sample, Target) %>%
  summarise(Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  mutate(Cq = ifelse(is.nan(Cq), NA, Cq))

ct_samples_rows <- ct_clean %>%
  pivot_wider(names_from = Target, values_from = Cq) %>%
  dplyr::select(-Sample)

print(sapply(ct_samples_rows, sd))

print(cor(ct_samples_rows$Actin, ct_samples_rows$GAPDH))
print(sd(ct_samples_rows$Actin - ct_samples_rows$GAPDH))

ct_clean_plot <- ct_clean %>%
  filter(Target %in% c("Actin", "GAPDH")) %>%
  ggplot(aes(x = Sample, y = Cq, colour = Target)) +
  geom_point() +
  theme_classic()

print(ct_clean_plot)

#Transfer Data path into a txt table for dCt analysis ----
selected_ct_1 <- ct_data_1 %>%
  dplyr::select("Group", "Target", "Content", "Sample", "Cq")

selected_ct_2 <- ct_data_2 %>%
  dplyr::select("Group", "Target", "Content", "Sample", "Cq")

combined_delta_data <- rbind(selected_ct_1, selected_ct_2) %>%
  filter(Sample != "NTC", Sample != "NRT", Cq != "NaN")

write.table(combined_delta_data, file = "combined_transfer_predation_data.txt", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)

#Delta Ct Analysis ----
delta_data <- read_Ct_long(path = "C:/Users/pheebsy/OneDrive/Documents/Research/Data analysis/combined_transfer_predation_data.txt",
                           sep = "\t",
                           dec = ".",
                           skip = 0,
                           add.column.Flag = FALSE,
                           column.Sample = 4,
                           column.Gene = 2,
                           column.Ct = 5,
                           column.Group = 1)

data.CtF.ready <- make_Ct_ready(data = delta_data,
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
                         do.tests = FALSE,
                         group.study = "Temp",
                         group.ref = "Naive")

results.zero.dCt <- RQ_dCt(data = data.dCt.exp,
                         do.tests = FALSE,
                         group.study = "Naive",
                         group.ref = "Naive")

print(results.pn.dCt)
print(results.vn.dCt)

#FCh comparison ----

naive_comparison <- results.pn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1"))%>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Predated vs Naive")
nrow(naive_comparison)
vnaive_comparison <- results.vn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "37 °C vs Naive")

zero_comparison <- results.zero.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
  dplyr::select(Gene, FCh) %>%
  mutate(Comparison = "Naive vs Naive")

plot_data <- bind_rows(naive_comparison, vnaive_comparison, zero_comparison)

plot_data <- plot_data %>%
  mutate(Comparison = fct_recode(Comparison,
                                 "Predated" = "Predated vs Naive",
                                 "37 °C" = "37 °C vs Naive",
                                 "Naive" = "Naive vs Naive"),
         Comparison = fct_relevel(Comparison,
                                  "Naive",
                                  "Predated",
                                  "37 °C"))

FCh.plot <- ggplot(plot_data, aes(x = Comparison, y = FCh, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  ylab("Relative Gene Expression")+
  xlab("")+
  scale_fill_manual(values = c("Lac 1" = "#5e3c99",
                               "Ure 5" = "#e66101",
                               "CAP 10" = "#b2abd2",
                               "COX 1" = "#fdb863"))+
  theme_classic()+
  theme(legend.text = element_text(face = "italic", size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))

print(FCh.plot)
ggsave("C:/Users/Pheebsy/OneDrive/Documents/Research/Results/transfer_relative_FCh.png", 
       plot = FCh.plot, 
       width = 10, 
       height = 8, 
       dpi = 600)

#LogFCh plot ----
naive_comparison_log <- results.pn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1"))%>%
  dplyr::select(Gene, log10FCh) %>%
  mutate(Comparison = "Predated vs Naive")
nrow(naive_comparison)
vnaive_comparison_log <- results.vn.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
  dplyr::select(Gene, log10FCh) %>%
  mutate(Comparison = "37 °C vs Naive")

zero_comparison_log <- results.zero.dCt %>%
  filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
  dplyr::select(Gene, log10FCh) %>%
  mutate(Comparison = "Naive vs Naive")

plot_data_log <- bind_rows(naive_comparison_log, vnaive_comparison_log, zero_comparison_log)

plot_log <- plot_data_log %>%
  mutate(Comparison = fct_recode(Comparison,
                                 "Predated" = "Predated vs Naive",
                                 "37 °C" = "37 °C vs Naive",
                                 "Naive" = "Naive vs Naive"),
         Comparison = fct_relevel(Comparison,
                                  "Naive",
                                  "Predated",
                                  "37 °C"))

log10FCh.plot <- ggplot(plot_log, aes(x = Comparison, y = log10FCh, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  ylab("Relative Gene Expression")+
  xlab("")+
  scale_fill_manual(values = c("Lac 1" = "#5e3c99",
                               "Ure 5" = "#e66101",
                               "CAP 10" = "#b2abd2",
                               "COX 1" = "#fdb863"),
                    labels = function(x) parse(text = paste0("italic('", x, "')")))+
  theme_classic()+
  theme(legend.text = element_text(face = "italic", size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))

print(log10FCh.plot)
ggsave("C:/Users/Pheebsy/OneDrive/Documents/Research/Results/transfer_relative_log10.png", 
       plot = log10FCh.plot, 
       width = 10, 
       height = 8, 
       dpi = 600)


#ddCt Analysis ----
ddct.pn.dct <- RQ_ddCt(data = data.dCt.exp, 
                       do.tests = TRUE,
                       group.study = "Predated",
                       group.ref = "Naive")

ddct.vn <- RQ_ddCt(data = data.dCt.exp, 
                   do.tests = FALSE, 
                   group.study = "Temp",
                   group.ref = "Naive")

print(ddct.pn.dct)

ddct_log2_FC_p <- data.frame (
  Gene = ddct.pn.dct$Gene,
  log2FCh = log2(ddct.pn.dct$FCh))

ddct_log2_p <- ddct_log2_FC_p %>% rename(Predated = log2FCh)

ddct_log2_FC_t <- data.frame(
  Gene = ddct.vn$Gene, 
  log2FCh = log2(ddct.vn$FCh))

ddct_log2_t <- ddct_log2_FC_t %>% rename(Temp = log2FCh)

ddct_log2 <- ddct_log2_p %>% mutate(Temp = ddct_log2_t$Temp)

print(ddct_log2)

plot_ddct_log2 <- pivot_longer(
  ddct_log2, 
  cols = c(Predated, Temp), 
  names_to = "Group", 
  values_to = "log2FCh"
)


ddct_plot <- ggplot(plot_ddct_log2, aes(x = Group, y = log2FCh, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ Gene, scales = "free_y")+
  #scale_y_continuous(transform = scales::pseudo_log_trans(base = 10))+
  scale_fill_manual(values = c("Lac 1" = "#5e3c99",
                               "Ure 5" = "#e66101",
                               "CAP 10" = "#b2abd2",
                               "COX 1" = "#fdb863"),
                    labels = function(x) parse(text = paste0("italic('", x, "')")))+
  theme_classic() +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(face = "italic", size = 20),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))+
  labs(y = "Relative Gene Expression")

print(ddct_plot)

print(result_ddct_log2_FC)
