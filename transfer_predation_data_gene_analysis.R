library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)
library(RefSeeker)

ct_data_1 <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/ACT_GAP_LAC_transfer_predation_samples_4_nov_2026.txt",
                      header = TRUE, sep = "\t") 

colnames(ct_data_1) <- c("Group", "Well", "Target", "Content", "Sample", "Cq")

ct_data_2 <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/URE_CAP_COX_transfer_predation_ct_values.txt",
                        header = TRUE, sep = "\t")

colnames(ct_data_2) <- c("Well", "Target", "Content", "Sample", "Group", "Cq")


#Transfer Data path into a txt table for dCt analysis ----
selected_ct_1 <- ct_data_1 %>%
  dplyr::select("Group", "Target", "Content", "Sample", "Cq")

selected_ct_2 <- ct_data_2 %>%
  dplyr::select("Group", "Target", "Content", "Sample", "Cq")

combined_delta_data <- rbind(selected_ct_1, selected_ct_2) %>%
  mutate(Cq = as.numeric(Cq)) %>%
  filter(Sample != "NTC", Sample != "NRT", !is.na(Cq))

write.table(combined_delta_data, file = "combined_transfer_predation_data.txt", 
            sep = "\t", row.names = FALSE,
            quote = FALSE)


#Gene Expression Analysis function ----
mean_cq <- function(data){
  data %>%
    group_by(Target, Sample, Group) %>%
    summarise(Cq = mean (Cq, na.rm = TRUE), .groups = "drop")
}

average_hk <- function(data, hk_genes) {
  data %>%
    filter(Target %in% hk_genes) %>%
    group_by(Sample) %>%
    summarise(hk_Ct = mean(Cq, na.rm = TRUE), .groups = "drop")
}

deltaCt <- function(data, hk_genes, target_gene){
  hk = average_hk (data, hk_genes)
  
  data %>%
    filter(Target %in% target_gene) %>%
    left_join(hk, by = "Sample") %>%
    mutate(dCt = Cq - hk_Ct)
}

ddct <- function(data, group.ref){
  group.ref.mean <- data %>%
    filter(Group == group.ref) %>%
    group_by(Target) %>%
    summarise(ref_dCt = mean(dCt, na.rm = TRUE), .groups = "drop")
  
  ddct.result <- data %>%
    #filter(Group == group.study) %>%
    left_join(group.ref.mean, by = "Target") %>%
    mutate(ddCt = dCt - ref_dCt)
}

fCh <- function(data){
  data %>%
  mutate(fCh = 2^(-ddCt))
}

log2fCh <- function(data) {
  data %>% 
    mutate(log2fCh = log2(fCh))
}

AP_dCt <- function(data, hk_genes, target_gene, group.ref){
  data %>%
    mean_cq() %>%
    deltaCt(hk_genes, target_gene)%>%
    ddct(group.ref) %>%
    fCh() %>%
    log2fCh()
}

# Delta Ct Manual ----
result_dCt <- AP_dCt(combined_delta_data, 
                            hk_genes = c("Actin", "GAPDH"), 
                            target_gene = c("Lac 1", "Ure 5", "COX 1", "CAP 10"),
                     group.ref = "Naive")
print(result_dCt)

#plot normalised count ----

plot_normalised_data <- result_dCt %>%
  dplyr::select(Target, Sample, dCt, Group)

group_cols <- c("Naive" = "#66c2a5",
                "Temp" = "#74a9cf", 
                "Predated" = "#a50f15")

normalised_plot <- plot_normalised_data %>% 
  ggplot(aes(x = Sample, y = dCt, fill = Group))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ Target, scales = "free_y")+
  scale_fill_manual(values = group_cols)+
  scale_x_discrete(labels = c("CECN 24h C1" = "Chunk", 
                              "CECN 24h F1" = "Filter paper",
                              "CECN 24h W1" = "Wash",
                              "CECN 48h F1" = "48h Filter paper",
                              "37CN 1" = "37 °C 1",
                              "37CN 2" = "37 °C 2",
                              "25CN 1" = "25 °C 1",
                              "25CN 2" = "25 °C 2",
                              "25Cn 3" = "25 °C 3"))+
  theme_classic() +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=14),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Normalised Expression")

print(normalised_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/normalised_counts_transfer_qpcr.png", 
       plot = normalised_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)

# plot log fold change ----

plot_log_data <- result_dCt %>%
  filter(Group != "Naive")%>%
  dplyr::select(Target, Sample, Group, log2fCh)

gene_cols <- c("Lac 1" = "#5e3c99",
                "Ure 5" = "#e66101",
                "CAP 10" = "#b2abd2",
                "COX 1" = "#fdb863")

log_plot <- plot_log_data %>%
  ggplot(aes(x = Sample, y = log2fCh, fill = Target))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_grid(Target ~ Group, scales = "free", space = "free_x")+
  scale_fill_manual(values = gene_cols)+
  scale_x_discrete(labels = c("CECN 24h C1" = "Chunk", 
                              "CECN 24h F1" = "Filter paper",
                              "CECN 24h W1" = "Wash",
                              "CECN 48h F1" = "48h Filter paper",
                              "37CN 1" = "37 °C 1",
                              "37CN 2" = "37 °C 2"))+
  theme_classic() +
  theme(strip.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 18, face = "italic"),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=14),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Gene Expression")
print(log_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/relative_gene_expression_transfer_qpcr.png", 
       plot = log_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)

# plots with old data ----

pick_predation_data <- read.table("D:/Documents/University/2026/Results/Data analysis/Code/test_predation_gene_analysis_data.txt",
                                  header = TRUE, sep = "\t")

combined_analysis <- rbind(result_dCt, pick_predation_data)

plot_log_data_all <- combined_analysis %>%
  filter(Group != "Naive")%>%
  dplyr::select(Target, Sample, Group, log2fCh)

gene_cols <- c("Lac 1" = "#5e3c99",
               "Ure 5" = "#e66101",
               "CAP 10" = "#b2abd2",
               "COX 1" = "#fdb863")

log_plot <- plot_log_data_all %>%
  ggplot(aes(x = Sample, y = log2fCh, fill = Target))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_grid(Target ~ Group, scales = "free", space = "free_x")+
  scale_fill_manual(values = gene_cols)+
  scale_x_discrete(labels = c("CECN 24h C1" = "Chunk", 
                              "CECN 24h F1" = "Filter paper",
                              "CECN 24h W1" = "Wash",
                              "CECN 48h F1" = "48h Filter paper",
                              "37CN 1" = "37 °C 1",
                              "37CN 2" = "37 °C 2",
                              "CNCE(25)" = "48h Picked", 
                              "CN(25)" = "25 °C",
                              "CN(37)" = "37 °C"))+
  theme_classic() +
  theme(strip.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 18, face = "italic"),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=14),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Relative Gene Expression")
print(log_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/relative_gene_expression_all_predation_tests.png", 
       plot = log_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)



plot_normalised_data_all <- combined_analysis %>%
  dplyr::select(Target, Sample, dCt, Group)

group_cols <- c("Naive" = "#66c2a5",
                "Temp" = "#74a9cf", 
                "Predated" = "#a50f15")

normalised_plot <- plot_normalised_data_all %>% 
  ggplot(aes(x = Sample, y = dCt, fill = Group))+
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ Target, scales = "free_y")+
  scale_fill_manual(values = group_cols)+
  scale_x_discrete(labels = c("CECN 24h C1" = "Chunk", 
                              "CECN 24h F1" = "Filter paper",
                              "CECN 24h W1" = "Wash",
                              "CECN 48h F1" = "48h Filter paper",
                              "37CN 1" = "37 °C 1",
                              "37CN 2" = "37 °C 2",
                              "25CN 1" = "25 °C 1",
                              "25CN 2" = "25 °C 2",
                              "25Cn 3" = "25 °C 3",
                              "CNCE(25)" = "48h Picked", 
                              "CN(25)" = "25 °C",
                              "CN(37)" = "37 °C"))+
  theme_classic() +
  theme(strip.text = element_text(face = "italic"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=14),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(y = "Normalised Expression")

print(normalised_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/normalised_counts_all_predation_tests.png", 
       plot = normalised_plot, 
       width = 10, 
       height = 8, 
       dpi = 600)
# Check stability, collapse into only genes, then run rs ----
# ct_clean <- ct_data_1 %>%
#   filter(Sample != "NTC", Sample != "NRT", Target != "Lac 1") %>%
#   mutate(Cq = as.numeric(Cq)) %>%
#   group_by(Sample, Target) %>%
#   summarise(Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
#   mutate(Cq = ifelse(is.nan(Cq), NA, Cq))
# 
# ct_samples_rows <- ct_clean %>%
#   pivot_wider(names_from = Target, values_from = Cq) %>%
#   dplyr::select(-Sample)
# 
# print(sapply(ct_samples_rows, sd))
# 
# print(cor(ct_samples_rows$Actin, ct_samples_rows$GAPDH))
# print(sd(ct_samples_rows$Actin - ct_samples_rows$GAPDH))
# 
# ct_clean_plot <- ct_clean %>%
#   filter(Target %in% c("Actin", "GAPDH")) %>%
#   ggplot(aes(x = Sample, y = Cq, colour = Target)) +
#   geom_point() +
#   theme_classic()
# 
# print(ct_clean_plot)

#LogFCh plot ----
# naive_comparison_log <- results.pn.dCt %>%
#   filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1"))%>%
#   dplyr::select(Gene, log10FCh) %>%
#   mutate(Comparison = "Predated vs Naive")
# nrow(naive_comparison)
# vnaive_comparison_log <- results.vn.dCt %>%
#   filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
#   dplyr::select(Gene, log10FCh) %>%
#   mutate(Comparison = "37 °C vs Naive")
# 
# zero_comparison_log <- results.zero.dCt %>%
#   filter(Gene %in% c("Lac 1", "Ure 5", "CAP 10", "COX 1")) %>%
#   dplyr::select(Gene, log10FCh) %>%
#   mutate(Comparison = "Naive vs Naive")
# 
# plot_data_log <- bind_rows(naive_comparison_log, vnaive_comparison_log, zero_comparison_log)
# 
# plot_log <- plot_data_log %>%
#   mutate(Comparison = fct_recode(Comparison,
#                                  "Predated" = "Predated vs Naive",
#                                  "37 °C" = "37 °C vs Naive",
#                                  "Naive" = "Naive vs Naive"),
#          Comparison = fct_relevel(Comparison,
#                                   "Naive",
#                                   "Predated",
#                                   "37 °C"))
# 
# log10FCh.plot <- ggplot(plot_log, aes(x = Comparison, y = log10FCh, fill = Gene)) +
#   geom_col(position = position_dodge(width = 0.7), width = 0.6) +
#   ylab("Relative Gene Expression")+
#   xlab("")+
#   scale_fill_manual(values = c("Lac 1" = "#5e3c99",
#                                "Ure 5" = "#e66101",
#                                "CAP 10" = "#b2abd2",
#                                "COX 1" = "#fdb863"),
#                     labels = function(x) parse(text = paste0("italic('", x, "')")))+
#   theme_classic()+
#   theme(legend.text = element_text(face = "italic", size = 20),
#         legend.title = element_text(size = 20),
#         axis.text = element_text(size=20),
#         axis.title = element_text(size=24))
# 
# print(log10FCh.plot)
# ggsave("C:/Users/User/OneDrive/Documents/Research/Results/transfer_relative_log10.png", 
#        plot = log10FCh.plot, 
#        width = 10, 
#        height = 8, 
#        dpi = 600)


#ddCt Analysis ----
# ddct.pn.dct <- RQ_ddCt(data = data.dCt.exp, 
#                        do.tests = TRUE,
#                        group.study = "Predated",
#                        group.ref = "Naive")
# 
# ddct.vn <- RQ_ddCt(data = data.dCt.exp, 
#                    do.tests = FALSE, 
#                    group.study = "Temp",
#                    group.ref = "Naive")
# 
# print(ddct.pn.dct)
# 
# ddct_log2_FC_p <- data.frame (
#   Gene = ddct.pn.dct$Gene,
#   log2FCh = log2(ddct.pn.dct$FCh))
# 
# ddct_log2_p <- ddct_log2_FC_p %>% rename(Predated = log2FCh)
# 
# ddct_log2_FC_t <- data.frame(
#   Gene = ddct.vn$Gene, 
#   log2FCh = log2(ddct.vn$FCh))
# 
# ddct_log2_t <- ddct_log2_FC_t %>% rename(Temp = log2FCh)
# 
# ddct_log2 <- ddct_log2_p %>% mutate(Temp = ddct_log2_t$Temp)
# 
# print(ddct_log2)
# 
# plot_ddct_log2 <- pivot_longer(
#   ddct_log2, 
#   cols = c(Predated, Temp), 
#   names_to = "Group", 
#   values_to = "log2FCh"
# )
# 
# 
# ddct_plot <- ggplot(plot_ddct_log2, aes(x = Group, y = log2FCh, fill = Gene)) +
#   geom_col(position = position_dodge(width = 0.7), width = 0.6) +
#   geom_hline(yintercept = 0) +
#   facet_wrap(~ Gene, scales = "free_y")+
#   #scale_y_continuous(transform = scales::pseudo_log_trans(base = 10))+
#   scale_fill_manual(values = c("Lac 1" = "#5e3c99",
#                                "Ure 5" = "#e66101",
#                                "CAP 10" = "#b2abd2",
#                                "COX 1" = "#fdb863"),
#                     labels = function(x) parse(text = paste0("italic('", x, "')")))+
#   theme_classic() +
#   theme(strip.text = element_text(face = "italic"),
#         legend.text = element_text(face = "italic", size = 20),
#         legend.title = element_text(size = 20),
#         axis.text = element_text(size=20),
#         axis.title = element_text(size=24))+
#   labs(y = "Relative Gene Expression")
# 
# print(ddct_plot)
# ggsave("C:/Users/User/OneDrive/Documents/Research/Results/transfer_relative_log2ddct.png", 
#        plot = ddct_plot, 
#        width = 10, 
#        height = 8, 
#        dpi = 600)
# 
# print(result_ddct_log2_FC)
