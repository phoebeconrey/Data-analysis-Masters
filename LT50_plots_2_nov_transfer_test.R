library(tidyr)
library(dplyr)
library(drc)
library(ggplot2)
library(ggbreak)

survival_data <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/24h_arrest_predation_transfer_type.txt",
                            header = TRUE, sep = "\t")

colnames(survival_data) <- c("Time", "Transfer", "E. coli", "C. neoformans", "S. cerevisiae")

data_long <- survival_data %>% pivot_longer(cols = -c(Time, Transfer),
                                            names_to = "Food", values_to = "Alive") %>%
  mutate(Food = factor (Food, levels = c("E. coli", "C. neoformans", "S. cerevisiae")))


starting_counts <- data_long %>%
  filter(Time == 0) %>%
  dplyr::select(Food, Transfer, Total = Alive)
data_fixed <- data_long %>%
  left_join(starting_counts, by = c("Food", "Transfer"))

adjusted_data <- data_fixed %>%
  mutate(Dead = Total - Alive,
         Survival_pct = (Alive/Total)*100)
mean_data <- adjusted_data %>% 
  group_by(Food, Time) %>%
  summarise(Survival_pct = mean(Survival_pct), .groups = "drop")

food_cols <- c("E. coli" = "#e66101",
               "C. neoformans" = "#5e3c99",
               "S. cerevisiae" = "#b2abd2")
sample_cols <- c("Filter" = "#FDB833",
                 "Chunk" = "#05204A",
                 "Wash" = "#DD1C1A")
#Normal Death curve ----
survival_plot <- ggplot() +
  geom_step(data = mean_data,
            aes(x = Time, y = Survival_pct, colour = Food), direction = "hv", linewidth = 1.5)+
  geom_point(data = adjusted_data, 
             aes (x = Time, y = Survival_pct, fill = Transfer), alpha = 0.3, size = 3) +
  geom_point(data = mean_data, 
             aes (x = Time, y = Survival_pct, colour = Food), size = 3)+
  labs(title = NULL, x = "Time (hours)",
       y = "% Survival", colour = "Food Source")+
  
  scale_colour_manual(name = "Food Source",values = food_cols,
                      labels = function(x) parse(text = paste0("italic('", x, "')")))+
  scale_fill_manual(name = "Transfer Type", values = sample_cols)+
  theme_classic()+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))
print(survival_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50_transfer_test_plot.png",
       plot = survival_plot, width = 10, height = 8, dpi = 600)


#Transfer Death Curve ----
trans_cols <- c("Filter" = "#e66101",
               "Chunk" = "#5e3c99",
               "Wash" = "#b2abd2")

transfer_data <- adjusted_data %>% filter(Time <=3)

mean_transfer_data <- transfer_data %>%
  filter(Time == 3) %>%
  group_by(Transfer) %>%
  summarise(mean_surv = mean(Survival_pct, na.rm = TRUE), .groups = "drop")



other_transfer_plot <- ggplot()+
  geom_jitter(data = transfer_data, 
              aes(x = Time, y = Survival_pct, colour = Food),
              width = 0.1, size = 2)+
  geom_hline(data = mean_transfer_data, 
                aes(yintercept = mean_surv), colour = "black", linewidth = 0.8)+
  facet_wrap(~ Transfer) +
  labs(x = "Time(hours", y = "% Survival")+
  scale_colour_manual(values = food_cols)+
  theme_classic()+
  theme(legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))
print(other_transfer_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50_transfer_survival_plot.png", 
       plot = other_transfer_plot, width = 10, height = 8, dpi = 600)

#transfer_plot <- adjusted_data %>%
  #filter(Time <=3) %>%
  #ggplot(aes(x = Time, y = Survival_pct, colour = Transfer)) +
  #geom_point(alpha = 0.4, size = 2)+
  #facet_wrap(~ Food)+
  #labs(title = NULL, x = "Time (hours)",
      # y = "% Survival", colour = "Transfer type")+
  #scale_colour_manual(values = trans_cols)+
  #theme_classic()
#print(transfer_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50_transfer_comparison_test_plot.png", plot = other_transfer_plot, width = 10, height = 8, dpi = 600)

lt50_model <- drm(Dead/Total ~ Time, weights = Total, curveid = Food,
                  data = adjusted_data, fct = LL.2(), type = "binomial")
ED(lt50_model, 50, interval = "delta")
compParm(lt50_model, "e")
print(lt50_model)
