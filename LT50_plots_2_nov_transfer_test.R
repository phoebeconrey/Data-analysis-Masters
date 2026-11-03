library(tidyr)
library(dplyr)
library(drc)
library(ggplot2)

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

mean_data_transfer <- adjusted_data %>% 
  group_by(Food, Transfer, Time) %>%
  filter(0 <= Time => 3) %>%
  summarise(Survival_pct = mean(Survival_pct), .groups = "drop")

food_cols <- c("E. coli" = "#e66101",
               "C. neoformans" = "#5e3c99",
               "S. cerevisiae" = "#b2abd2")

survival_plot <- ggplot() +
  geom_step(data = mean_data,
            aes(x = Time, y = Survival_pct, colour = Food), direction = "hv", linewidth = 1.5)+
  geom_point(data = adjusted_data, 
             aes (x = Time, y = Survival_pct, colour = Food), alpha = 0.3, size = 3) +
  geom_point(data = mean_data, 
             aes (x = Time, y = Survival_pct, colour = Food), size = 3)+
  labs(title = NULL, x = "Time (hours)",
       y = "% Survival", colour = "Food Source")+
  
  scale_colour_manual(values = food_cols,
                      labels = function(x) parse(text = paste0("italic('", x, "')")))+
  theme_classic()+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))
print(survival_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50_transfer_test_plot.png", plot = survival_plot, width = 10, height = 8, dpi = 600)

trans_cols <- c("Filter" = "#e66101",
               "Chunk" = "#5e3c99",
               "Wash" = "#b2abd2")

transfer_plot <- ggplot() +
  geom_step(data = mean_data_transfer,
            aes(x = Time, y = Survival_pct, colour = Transfer), direction = "hv", linewidth = 1.5)+
  geom_point(data = adjusted_data, 
             aes (x = Time, y = Survival_pct, colour = Transfer), alpha = 0.3, size = 3) +
  geom_point(data = mean_data_transfer, 
             aes (x = Time, y = Survival_pct, colour = Transfer), size = 3)+
  facet_wrap(~ Food)+
  labs(title = NULL, x = "Time (hours)",
       y = "% Survival", colour = "Transfer type")+
  
  scale_colour_manual(values = trans_cols,
                      labels = function(x) parse(text = paste0("italic('", x, "')")))+
  theme_classic()+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))
print(transfer_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50_transfer_comparison_test_plot.png", plot = transfer_plot, width = 10, height = 8, dpi = 600)

lt50_model <- drm(Dead/Total ~ Time, weights = Total, curveid = Food,
                  data = adjusted_data, fct = LL.2(), type = "binomial")
ED(lt50_model, 50, interval = "delta")
compParm(lt50_model, "e")
print(lt50_model)
