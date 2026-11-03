library(tidyr)
library(dplyr)
library(drc)
library(ggplot2)

survival_data <- read.table("C:/Users/User/OneDrive/Documents/Research/Data analysis/first_predation_LT50_data.txt",
                            header = TRUE, sep = "\t")

colnames(survival_data) <- c("Time", "Replicate", "E. coli", 
                             "C. neoformans", "S. cerevisiae", 
                             "Starved")

data_long <- survival_data %>% pivot_longer(cols = -c(Time, Replicate), 
                                            names_to = "Treatment", values_to = "Alive")
starting_counts <- data_long %>%
  filter(Time == 0) %>%
  dplyr::select(Treatment, Replicate, Total = Alive)
data_fixed <- data_long %>%
  left_join(starting_counts, by = c("Treatment", "Replicate"))

adjusted_data <- data_fixed %>% 
  mutate(Dead = Total - Alive,
         Survival_pct = (Alive/Total)*100)
mean_data <- adjusted_data %>%
  group_by(Treatment, Time) %>%
  summarise(Survival_pct = mean(Survival_pct), .groups = "drop")

survival_plot <- ggplot()+
 # geom_step(data = adjusted_data, 
     #       aes(x = Time, y = Survival_pct, color = Treatment), alpha = 0.3, direction = "hv") +
  geom_step(data = mean_data,
            aes(x = Time, y = Survival_pct, color = Treatment), direction = "hv", linewidth = 1.5)+
  geom_point(data = mean_data, 
              aes (x = Time, y = Survival_pct, color = Treatment), size = 3) +
  geom_point(data = adjusted_data, 
             aes (x = Time, y = Survival_pct, color = Treatment), alpha = 0.3, size = 3) +
  labs(title = NULL, x = "Time (hours)", 
       y = "% Survival", colour = "Predation Group") +
  scale_color_manual(values = c("E. coli" = "#e66101",
                                "C. neoformans" = "#5e3c99",
                                "S. cerevisiae" = "#b2abd2",
                                "Starved" = "#fdb863"),
                     labels = c(expression(italic("C. neoformans")),
                                expression(italic("E. coli")),
                                expression(italic("S. cerevisiae")),
                                libr"Starved"))+
  theme_classic()+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24))
print(survival_plot)
ggsave("C:/Users/User/OneDrive/Documents/Research/Results/LT50 Plot.png", plot = survival_plot, width = 10, height = 8, dpi = 600)

lt50_model <- drm(Dead/Total ~ Time, weights = Total, curveid = Treatment,
                  data = adjusted_data, fct = LL.2(), type = "binomial")
ED(lt50_model, 50, interval = "delta")
compParm(lt50_model, "e")
print(ED)
print(compParm)