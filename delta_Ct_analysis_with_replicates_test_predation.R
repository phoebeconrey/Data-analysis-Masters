library(RQdeltaCT)
library(tidyverse)
library(coin)
library(ggsignif)

data <- "C:/Users/User/OneDrive/Documents/Research/Data analysis/test_predation_ct_values_r_replicates.txt"

data.Ct <- read_Ct_long(path = data, 
                        sep = "\t",
                        dec = ".",
                        skip = 0,
                        add.column.Flag = FALSE,
                        column.Sample = 2,
                        column.Gene = 3,
                        column.Ct = 4,
                        column.Group = 1)
str(data.Ct)

data.CtF.ready <- make_Ct_ready(data = data.Ct,
                                imput.by.mean.within.groups = TRUE)
data.dCt.exp <- delta_Ct(data = data.CtF.ready,
                         normalise = TRUE,
                         ref = c("GAPDH 3", "Act 1"),
                         transform = TRUE)

results.dCt <- RQ_dCt(data = data.dCt.exp,
                      do.tests = FALSE, 
                      group.study = "Predated",
                      group.ref = "Naive")
head(as.data.frame(arrange(results.dCt, MW_test_p)))
FCh.plot <- FCh_plot(data = results.dCt, 
                     use.p = FALSE, 
                     mode = "depends",
                     FCh.threshold = 2, 
                     angle = 20)
