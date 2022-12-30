#
# Auxiliary functions for plotting results of MEA analysis
# 01.12.2022; AB
#

set.seed(42)

library(tidyverse)
library(purrr)
library(janitor)
library(patchwork)


## loc pal 3 way
pal_loc.3way <- c(Projections = "#82b21e", Cytoplasm = "#006beb", 
                  NS = "gray50", Nucleus = "#236477")
## loc pal 2 way
pal_loc.2way <- c(Neurite = "#82b21e", NS = "gray50", Soma = "#006beb")

## for dnCaf1 and GFP control
pal_mut <- c(GFP = "steelblue3", NS = "gray50", dnCaf1 = "orangered3")
## for variour m6A perturbations
pal_mut.msa <- c("#f18ea3", "grey50", "#481f31")
pal_msa <- c("#009cb8", "grey50", "#d64ecf")

##
ggplot2::theme_set(theme_bw(base_size = 14) +
                     theme(panel.border = element_blank(),
                           axis.line = element_line(size = 0.25, color = "black"), 
                           axis.ticks = element_line(size = 0.25, color = "black"),
                           axis.text = element_text(color = "black"),
                           legend.background = element_blank()
                     )
)


## Auxiliary function to perform a permutation test
perm_test <- function(tab, col_to_perm, col_to_test, n_perm = 5000) {
  grps <- unique(tab[[col_to_perm]])
  purrr::map_dfr(1:n_perm,
                 function(iter) {
                   tab[[col_to_perm]] <- sample(grps, 
                                                dim(tab)[1],
                                                replace = T)
                   
                   tibble(n_iter = iter,
                          p.value = wilcox.test(formula = !!sym(col_to_test) ~ !!sym(col_to_perm), 
                                                data = tab)$p.value)
                 })
}