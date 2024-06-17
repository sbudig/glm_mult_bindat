# Setup -------------------------------------------------------------------

# Functions
source("bin_figures_source.R")

# required packages

library(tidyverse)
library(patchwork)
library(reshape2)

# Regression model --------------------------------------------------------

## Figure 1 - FWERW -------------------------------------------------------

# read and append all intermediate results from FWER of regression models
bin_reg_FWER <- r_a_a_bin_reg_FWER_csvs("..\\Simulations\\Sim_Regression\\intermediate_results")

# summarise and calculate percentages for FWER of regression models intermediate results
bin_reg_FWER <- s_a_cp_bin_reg_FWER(bin_reg_FWER)

# column percenti indicates the percentage of simulations,
# where at least one rejection has taken place

bin_reg_FWER$model <- as.factor(bin_reg_FWER$model)
bin_reg_FWER$model <- relevel(bin_reg_FWER$model, ref = "glm")

J_FWERW_label <- c(expression(paste(J[0], ": ", "3")),
                   expression(paste(J[0], ": ", "5")),
                   expression(paste(J[0], ": ", "10")))

N_FWERW_label <- c(expression(paste("N", ": ", "100")),
                   expression(paste("N", ": ", "200")),
                   expression(paste("N", ": ", "400")))

Cor_FWERW_label <- c(expression(paste("cor", ": ", "0")),
                     expression(paste("cor", ": ", "0.5")),
                     expression(paste("cor", ": ", "0.9")))

bin_reg_FWER$J <- factor(bin_reg_FWER$J, labels = J_FWERW_label)
bin_reg_FWER$N <- factor(bin_reg_FWER$N, labels = N_FWERW_label)
bin_reg_FWER$cor <-
  factor(bin_reg_FWER$cor, labels = Cor_FWERW_label)

plot_reg_FWER <- ggplot(bin_reg_FWER, aes(x = pis, y = percenti)) +
  geom_point(aes(colour = method, shape = model), size = 3) +
  facet_grid(cor ~ N + J, labeller = "label_parsed") +
  ylab(expression("FWERW" ~ (alpha))) +
  xlab(expression("Probability of occurence" ~ (pi))) +
  geom_hline(yintercept = 0.05) +
  geom_hline(yintercept = 0.04580994, linetype = "dashed") +
  geom_hline(yintercept = 0.05445453, linetype = "dashed") +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  scale_color_discrete(name = "Method",
                       labels = c("Bonferroni", "Bootstrap", "MMM")) +
  theme_bw(base_size = 18) + theme(legend.position = "bottom") +
  scale_x_continuous(labels = prettyZero, guide = guide_axis(angle = 90)) +
  scale_y_continuous(labels = prettyZero, breaks = c(0, 0.02, 0.04, 0.05, 0.06))

#plot_reg_FWER

ggsave(
  ".\\Budig_Figure_1.eps",
  plot = plot_reg_FWER ,
  width = 14,
  height = 10,
  device = "eps",
  dpi = 800
)

## Figure 2 - Power -------------------------------------------------------

# read and append all intermediate results from FWER of regression models
bin_reg_POWER <- r_a_a_bin_reg_Power_csvs("..\\Simulations\\Sim_Regression\\intermediate_results")

# summarise and calculate percentages for FWER of regression models intermediate results
bin_reg_POWER <- s_a_cp_bin_reg_Power(bin_reg_POWER)

# columns
# rej_p     # How many rejected per Dataset in percent
# reji_p    # At least one rejected per Dataset in percent
# rejcor_p  # how often all h0 are true in percent
#
# pow_p     # How many differences found per Dataset in percent
# powi_p   # At least one difference found per Dataset in percent
# powcor_p  # How often find all differences in percent

bin_reg_POWER$model <-
  factor(bin_reg_POWER$model, levels = c("glm", "bglm"))

J_label <- c(expression(paste(J[0], ": ", "2", ", ", J[A], ": ", "8")),
             expression(paste(J[0], ": ", "8", ", ", J[A], ": ", "2")))

N_label <- c(expression(paste("N", ": ", "100")),
             expression(paste("N", ": ", "200")),
             expression(paste("N", ": ", "400")))

Cor_label <- c(expression(paste("cor", ": ", "0")),
               expression(paste("cor", ": ", "0.5")),
               expression(paste("cor", ": ", "0.9")))

bin_reg_POWER$J0_JA <- factor(bin_reg_POWER$J0_JA, labels = J_label)
bin_reg_POWER$N <- factor(bin_reg_POWER$N, labels = N_label)
bin_reg_POWER$cor <- factor(bin_reg_POWER$cor, labels = Cor_label)


plot_reg_POWER <-
  ggplot(bin_reg_POWER, aes(x = slope, y = powi_p)) +
  geom_line(aes(group = interaction(method, model), colour = method)) +
  geom_point(aes(colour = method, shape = model), size = 3) +
  facet_grid(cor ~ N + J0_JA, labeller = "label_parsed")  +
  ylab("Power") +
  xlab(expression("Slope" ~ (beta["j1"]))) +
  #scale_color_manual(values = c("#117733","#661100","#0072B2"),guide_legend(title="N. of Trials"))+
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = prettyZero, guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = prettyZero) +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  scale_color_discrete(name = "Method",
                       labels = c("Bonferroni", "Bootstrap", "MMM"))

#plot_reg_POWER

ggsave(
  ".\\Budig_Figure_2.eps",
  plot = plot_reg_POWER ,
  width = 14,
  height = 10,
  device = "eps",
  dpi = 800
)

## Figure 3 - FWERS -------------------------------------------------------

plot_reg_FWERS <-
  ggplot(bin_reg_POWER, aes(x = slope, y = reji_p)) +
  geom_point(aes(colour = method, shape = model), size = 3) +
  facet_grid(cor ~ N + J0_JA, labeller = "label_parsed") +
  # geom_line(aes(group = interaction(model,method), colour = method)) +
  theme_bw(base_size = 18) +
  ylab(expression("FWERS" ~ (alpha))) +
  xlab(expression("Slope" ~ (beta["1j"]))) +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = prettyZero, guide = guide_axis(angle = 45)) +
  scale_y_continuous(labels = prettyZero) +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  scale_color_discrete(name = "Method",
                       labels = c("Bonferroni", "Bootstrap", "MMM"))

#plot_reg_FWERS

ggsave(
  ".\\Budig_Figure_3.eps",
  plot = plot_reg_FWERS ,
  width = 14,
  height = 10,
  device = "eps",
  dpi = 800
)

# Group Comparison Model --------------------------------------------------

# read and append all intermediate results from group comparison models
bin_gc_ss <- r_a_a_bin_gc_ss_csvs("..\\Simulations\\Sim_Group_comparison\\intermediate_results")

all_gc_ss <- s_a_cp_bin_gc_ss(bin_gc_ss)

all_gc_ss_wide <-  l_t_w_bin_gc_ss(all_gc_ss)

all_gc_ss_wide$model <- factor(all_gc_ss_wide$model)
all_gc_ss_wide$model <- relevel(all_gc_ss_wide$model, ref = "glm")
all_gc_ss_wide$Ncor <-
  factor(all_gc_ss_wide$N):factor(all_gc_ss_wide$cor)
all_gc_ss_wide$Ncor <-
  factor(
    all_gc_ss_wide$Ncor,
    levels = c("100:0", "200:0", "100:0.5", "200:0.5", "100:0.9", "200:0.9")
  )

## Figure 4 - FWER --------------------------------------------------------

fwer_mmm_bt <-
  ggplot(all_gc_ss_wide, aes(x = reji_p_mmm, y = reji_p_bt)) +
  geom_point(aes(shape = model, color = Ncor), size = 3) +
  #scale_color_manual(values = c("#117733","#661100","#0072B2"),guide_legend(title="N. of Trials"))+
  theme_bw(base_size = 18) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "solid") +
  geom_hline(yintercept = 0.05445453, linetype = "dashed") +
  geom_hline(yintercept = 0.04580994 , linetype = "dashed") +
  geom_vline(xintercept = 0.05, linetype = "solid") +
  geom_vline(xintercept = 0.05445453, linetype = "dashed") +
  geom_vline(xintercept = 0.04580994 , linetype = "dashed") +
  #scale_y_continuous(breaks = c(0,0.02,0.04,0.05,0.06))+
  #scale_x_continuous(breaks = c(0,0.02,0.04,0.05,0.06))+
  xlim(0, 0.06) + ylim(0, 0.06) +
  scale_color_manual(values = c(
    "#f6e8c3",
    "#c7eae5",
    "#d8b365",
    "#5ab4ac",
    "#8c510a",
    "#01665e"
  )) +
  xlab("FWER MMM") + ylab("FWER Bootstrap") +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  guides(color = guide_legend(title = "N:cor"), shape = guide_legend(order = 1)) + theme(legend.position =
                                                                                           "bottom")

fwer_mmm_bonf <-
  ggplot(all_gc_ss_wide, aes(x = reji_p_mmm, y = reji_p_bonf)) +
  geom_point(aes(shape = model, color = Ncor), size = 3) +
  #scale_color_manual(values = c("#117733","#661100","#0072B2"),guide_legend(title="N. of Trials"))+
  theme_bw(base_size = 18) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 0.05, linetype = "solid") +
  geom_hline(yintercept = 0.05445453, linetype = "dashed") +
  geom_hline(yintercept = 0.04580994 , linetype = "dashed") +
  geom_vline(xintercept = 0.05, linetype = "solid") +
  geom_vline(xintercept = 0.05445453, linetype = "dashed") +
  geom_vline(xintercept = 0.04580994 , linetype = "dashed") +
  xlim(0, 0.06) + ylim(0, 0.06) +
  scale_color_manual(values = c(
    "#f6e8c3",
    "#c7eae5",
    "#d8b365",
    "#5ab4ac",
    "#8c510a",
    "#01665e"
  )) +
  xlab("FWER MMM") + ylab("FWER Bonferroni") +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  guides(color = guide_legend(title = "N:cor"), shape = guide_legend(order = 1)) + theme(legend.position =
                                                                                           "bottom")

plot_gc_FWER <-
  fwer_mmm_bt + fwer_mmm_bonf + plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

#plot_gc_FWER

ggsave(
  ".\\Budig_Figure_4.eps",
  plot = plot_gc_FWER,
  width = 14,
  height = 7.8,
  device = "eps",
  dpi = 800
)

## Figure 5 - Power -------------------------------------------------------

power_mmm_bt <-
  ggplot(all_gc_ss_wide, aes(x = powi_p_mmm, y = powi_p_bt)) +
  geom_point(aes(shape = model, color = Ncor), size = 3) +
  #scale_color_manual(values = c("#117733","#661100","#0072B2"),guide_legend(title="N. of Trials"))+
  theme_bw(base_size = 18) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c(
    "#f6e8c3",
    "#c7eae5",
    "#d8b365",
    "#5ab4ac",
    "#8c510a",
    "#01665e"
  )) +
  xlab("Power MMM") + ylab("Power Bootstrap") +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  guides(color = guide_legend(title = "N:cor"), shape = guide_legend(order = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme(legend.position = "bottom")

power_mmm_bonf <-
  ggplot(all_gc_ss_wide, aes(x = powi_p_mmm, y = powi_p_bonf)) +
  geom_point(aes(shape = model, color = Ncor), size = 3) +
  #scale_color_manual(values = c("#117733","#661100","#0072B2"),guide_legend(title="N. of Trials"))+
  theme_bw(base_size = 18) +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = c(
    "#f6e8c3",
    "#c7eae5",
    "#d8b365",
    "#5ab4ac",
    "#8c510a",
    "#01665e"
  )) +
  xlab("Power MMM") + ylab("Power Bonferroni") +
  scale_shape_discrete(name = "Model class", labels = c("GLM", "BGLM")) +
  guides(color = guide_legend(title = "N:cor"), shape = guide_legend(order = 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme(legend.position = "bottom")

plot_gc_POWER <-
  power_mmm_bt + power_mmm_bonf + plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

#plot_gc_POWER

ggsave(
  ".\\Budig_Figure_5.eps",
  plot = plot_gc_POWER,
  width = 14,
  height = 7.8,
  device = "eps",
  dpi = 800
)


# Example -----------------------------------------------------------------

## Figure 6 - NTP bioassay methyleugenol ----------------------------------

load(file = "..\\Example\\data\\miceF.rda")
miceF$dose <- miceF$group
miceF$dose[miceF$dose == 0] <- 0
miceF$dose[miceF$dose == 1] <- 37
miceF$dose[miceF$dose == 2] <- 75
miceF$dose[miceF$dose == 3] <- 150

cst <- colSums(miceF[, 4:92])
tt5 <- as.data.frame(miceF[, names(cst[cst > 5])])
tt5$id <- seq(1, 200, 1)
tt5$dose <- miceF$dose
long <- melt(tt5, id.vars = c("dose", "id"))
long$value <- as.factor(long$value)

newdf <- rbind(head(tt5, 5), tail(tt5, 4))
knitr::kable(newdf)

plot_ntpexample <- ggplot(long, aes(x = id, y = variable)) +
  geom_tile(aes(fill = value)) +
  theme_bw(base_size = 15.3) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    #     strip.background = element_rect(fill="lightblue"),
    # axis.text=element_text(size=16),
    # axis.title.y = element_text(size = 16),
    # strip.text = element_text(size=16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(
      fill = "white",
      colour = "grey",
      size = 0.5
    ),
    axis.line.y = element_line(colour = "black")
  ) +
  facet_grid( ~ dose, scale = "free_x", labeller = label_both) +
  ylab("Tumour") +
  scale_fill_manual(values = c("white", "black"))

#plot_ntpexample

ggsave(
  ".\\Budig_Figure_6.eps",
  plot = plot_ntpexample,
  width = 12,
  height = 4.5,
  device = "eps",
  dpi = 800
)

# Appendix ----------------------------------------------------------------

## Figure 7 - Park Correlation --------------------------------------------

bin_park_cor <-
  read.csv(
    "..\\Simulations\\Sim_Park_cor_bin_dat\\intermediate_results\\bin_cor_endpoints.csv"
  )
bin_park_cor <-
  bin_park_cor %>% pivot_longer(cols = starts_with("cor"),
                                names_to = "type",
                                values_to =  "corr")


pi_label <- c(expression(paste(pi, ": ", "0.1")),
              expression(paste(pi, ": ", "0.5")))

N_label <- c(expression(paste("N", ": ", "100")),
             expression(paste("N", ": ", "200")),
             expression(paste("N", ": ", "400")))

bin_park_cor$pi0 <- factor(bin_park_cor$pi0, labels = pi_label)
bin_park_cor$ntrt <- factor(bin_park_cor$ntrt, labels = N_label)

bin_park_cor <- bin_park_cor %>%
  filter(type != "cor_s_s1s2" &
           type != "cor_s_s1s3" & type != "cor_s_s2s3")
bin_park_cor <- bin_park_cor %>%
  mutate(
    newtype = case_when(
      type == "cor_m_s1s2" ~ "MMM_j1j2",
      type == "cor_m_s1s3" ~ "MMM_j1j3",
      type == "cor_m_s2s3" ~ "MMM_j2j3",
      type == "cor_p_s1s2" ~ "Pearson_j1j2",
      type == "cor_p_s1s3" ~ "Pearson_j1j3",
      type == "cor_p_s2s3" ~ "Pearson_j2j3"
    )
  )

plot_park_cor <- ggplot(bin_park_cor, aes(x = newtype, y = corr)) +
  geom_boxplot() +
  #geom_jitter(width = 0.05, height = 0)+
  facet_grid(ntrt ~ pi0, labeller = "label_parsed") +
  theme_bw(base_size = 15.4) +
  xlab("Correlation measure and endpoint combination") +
  ylab("Correlation") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  geom_hline(yintercept = 0.9) +
  geom_hline(yintercept = 0.5)

#plot_park_cor

ggsave(
  ".\\Budig_Figure_7.eps",
  plot = plot_park_cor,
  width = 12,
  height = 6.5,
  device = "eps",
  dpi = 800
)

