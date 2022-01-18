library(tidyverse)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)

#### n < 50 (Ecos, Hotspots, Regions) ####
setwd("")
files <- list.files(getwd(), pattern = '.csv', full.names = F)

df <- read_csv("countries-stats.csv") %>%
  select(-Region) %>%
  data.frame()
df[is.na(df)] <- 0

## one-sample Shapiro-Wilk Normality Test
for (j in names(df)){
  sink("Shapiro-Wilk.txt", append = T)
  cat("=============================\n")
  cat(paste(j, "\n"))
  cat("=============================\n")

  mww <- capture.output(shapiro.test(df[[j]]))

  cat(mww, file = "Shapiro-Wilk.txt", sep = "\n", append = T)
  sink()
}

## Mann-Whitney-Wilcox Non-Parametric Test
def <- c("m","f","d","a")
ch <- "_ch"
wc <- "_wc"

for (j in unique(def)){
  sink("Mann-Whitney-Wilcox.txt", append = T)
  cat("=============================\n")
  cat(paste(j, ch, ', ', j, wc, sep = '', "\n"))
  cat("=============================\n")

  x <- df[[paste(j, ch, sep = '')]]
  y <- df[[paste(j, wc, sep = '')]]

  t <- capture.output(wilcox.test(x, y))

  cat(t,file="Mann-Whitney-Wilcox.txt", sep = "\n", append = T)
  sink()
}

## Plots
df <- df %>%
  gather("m_ch", "f_ch", "d_ch", "a_ch",
         "m_wc", "f_wc", "d_wc", "a_wc",
         key = "text",
         value = "value") %>%
  mutate(value = round(as.numeric(value),0))

# Bar charts
lab <- c("ML (CH)", "FAO (CH)", "Dry (CH)", "AI (CH)",
         "ML (WC)", "FAO (WC)", "Dry (WC)", "AI (WC)")
df$text <- lab

ggplot(df, aes(df$text, df$value, fill = df$text)) +
  geom_bar(stat="identity") +
  theme_ipsum() +
  ylab("Extent (Km^2)") +
  xlab("Definition (Dataset)") +
  theme(
    axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    panel.background = element_rect(fill = "ghostwhite"),
    panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
    panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))

# Boxplot
ggplot(df, aes(df$text, df$value, fill = df$text)) +
  geom_boxplot() +
  theme_ipsum() +
  ylab("Extent (Km^2)") +
  xlab("Definition (Dataset)") +
  theme(
    axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    panel.background = element_rect(fill = "ghostwhite"),
    panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
    panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))

# Density Plots
lab <- c("Murphy & Lugo (CHELSA)", "FAO (CHELSA)", "Dryflor (CHELSA)", "Aridity (CHELSA)",
         "Murphy & Lugo (Worldclim)", "FAO (Worldclim)", "Dryflor (Worldclim)", "Aridity (Worldclim)")
df$text <- lab

ggplot(df, aes(df$value, group = df$text, fill = df$text)) +
  geom_density() +
  theme_ipsum() +
  xlab("Extent (Km^2)") +
  ylab("Density") +
  facet_wrap(~df$text, nrow = 2) +
  theme(
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    panel.background = element_rect(fill = "ghostwhite"),
    panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
    panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"),
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank())

# Stacked Density Plot
ggplot(df, aes(df$value, group = df$text, fill = df$text)) +
  geom_density(alpha = 0.4) +
  theme_ipsum() +
  xlab("Extent (Km^2)") +
  ylab("Density") +
  labs(fill = "Definition (Dataset)") +
  theme(
    axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(color = "darkblue", size = 16),
    legend.background = element_rect(fill = "white"),
    legend.position = c(0.8, 0.7),
    panel.background = element_rect(fill = "ghostwhite"),
    panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
    panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))

# Violin plot
ggplot(df, aes(df$text, df$value)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Distribution of Tropical Dry Forest Extent at the Macro-Regional Level",
       x = "Definition (Dataset)",
       y = "Extent (Km^2)") +
  scale_fill_brewer(palette = "Spectral") +
  theme_minimal()


# #### n > 50 (Countries) ####
# ## One-sample Kolmogorov-Smirnov Normality Test
# for (j in names(df)){
#   sink("ks-one_sample.txt", append = T)
#   cat("=============================\n")
#   cat(paste(j, "\n"))
#   cat("=============================\n")
#   ks1 <- capture.output(ks.test(df[[j]],
#                                 "pnorm",
#                                 mean=mean(df[[j]]),
#                                 sd=sd(df[[j]])))
#   cat(ks1, file = "ks-one_sample.txt", sep = "\n", append = T)
# 
#   # QQ1 <- ggplot(df, aes(sample = df[[j]])) +
#   #   stat_qq() +
#   #   stat_qq_line(col = "red")
#   # gQQ1 <- ggplotGrob(QQ1)
#   # 
#   # ed <- ecdf(df[[j]])
#   # maxdiffidx <- which.max(abs(ed(df[[j]])-pnorm(df[[j]])))
#   # maxdiffat <- df[[j]][maxdiffidx]
#   # 
#   # KS1 <- ggplot(df, aes(df[[j]])) +
#   #   stat_ecdf() +
#   #   stat_function(fun = pnorm, colour = "red") +
#   #   geom_vline(xintercept = maxdiffat, lty = 2)
#   # gKS1 <- ggplotGrob(KS1)
#   # 
#   # Stack <- ggplot(df, aes(df[[j]], group = df[[j]], fill = df[[j]])) +
#   #   geom_density(adjust=1.5) +
#   #   theme_ipsum()
#   # gStack <- ggplotGrob(Stack)
#   # 
#   # aspect_ratio <- 1.5
#   # 
#   # ggsave(grid::grid.draw(cbind(gQQ1, gKS1)),
#   #        file = paste("qqplots/", "-", j,
#   #                     ".png",
#   #                     sep=""),
#   #        height = 4, width = 4*aspect_ratio, scale = 2)
#   # 
#   # ggsave(gStack,
#   #        file = paste("stacked/", "-", j,
#   #                     ".png",
#   #                     sep=""),
#   #        height = 4, width = 4*aspect_ratio, scale = 2)
# 
#   sink()
# }
# 
# ## Paired T-Test
# for (i in files) {
#   df <- read_csv(i) %>%
#     select(-Region) %>%
#     data.frame()
#   df[is.na(df)] <- 0
#   sink("paired-t_test.txt", append=T)
#   def <- c("m","f","d","a")
#   ch <- "_ch"
#   wc <- "_wc"
#   for (j in unique(def)){
#     cat("=============================\n")
#     cat(paste(i,j,ch,j,wc,"\n"))
#     cat("=============================\n")
# 
#     x <- df[[paste(j,ch,sep='')]]
#     y <- df[[paste(j,wc,sep='')]]
# 
#     t <- capture.output(t.test(x, y, paired = T))
# 
#     cat(t,file="paired-t_test.txt",sep="\n",append=T)
#   }
#   sink()
# }
# 
# ## Plots
# df <- df %>%
#   gather("m_ch", "f_ch", "d_ch", "a_ch",
#          "m_wc", "f_wc", "d_wc", "a_wc",
#          key = "text",
#          value = "value") %>%
#   mutate(value = round(as.numeric(value),0))
# 
# # Bar charts
# lab <- c("ML (CH)", "FAO (CH)", "Dry (CH)", "AI (CH)",
#          "ML (WC)", "FAO (WC)", "Dry (WC)", "AI (WC)")
# df$text <- lab
# 
# ggplot(df, aes(df$text, df$value, fill = df$text)) +
#   geom_bar(stat="identity") +
#   theme_ipsum() +
#   ylab("Extent (Km^2)") +
#   xlab("Definition (Dataset)") +
#   theme(
#     axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
#     axis.title.y = element_text(size = 16),
#     axis.title.x = element_text(size = 16),
#     legend.position = "none",
#     panel.background = element_rect(fill = "ghostwhite"),
#     panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
#     panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))
# 
# # Boxplot
# ggplot(df, aes(df$text, df$value, fill = df$text)) +
#   geom_boxplot() +
#   theme_ipsum() +
#   ylab("Extent (Km^2)") +
#   xlab("Definition (Dataset)") +
#   theme(
#     axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
#     axis.title.y = element_text(size = 16),
#     axis.title.x = element_text(size = 16),
#     legend.position = "none",
#     panel.background = element_rect(fill = "ghostwhite"),
#     panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
#     panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))
# 
# # Density Plots
# lab <- c("Murphy & Lugo (CHELSA)", "FAO (CHELSA)", "Dryflor (CHELSA)", "Aridity (CHELSA)",
#          "Murphy & Lugo (Worldclim)", "FAO (Worldclim)", "Dryflor (Worldclim)", "Aridity (Worldclim)")
# df$text <- lab
# 
# ggplot(df, aes(df$value, group = df$text, fill = df$text)) +
#   geom_density() +
#   theme_ipsum() +
#   xlab("Extent (Km^2)") +
#   ylab("Density") +
#   facet_wrap(~df$text, nrow = 2) +
#   theme(
#     axis.title.y = element_text(size = 16),
#     axis.title.x = element_text(size = 16),
#     panel.background = element_rect(fill = "ghostwhite"),
#     panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
#     panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"),
#     legend.position="none",
#     panel.spacing = unit(0.1, "lines"),
#     axis.ticks.x=element_blank())
# 
# # Stacked Density Plot
# ggplot(df, aes(df$value, group = df$text, fill = df$text)) +
#   geom_density(alpha = 0.4) +
#   theme_ipsum() +
#   xlab("Extent (Km^2)") +
#   ylab("Density") +
#   labs(fill = "Definition (Dataset)") +
#   theme(
#     axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
#     axis.title.y = element_text(size = 16),
#     axis.title.x = element_text(size = 16),
#     legend.text = element_text(size = 12),
#     legend.title = element_text(color = "darkblue", size = 16),
#     legend.background = element_rect(fill = "white"),
#     legend.position = c(0.8, 0.7),
#     panel.background = element_rect(fill = "ghostwhite"),
#     panel.grid.major = element_line(size = 1 , linetype = "solid" , colour = "white"),
#     panel.grid.minor = element_line(size = 1 , linetype = "solid" , colour = "white"))