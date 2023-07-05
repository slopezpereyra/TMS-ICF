source("r/loader.r")
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(tibble)
library(robustbase)
library(lmPerm)
library(MASS)

df <- read_csv("JLRA_DF.csv")
df <- filter(df, ISI != 0)
df$ISI <- as.factor(df$ISI)
df$Group  <- as.factor(df$Group)

plot_analysis <- function(df) {
    p <- ggplot(df, aes(x = isi, y = sra, group = label, color = label)) +
        geom_line() +
        geom_point(aes(x = isi, y = sra, group = label)) +
        geom_hline(yintercept = 1, linetype = "dashed")
    p
}



# ------------ 

HCSWD <- df %>% filter(Label == "HC SWD") %>% filter(Subject != 8)
HCBL <- df %>% filter(Label == "HC BL")
MDDSWD <- df %>% filter(Label == "MDD SWD") %>% filter(Subject != 8)
MDDBL <- df %>% filter(Label == "MDD BL")

# Is median before SWD greater? 
wilcox.test(HCBL$EMGPeakToPeak, HCSWD$EMGPeakToPeak, paired=TRUE, 
            alternative="greater")
# Is median after SWD lower?
wilcox.test(MDDBL$EMGPeakToPeak, MDDSWD$EMGPeakToPeak, paired=TRUE, 
            alternative="less")

wilcox.test(HCBL$EMGPeakToPeak, MDDBL$EMGPeakToPeak,
            alternative="greater")

wilcox.test(HCBL$RRA, HCSWD$RRA, paired=TRUE)
# Is median after SWD lower?
wilcox.test(MDDBL$RRA, MDDSWD$RRA, paired=TRUE, alternative="less")

wilcox.test(HCSWD$RRA, MDDSWD$RRA)

library(coin)

df
X <- filter(df, Group == 1) %>% filter(Subject != 8)

median_test(RRA~as.factor(Label), data=X)
oneway_test(RRA~as.factor(Label),data=X)

# Calculate mean and standard deviation by group using dplyr
summary_stats <- df %>%
  group_by(Label) %>%
  summarize(mean = mean(EMGPeakToPeak), sd = sd(EMGPeakToPeak), 
            median = median(EMGPeakToPeak))

print(summary_stats)

library(effectsize)
X <- filter(df, Group == 2)
X
cohens_d(RRA ~ Type, data = X)


mean_values <- aggregate(RRA ~ Label + ISI, 
                         data = df %>% filter(Group == 1), 
                         FUN = mean)
ggplot(mean_values, aes(x = ISI, y = RRA, fill=Label)) + 
    geom_bar(stat="identity", position="dodge") + 
    labs(x = "ISI", y = "RRA: HC BL vs MDD SWD")


mean_values <- aggregate(RRA ~ Label, 
                         data = df %>% filter(Group == 2), 
                         FUN = mean)
ggplot(mean_values, aes(x = Label, y = RRA, fill=Label)) +
        geom_bar(stat="identity", position="dodge") + 
        labs(x = "Group", y = "Avg. EMG Peak to peak: MDDs Overall")





























