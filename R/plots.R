source("R/loader.R")
library(ggplot2)

plot_analysis <- function(df) {
    p <- ggplot(df, aes(x = ISI, y = SRA, group = Label, color = Label)) +
        geom_line() +
        geom_point(aes(x = ISI, y = SRA, group = Label)) +
        geom_hline(yintercept = 1, linetype = "dashed")
    p
}
