source("R/loader.R")
library(ggplot2)

analyze_subject <- function(df, subject, session) {
    df <- df %>%
        get_subject_data(subject) %>%
        subset(Session == session)
    t <- df %>%
        subset(ISI == -1) %>%
        pull(EMGPeakToPeak)
    d_df <- df %>%
        subset(ISI != -1) %>%
        subset(ISI != 0)
    dmeans <- c()
    for (isi in (unique(sort(d_df$ISI)))) {
        ind <- which(df$ISI == isi)
        if (nrow(df) == 120) {
            ind <- ind + 1
        }
        d <- df[ind, ] %>% pull(EMGPeakToPeak)
        avg_emg <- mean(d)
        dmeans <- append(dmeans, avg_emg)
    }
    ras <- dmeans / mean(t)
    n <- length(ras)
    result <- tibble(
        Subject = rep(subject, n),
        Label = rep(df$Label[1]),
        ISI = unique(sort(d_df$ISI)),
        SRA = ras
    )
    return(result)
}

analysis <- function(df) {
    result <- tibble()
    for (subject in unique(df$Subject)) {
        result <- result %>% bind_rows(analyze_subject(df, subject, 1), 
                                        analyze_subject(df, subject, 2))
    }
    return(result)
}

group_type_analysis <- function(df, label) {
    df <- df %>% subset(Label == label)
    means <- c()
    for (isi in unique(sort(df$ISI))) {
        isi_df <- df %>% subset(ISI == isi)
        means <- append(means, mean(isi_df$SRA))
    }
    n <- length(unique(isi))
    result <- tibble(
        Label = rep(label, n),
        ISI = unique(sort(df$ISI)),
        RA = means
    )
    return(result)
}

full_group_analysis <- function(df) {
    result <- bind_rows(group_type_analysis(df, "HC BL"),
        group_type_analysis(df, "HC SWD"), 
        group_type_analysis(df, "MDD BL"),
        group_type_analysis(df, "MDD SWD"))

    return(result)
}

plot_analysis <- function(df) {
    p <- ggplot(df, aes(x = ISI, y = RA, group = Label, color = Label)) +
        geom_line() +
        geom_point(aes(x=ISI, y=RA, group=Label))
    p
}


