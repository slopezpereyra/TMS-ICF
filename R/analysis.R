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
        print(isi)
        if (nrow(df) == 120) {
            ind <- ind + 1
        }
        d <- df[ind, ] %>% pull(EMGPeakToPeak)
        print(d)
        avg_emg <- mean(d)
        print(avg_emg)
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
        s1 <- df %>% analyze_subject(subject, 1)
        s2 <- df %>% analyze_subject(subject, 2)
        result <- result %>%
            bind_rows(s1) %>%
            bind_rows(s2)
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
    a <- group_type_analysis(df, "HC BL")
    b <- group_type_analysis(df, "HC SWD")
    c <- group_type_analysis(df, "MDD BL")
    d <- group_type_analysis(df, "MDD SWD")
    result <- a %>%
        bind_rows(b) %>%
        bind_rows(c) %>%
        bind_rows(d)
    return(result)
}
df <- load_data() %>% clean_data()
an <- analysis(df)
an

gan <- group_type_analysis(an, "HC BL")
gan
View(gan)
final <- full_group_analysis(an)
View(final)


plot_analysis <- function(df) {
    p <- ggplot(df, aes(x = ISI, y = RA, group = Label, color = Label)) +
        geom_line()
    p
}


plot_analysis(final)

a
a$Type
