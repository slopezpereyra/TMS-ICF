library(tidyr)
library(readr)
library(dplyr)
library(tibble)

# GLOBALS
# -1 = test pulse, 0 conditioned pulse
ISI_ORDER_120 <- c(
    10, 0, 8, 0, -1, 10, 0, -1, 20, 0, -1, 15, 0, 10,
    0, 20, 0, 8, 0, 10, 0, 15, 0, -1, 8, 0, -1, 10,
    0, -1, 5, 0, -1, 5, 0, 10, 0, 8, 0, 4, 0, 8, 0,
    -1, 4, 0, -1, 4, 0, 5, 0, -1, 20, 0, -1, -1, -1,
    5, 0, 5, 0, 10, 0, -1, -1, -1, 4, 0, -1, 20, 0, 8,
    0, 8, 0, 5, 0, -1, 15, 0, -1, 15, 0, -1, 5, 0, 15,
    0, -1, 15, 0, 15, 0, 4, 0, 8, 0, 20, 0, 4, 0, 20,
    0, 4, 0, -1, -1, 10, 0, 20, 0, -1, 15, 0, 5, 0, 20, 0, 4, 0
)

ISI_ORDER_72 <- c(
    10, 8, -1, 10, -1, 20, -1, 15, 10, 20, 8, 10, 15, -1, 8, -1, 10, -1,
    5, -1, 5, 10, 8, 4, 8, -1, 4, -1, 4, 5, -1, 20, -1, -1, -1, 5, 5,
    10, -1, -1, -1, 4, -1, 20, 8, 8, 5, -1, 15, -1, 15, -1, 5, 15, -1,
    15, 15, 4, 8, 20, 4, 20, 4, -1, -1, 10, 20, -1, 15, 5, 20, 4
)

WD <- getwd()
DATA_DIR <- paste(WD, "Data/", sep = "/")
PARTS_INFO <- read_csv("TMS_ptinfo.csv")
PARTS_INFO

#' Use .txt file name to find subject number and return it as numeric.
#'
#' @param x String with file name
#' @return Numeric value representing subject ID.
get_subject_from_fn <- function(filename) {
    filename_no_letters <- gsub("[^0-9.-]", "", filename) # Remove all letters
    str_length <- nchar(filename_no_letters)
    subject_number <- substring(filename_no_letters, 1, str_length - 2) # Remove session number and dot
    return(as.numeric(subject_number))
}

get_session_from_fn <- function(filename) {
    filename_no_letters <- gsub("[^0-9.-]", "", filename) # Remove all letters
    str_length <- nchar(filename_no_letters)
    subject_session <- substring(filename_no_letters, str_length - 1, str_length) # Remove session number and dot
    return(as.numeric(subject_session))
}

get_subj_group <- function(subject) {
    g <- PARTS_INFO$group[which(PARTS_INFO$patientid == subject)]
    return(g)
}

get_session_type <- function(subject, session) {
    if (session == 1) {
        type <- PARTS_INFO$session1[which(PARTS_INFO$patientid == subject)]
    } else if (session == 2) {
        type <- PARTS_INFO$session2[which(PARTS_INFO$patientid == subject)]
    }
    return(type)
}

get_session_as_number <- function(session_name) {
    case_when(
        session_name == "Session 1" ~ 1,
        session_name == "Session 2" ~ 2,
        session_name == "Session 3" ~ 2
    )
}

# Horrible but serves its purpose and will
# only be called once in human history *cringes*
get_label <- function(df) {
    if ((df$Group[1] == 1) && (df$Type == "SWD")) {
        return("HC SWD")
    }
    if ((df$Group[1] == 1) && (df$Type == "BL")) {
        return("HC BL")
    }
    if ((df$Group[1] == 2) && (df$Type == "SWD")) {
        return("MDD SWD")
    }
    if ((df$Group[1] == 2) && (df$Type == "BL")) {
        return("MDD BL")
    }
}

load_data <- function() {
    df <- tibble()
    for (file in list.files(DATA_DIR)) {
        print(file)
        path <- paste(DATA_DIR, file, sep = "/")

        t <- read_delim(path,
            delim = "\t", comment = "#", na = c("", "NA", "(null)"),
            col_types = cols_only(
                `Sample Name` = col_character(),
                `Session Name` = col_character(),
                `EMG Peak-to-peak 1` = col_number()
            )
        )
        ind <- grep("ICF", t$`Sample Name`)
        t <- t[ind[1]:ind[2], ]
        session <- ifelse(t$`Session Name`[1] == "Session 1", 1, 2)
        subject <- get_subject_from_fn(file)
        t <- t %>%
            add_column(Subject = rep(subject)) %>%
            add_column(Type = rep(get_session_type(subject, session))) %>%
            add_column(Group = rep(get_subj_group(subject))) %>%
            add_column(ISI = create_isi_col(t))

        if (is.null(t$ISI)) {
            next
        }
        t <- t %>%
            add_column(Label = rep(get_label(t)))
        df <- bind_rows(df, t)
    }
    colnames(df) <- c(
        "Sample", "Session", "EMGPeakToPeak",
        "Subject", "Type", "Group", "ISI", "Label"
    )
    return(df)
}

clean_data <- function(df) {
    # Fix Sam's ICF coding...
    start_patterns <- c("First", "FIRST", "first", "start", "Start", "START")
    end_patterns <- c("END", "Last", "LAST", "last", "End")
    df$Sample <- df$Sample %>%
        recode_by_patterns("ICF", start_patterns, "ICF Start") %>%
        recode_by_patterns("ICF", end_patterns, "ICF End") %>%
        recode_by_patterns("MEP", start_patterns, "MEP Start") %>%
        recode_by_patterns("MEP", end_patterns, "MEP End")

    df <- df %>% mutate(Session = get_session_as_number(Session))
    return(df)
}

recode_by_patterns <- function(vec, pattern, subpatterns, new_value) {
    # Find index on elements containing pattern
    pat_indexes <- which(vec %in% grep(pattern, vec, value = TRUE))
    # print(indexes)# Correct index for each value with ICF
    for (subp in subpatterns)
    {
        subpat_indexes <- which(vec %in% grep(subp, vec[pat_indexes], value = TRUE))
        vec <- replace(vec, subpat_indexes, new_value)
    }
    return(vec)
}

#' Wrapper
get_subject_data <- function(df, subject) {
    df <- df %>% subset(Subject == subject)
    return(df)
}

create_isi_col <- function(df) {
    len <- nrow(df)
    vec <- c()
    if (len == 120) # If this ICF has 120 pulses.
        {
            vec <- ISI_ORDER_120
        } else if (len == 72) # If this ICF has 72 pulses.
        {
            vec <- ISI_ORDER_72
        } else if (len == 118) {
        vec <- ISI_ORDER_120[1:118] # Special case (see Elly's mails)
    } else {
        warning(paste("ICF INTEGRITY WARNING: ", len, " on subject ", df$Subject[1]))
        return
    }
    return(vec)
}

integrity_test <- function(df) {
    wrongs <- c()
    for (subject in unique(df$Subject)) {
        t <- subset(df, Subject == subject)
        if (nrow(t) != 240 && nrow(t) != 144 && nrow(t) != 192) {
            wrongs <- append(wrongs, subject)
        }
    }
    return(wrongs)
}

#df <- load_data() %>% clean_data()
#write_csv(df, "df.csv")
#View(df)
#t <- read_delim("SWIP_018_S1.txt",
#    delim = "\t", comment = "#", na = c("", "NA", "(null)"),
#    col_types = cols_only(
#        `Sample Name` = col_character(),
#        `Session Name` = col_character(),
#        `EMG Peak-to-peak 1` = col_number()
#    )
#)
#
#View(t)
#View(get_subject_data(df, 9))
#df
