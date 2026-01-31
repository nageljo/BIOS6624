## import libraries
library(tidyverse)
library(janitor)
library(naniar)

## read in dataframe and clean names
cort <- (
    read_csv(
        './Project0/Data/Project0_Clean_v2.csv',
        show_col_types = FALSE
    )
    %>% clean_names()
)

## data management
(
    cort
    %>% select(-c(booklet_sample_interval_decimal_time_mins,
                  me_ms_sample_interval_decimal_time_mins))
    ## create time since wake columns
    %>% group_by(subject_id, collection_date)
    %>% arrange(collection_sample, .by_group = TRUE)
    %>% fill(sleep_diary_reported_wake_time, .direction = "down")
    %>% ungroup()
    %>% mutate(
        .keep = "none",
        mems_time_since_wake =
            as.numeric(me_ms_clock_time - sleep_diary_reported_wake_time) / 3600,
        booklet_time_since_wake =
            as.numeric(booket_clock_time - sleep_diary_reported_wake_time) / 3600,
    )
)

