## import libraries
library(tidyverse)
library(janitor)
library(naniar)

## data management
(
    ## read in dataframe
    read_csv(
        './Project0/Data/Project0_Clean_v2.csv',
        show_col_types = FALSE
    )
    
    ## clean names
    %>% clean_names()
    
    ## create time since wake columns
    %>% group_by(subject_id, daynumb)
    %>% arrange(collection_sample, .by_group = TRUE)
    %>% fill(sleep_diary_reported_wake_time, .direction = "down")
    %>% ungroup()
    %>% mutate(
        .after = sleep_diary_reported_wake_time,
        mems_time_since_wake =
            as.numeric(me_ms_clock_time - sleep_diary_reported_wake_time) / 3600,
        booklet_time_since_wake =
            as.numeric(booket_clock_time - sleep_diary_reported_wake_time) / 3600,
    )
    
    ## create 15 and 30 min window indicator columns
    %>% mutate(
        mems_adherence_15 = case_when(
            collection_sample %in% c(1, 3) ~ NA_real_,
            collection_sample == 2 ~ if_else(abs(mems_time_since_wake - 30/60) <= 7.5/60, 1, 0),
            collection_sample == 4 ~ if_else(abs(mems_time_since_wake - 10) <= 7.5/60, 1, 0)
        ),
        mems_adherence_30 = case_when(
            collection_sample %in% c(1, 3) ~ NA_real_,
            collection_sample == 2 ~ if_else(abs(mems_time_since_wake - 30/60) <= 15/60, 1, 0),
            collection_sample == 4 ~ if_else(abs(mems_time_since_wake - 10) <= 15/60, 1, 0)
        ),
        booklet_adherence_15 = case_when(
            collection_sample %in% c(1, 3) ~ NA_real_,
            collection_sample == 2 ~ if_else(abs(booklet_time_since_wake - 30/60) <= 7.5/60, 1, 0),
            collection_sample == 4 ~ if_else(abs(booklet_time_since_wake - 10) <= 7.5/60, 1, 0)
        ),
        booklet_adherence_30 = case_when(
            collection_sample %in% c(1, 3) ~ NA_real_,
            collection_sample == 2 ~ if_else(abs(booklet_time_since_wake - 30/60) <= 15/60, 1, 0),
            collection_sample == 4 ~ if_else(abs(booklet_time_since_wake - 10) <= 15/60, 1, 0)
        )
    )
) -> cort

