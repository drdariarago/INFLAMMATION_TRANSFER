### Phospho normalize
## Use reference genes to normalize phospho data

library(readxl)
library(tidyverse)
library(magrittr)

# Import phospho data 

results_path <-
  here::here('data/proteomics/20210202_second_run/MaxQuant_table_all4exp_Jan_2021.xlsx')
  
results_list <- 
  results_path %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>%
  map(.f = readxl::read_xlsx, path = results_path)

results_list_with_id <-
  results_list %>% 
  map(
    .f = ~ mutate(
      .data = .x, 
      sps_id = 
        Leading.proteins %>% 
        str_extract( string = , pattern = "[[:alnum:]]*(?=_MOUSE$)" ) %>% 
        paste(
          ., 
          paste(Amino.acid, Positions.within.proteins, sep = ""),
          "", sep = ";")
    ))

# Normalise using standard mouse proteins
library(PhosR)
data('SPSs')

results_mat_list <-
  # Convert to site by sample matrices
  results_list_with_id %>% 
  map(
    ~ select(.x, sps_id, starts_with(match = "Reporter") ) %>% 
      set_colnames( gsub( x = colnames(.), pattern = "Reporter Intensity_", replacement = "")) %>% 
      column_to_rownames( var = "sps_id")
  ) %>% 
  # Filter sites with low phosphorilation
  map(
    ~ selectGrps(
      mat = .x, percent = .6, n = 1, 
      grps = str_extract( string = colnames(.x), pattern = "^[:alnum:]")
    )
  ) %>% 
  # Perform median scaling
  map(
    ~ medianScaling(
      mat = .x, scale = TRUE, reorder = FALSE,
      grps = str_extract( string = colnames(.x), pattern = "^[:alnum:]" )
    )
  )

# Annotate the position of control sites and the design matrix 

results_ctl_list <-
  results_mat_list %>% 
  map(
    .f = ~ which( rownames(.x) %in% SPSs)
  )

results_M_list <-
  results_mat_list %>% 
  map(
    ~ data.frame( 
      row.names = colnames(.x) %>% janitor::make_clean_names(string = .),
      ctrl = grepl("Control", colnames(.x)) %>% as.numeric, 
      lps = grepl("LPS", colnames(.x)) %>% as.numeric
    ) %>%
      as.matrix
  )

# Compile results into list grouped by treatment
RUV_data_list <-
  transpose(
    list(
      M = results_M_list, 
      ctl = results_ctl_list, 
      mat = results_mat_list
    )
  )

RUV_results <-
  RUV_data_list %>% 
  map(
    .f = ~ RUVphospho(
      mat = .x$mat,
      M = .x$M,
      ctl = .x$ctl, 
      k = 5
    )
  )

# Check variance within treatments
data.frame(
  ctrl_var_pre = map(RUV_data_list %>% transpose %$% mat, ~ .x[,1:4] %>% apply(., 1, var)),
  ctrl_var_post = map(RUV_results, ~ .x[,1:4] %>% apply(., 1, var) %>% summary),
  trt_var_pre = map(RUV_data_list %>% transpose %$% mat, ~ .x[,5:8] %>% apply(., 1, var)),
  trt_var_post = map(RUV_results, ~ .x[,5:8] %>% apply(., 1, var))
)

# Save RUV normalized results

write_rds(x = RUV_results, file = "results/phospho_normalize/normalized_results_list.rds")
