
# Christine



####################
# Expression score functions, rolling window
####################
# Below is the function to calculate the expression score for a rolling window of a fixed size (e.g. 10kb, in 1kb steps),
# written for the 6 gate reporter hopping
# First function is for one cell line, then a wrapper for multiple cell lines and a wrapper for strand specific data.

# The TIB with integrations needs to contain a column mTurq_expr_gate with the median expression of the sorted gate

expr_score_fixedbin_CL = function(TIB, binsize = 1000, N_bins_window = 10, 
                                  calculation_window =  c(34500000, 35000000), 
                                  CL = "23_34A", EXP = all_exp, STRAND = c("+", "-")){
  tib_for_analysis = filter(TIB, hopped & cell_line %in% CL & experiment %in% EXP)

  # #for now guestimate the fluorescence per gate, later extract this from the SB transfected FACS data (of all experiments combined?)
  # expr_per_gate_tib = tibble(population = paste0("P", 1:6),
  #                            expression = c(3E4, 1E4, 4E3, 2E3, 8E2, 3E2))
  # 
  
  #don't filter by strand, we want to normalize to the total number of integrations of the population
  #calculate total integrations
  # total_ints_ctrl = filter(tib_for_analysis, population == "ctrl") %>% nrow()
  total_ints = tib_for_analysis %>% group_by(population) %>% summarise(total_ints = n(), .groups = 'drop_last')
  # total_ints_expr = left_join(total_ints, expr_per_gate_tib, by = join_by("population")) 
  
  ## BOUNDARIES BASED ON INTERVAL SIZE
  boundaries_all = seq(from = calculation_window[1], to = calculation_window[2], by = binsize)
  boundaries_all_tib = tibble(interval_number = 1:length(boundaries_all), start_interval = boundaries_all) %>%
    arrange(start_interval) %>%
    mutate(mid_interval = (start_interval + lead(start_interval, 1))/2)
  
  # number of integrations per interval
  cut_by_intervals = tib_for_analysis %>% 
    filter(strand %in% STRAND) %>% 
    filter(chr == "chr3" & start >= calculation_window[1] & start <= calculation_window[2])  %>%
    arrange(start) %>%
    #per integration: in which 1-ctrl interval is it located?
    mutate(interval_number = findInterval(start, boundaries_all)) %>%
    #per interval: how many integrations are there in each population
    group_by(mTurq_expr_gate, population, interval_number) %>%
    summarise(N_ints = n(), .groups = 'drop_last')
  
  # calculate weighted expression factor per gate
  expr_score_by_interval_long = cut_by_intervals %>%
    left_join(total_ints, by = join_by('population' )) %>%
    mutate(N_norm = N_ints / total_ints) %>%
    mutate(expr_norm = N_norm * mTurq_expr_gate ) 
  
  # add up the gates and divide by sum of the weighing factors
  # then calculate rolling mean
  expr_score_tib = expr_score_by_interval_long %>%
    filter(population != "ctrl") %>%
    group_by(interval_number) %>%
    summarise(sum_fluo_weighted = sum(expr_norm),
              sum_weighfactors = sum(N_norm),
              total_ints_bin = sum(N_ints),
              .groups = 'drop_last') %>%
    mutate(expr_score = sum_fluo_weighted / sum_weighfactors) %>%
    right_join(boundaries_all_tib, by = "interval_number") %>%
    arrange(interval_number) %>%
    mutate(expr_score_window = RcppRoll::roll_mean(expr_score, n = N_bins_window, align = "center", fill = NA, na.rm = T),
           total_ints_window = RcppRoll::roll_sum(total_ints_bin, n = N_bins_window, align = "center", fill = 0, na.rm = T))
  
  expr_score_tib
}

expr_score_function = function(TIB, binsize = 1000, N_bins_window = 10, 
                               calculation_window =  c(34500000, 35000000), 
                               CL = "23_34A", EXP = all_exp, STRAND = c("+", "-")){
  expr_ls = lapply(CL, function(cl){
    expr_cl = expr_score_fixedbin_CL(TIB, binsize = binsize, N_bins_window = N_bins_window,
                                     calculation_window = calculation_window, 
                                     CL = cl, EXP = EXP, STRAND = STRAND)
    expr_cl
  })
  names(expr_ls) = CL
  bind_rows(expr_ls, .id = 'cell_line')
}

expr_score_function_stranded = function(TIB, binsize = 1000, N_bins_window = 10, 
                                        calculation_window =  c(34500000, 35000000), 
                                        CL = "23_34A", EXP = all_exp){
  bind_rows("+" = expr_score_function(TIB, binsize = binsize, N_bins_window = N_bins_window, calculation_window = calculation_window, 
                                      CL = CL, EXP = EXP, STRAND = "+"), 
            "-" = expr_score_function(TIB, binsize = binsize, N_bins_window = N_bins_window, calculation_window = calculation_window, 
                                      CL = CL, EXP = EXP, STRAND = "-"), 
            .id = 'strand')}

####################
# Expression score functions, arbitrary fixed bins
####################
#Give tibble with the borders of arbitrary bins (column names start_interval, end_interval)
# the arbitrary bins need to be non-overlapping
# The TIB with integrations needs to contain a column mTurq_expr_gate with the median expression of the sorted gate


expr_score_arbitrary_bins_CL = function(TIB, 
                                          bins = bins_oi, #tibble, with column names start_interval, end_interval
                                          # binsize = 1000, 
                                          # N_bins_window = 10, 
                                          # calculation_window =  c(34500000, 35000000), 
                                          CL = "23_34A", 
                                          EXP = all_exp, 
                                          STRAND = c("+", "-")){
    tib_for_analysis = filter(TIB, hopped & cell_line %in% CL & experiment %in% EXP)
    
    #don't filter by strand, we want to normalize to the total number of integrations of the population
    #calculate total integrations
    total_ints = tib_for_analysis %>% 
      group_by(population) %>%
      summarise(total_ints = n(), .groups = 'drop_last')
    
    ## BOUNDARIES BASED ON GIVEN BINS
    boundaries_all = unique(sort(c(bins$start_interval, bins$end_interval))) #get all the unique boundaries of our bins 
    boundaries_all_tib = tibble(interval_number = 1:length(boundaries_all), start_interval = boundaries_all) %>%
      arrange(start_interval) %>%
      mutate(mid_interval = (start_interval + lead(start_interval, 1))/2) %>%
      left_join(bins,by = join_by('start_interval' )) #to add names of defined intervals
    
    calculation_window = c(min(boundaries_all_tib$start_interval),
                           max(boundaries_all_tib$start_interval)) #the last interval is artificial, just to set an endpoint
    
    # number of integrations per interval
    cut_by_intervals = tib_for_analysis %>% 
      filter(strand %in% STRAND) %>% 
      filter(chr == "chr3" & start >= calculation_window[1] & start <= calculation_window[2])  %>%
      arrange(start) %>%
      #per integration: in which 1-ctrl interval is it located?
      mutate(interval_number = findInterval(start, boundaries_all)) %>%
      #per interval: how many integrations are there in each population
      group_by(mTurq_expr_gate , population, interval_number) %>% 
      dplyr::summarise(N_ints = n(), .groups = 'drop_last')
    
    # calculate weighted expression factor per gate
    expr_score_by_interval_long = cut_by_intervals %>%
      left_join(total_ints, by = join_by('population' )) %>%
      mutate(N_norm = N_ints / total_ints) %>%
      mutate(expr_norm = N_norm * mTurq_expr_gate ) 
    
    # add up the gates and divide by sum of the weighing factors
    # !don't calculate a rolling mean!
    expr_score_tib = expr_score_by_interval_long %>%
      filter(population != "ctrl") %>%
      group_by(interval_number) %>%
      dplyr::summarise(sum_fluo_weighted = sum(expr_norm),
                       sum_weighfactors = sum(N_norm), .groups = 'drop_last') %>%
      mutate(expr_score = sum_fluo_weighted / sum_weighfactors) %>%
      right_join(boundaries_all_tib, by = "interval_number") %>%
      filter(!is.na(end_interval)) #remove the extra 'interval' at the end
    
    expr_score_tib
}

#######################
#bootstrapping
######################
#HINT: run this in a markdown block where you silence messages, otherwise you get all the dplyr messages N_resamples times...
#important: the tibble needs to be ungrouped, unless you want to create the samples separately for specific groups
create_bootstrap_tibble = function(FUN, TIB, N_resamples = 200){
  # suppressMessages(
  scores_ls = lapply(1:N_resamples, function(it_i){
    resample = TIB %>%
      slice_sample(prop = 1, replace = TRUE)
    FUN(resample)
  })
  # )
  names(scores_ls) = 1:N_resamples
  scores_tib = bind_rows(scores_ls, .id = 'random_sample')
}

#turn the long tibble of N resamples into a wide tibble with summary statistics
#NB: this function expects there to be a column 'interval_number'  and 'mid_interval', 
#you may have to add those to the function calculation your (enrichment) score
bootstrap_summary_fun = function(BOOTSTIB, score_column = 'expr_score_window', conf_int = 0.75){
  lower_bound_confint = (1-conf_int)/2
  upper_bound_confint = 1-lower_bound_confint
  
  BOOTSTIB %>%
    group_by(interval_number, mid_interval) %>%
    #note: .data[[score_column]] turns the string in score_column into a column name that dplyr can use (e.g. for summarize, mutate or filter)
    summarize(expr_score_window_mean = mean(.data[[score_column]], na.rm = T),
              expr_score_window_median = median(.data[[score_column]], na.rm = T),
              expr_score_window_lower = quantile(.data[[score_column]], lower_bound_confint, na.rm = T),
              expr_score_window_upper = quantile(.data[[score_column]], upper_bound_confint, na.rm = T),
              N_samples_notNA = sum(!is.na(.data[[score_column]]))
    ) 
}

bootstrap_summary_fun_CL = function(BOOTSTIB, score_column = 'expr_score_window', conf_int = 0.75){
  lower_bound_confint = (1-conf_int)/2
  upper_bound_confint = 1-lower_bound_confint
  
  BOOTSTIB %>%
    group_by(cell_line, interval_number, mid_interval) %>%
    #note: .data[[score_column]] turns the string in score_column into a column name that dplyr can use (e.g. for summarize, mutate or filter)
    summarize(expr_score_window_mean = mean(.data[[score_column]], na.rm = T),
              expr_score_window_median = median(.data[[score_column]], na.rm = T),
              expr_score_window_lower = quantile(.data[[score_column]], lower_bound_confint, na.rm = T),
              expr_score_window_upper = quantile(.data[[score_column]], upper_bound_confint, na.rm = T),
              N_samples_notNA = sum(!is.na(.data[[score_column]]))
    ) 
}

bootstrap_summary_fun_str = function(BOOTSTIB, score_column = 'expr_score_window', conf_int = 0.75){
  lower_bound_confint = (1-conf_int)/2
  upper_bound_confint = 1-lower_bound_confint
  
  BOOTSTIB %>%
    group_by(strand, interval_number, mid_interval) %>%
    #note: .data[[score_column]] turns the string in score_column into a column name that dplyr can use (e.g. for summarize, mutate or filter)
    summarize(expr_score_window_mean = mean(.data[[score_column]], na.rm = T),
              expr_score_window_median = median(.data[[score_column]], na.rm = T),
              expr_score_window_lower = quantile(.data[[score_column]], lower_bound_confint, na.rm = T),
              expr_score_window_upper = quantile(.data[[score_column]], upper_bound_confint, na.rm = T)
    ) 
}
