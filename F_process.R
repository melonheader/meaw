#
# Auxiliary functions to process electrophysiology readings from MEA plates
# 30.11.2022; AB
#

set.seed(42)

library(tidyverse)
library(purrr)
library(janitor)



parse_meta <- function(path_to) {
  ##
  mea_tab <- read_lines(path_to)
  ## extract well annotation data
  info.sep <- which(str_detect(mea_tab, "Well Information"))
  well.ann_slice <- mea_tab[c(info.sep + 1, info.sep + 5)]
  purrr::map(well.ann_slice,
             ~ str_split_fixed(.x, ",", Inf) %>% as.character()) %>% 
    purrr::reduce(., cbind) %>% 
    `colnames<-`(c("Well", "Treatment")) %>% 
    as_tibble() -> well.ann_raw
  return(well.ann_raw)
}

parse_mea <- function(path_to_spike_list) {
  ##
  mea_tab <- read_lines(path_to_spike_list)
  ## extract well annotation data
  info.sep <- which(str_detect(mea_tab, "Well Information"))
  well.ann_slice <- mea_tab[c(info.sep + 1, info.sep + 5)]
  purrr::map(well.ann_slice,
             ~ str_split_fixed(.x, ",", Inf) %>% as.character()) %>% 
    purrr::reduce(., cbind) %>% 
    `colnames<-`(c("Well", "Treatment")) %>% 
    as_tibble() -> well.ann_raw
  ##
  mea_tab <- read_csv(path_to_spike_list, show_col_types = F)
  ## extract main info
  tab.proc <- mea_tab[1:(info.sep - 1), -c(1:2)] %>% 
    mutate(across(!contains("lectrode"), ~ as.numeric(.x)),
           Well = str_split_fixed(Electrode, "_", 2)[, 1],
           Electrode = str_split_fixed(Electrode, "_", 2)[, 2]) %>% 
    inner_join(well.ann_raw, ., by = "Well") %>% 
    dplyr::select(Treatment, Well, Electrode, everything())
  ## out
  return(tab.proc)
}
parse_file.list <- function(wd.path, sep = ")_") {
  tibble(file_name = list.files(wd.path, pattern = ".csv")) %>% 
    mutate(data_type = str_remove(str_split_fixed(file_name, fixed(sep), Inf)[, 2], ".csv"),
           div = as.numeric(str_sub(str_split_fixed(file_name, "_DIV", 2)[, 2], 1, 2)),
           date_measured = str_split_fixed(file_name, "_", 3)[, 2]) %>% 
    dplyr::select(file_name, date_measured, div, data_type)
}


# 1 construct a retarded object with new characteristics
construct_s <- function(l.eltd, beg = NULL, end = NULL, min_spikes = 5, ids = NULL, bin_interval = 5) {
  
  # compute global time range (all electrodes)
  spike_range <- range(unlist(l.eltd))
  
  # if time borders are not provided, use the natural borders within the well
  if (is.null(end)) {
    end <- spike_range[2]
  } else {
    ## filtering maximum time
    l.eltd <- purrr::map2(l.eltd,
                          end,
                          function(v.eltd, time.max) {
                            v.eltd[v.eltd < time.max]
                          }
    )
  }
  
  if (is.null(beg)) {
    beg <- spike_range[1]
  } else {
    ## filtering minimum time
    l.eltd <- purrr::map2(l.eltd,
                          beg,
                          function(v.eltd, time.min) {
                            v.eltd[v.eltd > time.min]
                          }
    )
  }
  
  ## Remove any channels that has zero or less than a threshold of spikes, which can happen if
  ## the beg, end range is too narrow, or if a datafile is empty
  pass_index <- which(sapply(l.eltd, length) >= min_spikes)
  l.eltd <- l.eltd[pass_index]
  ## check if there are ny electrodes left active after filtering
  if (length(l.eltd) == 0) {
    return(NULL)
  }
  
  if (!is.null(ids)) {
    l.eltd <- l.eltd[ids]
  }
  
  ## Estimate the average number of spikes over the whole measurement period
  nspikes <- sapply(l.eltd, length)
  meanfiringrate <- nspikes / (end - beg)
  
  ## here comes another function
  frates <- estimate_frate(l.eltd, beg, end, time_interval = bin_interval)
  
  ## Construct output
  res <- list(channels = names(l.eltd), 
              spikes = l.eltd, 
              nspikes = nspikes,
              NCells = length(l.eltd), 
              meanfiringrate = meanfiringrate,
              firing_rates = frates,
              rec_time = end - beg)
  
  ## STTC correlation index can be implemented later
  # class(res) <- "spike.list"
  # if (length(corr_breaks) == 1) {
  #   res$corr <- NULL
  # } else {
  #   res$corr <- .corr_index(res, corr_breaks)
  # }
  return(res)
}

# -----------------------------------------------------------------------------
estimate_frate <- function(l.eltd_filt,
                           beg, end,
                           time_interval = 1, # time bin of 1sec.
                           clip = FALSE,
                           frate_min = 0,
                           frate_max = 20
) {
  
  if (is_null(l.eltd_filt)) {
    return(NULL)
  }
  
  time_breaks <- seq(from = beg, to = end, by = time_interval)
  nbins <- length(time_breaks)
  
  l.eltd.frates <- purrr::map(l.eltd_filt, 
                              function(eltd) {
                                purrr::map_int(eltd,
                                               function(spike) {
                                                 ## add count into the bin
                                                 as.integer((spike - time_breaks[1]) / time_interval)
                                               }
                                ) -> cts.raw
                                cts <- table(cts.raw)
                                ##
                                bins <- numeric()
                                bins <- rep(0, length = nbins + 1)
                                names(bins) <- 0:nbins
                                ## update counts
                                bins[names(cts)] <- unname(cts)
                                ##
                                return(bins)
                              }
  )
  
  ## Clip frequencies if option selected
  if (clip) {
    purrr::map(l.eltd.frates,
               ~ pmin(pmax(.x, frate_min), frate_max)) -> l.eltd.frates
  }
  
  ## compuate average frequency per bin
  ### First turn results into a table for brevitys
  purrr::reduce(l.eltd.frates, cbind) %>% as.data.frame() -> t.eltd.frates
  colnames(t.eltd.frates) <- names(l.eltd.frates)
  ### compute averages
  eltd.frates_avrgs <- apply(t.eltd.frates, 1, mean)
  
  ## prepare output
  res <- list(frate = t.eltd.frates,
              frate.avg = eltd.frates_avrgs,
              time_bins = time_breaks,
              time_interval = time_interval)
  return(res)
}


# -----------------------------------------------------------------------------
filter_spike.list <- function(s_list, 
                              eltd_min_rate = (1 / 60),
                              eltd_max_rate = 25,
                              well_min_rate = 4) {
  ## 
  if (is_null(s_list)) {
    return(NULL)
  }
  ## check for the number of active electrodes in a well and discard if too little
  if (s_list$NCells < well_min_rate) {
    return(NULL)
  }
  
  ## 
  low <- which(s_list$meanfiringrate < eltd_min_rate)
  high <- which(s_list$meanfiringrate > eltd_max_rate)
  extrms <- c(low, high)
  ##drop electrodes from the object
  if (!is_empty(extrms)) {
    s_list$channels <- s_list$channels[!s_list$channels %in% extrms]
    ## to keep
    eltd_keep <- s_list$channels
    s_list$spikes <- s_list$spikes[names(s_list$spikes) %in% eltd_keep]
    s_list$nspikes <- s_list$nspikes[names(s_list$nspikes) %in% eltd_keep]
    s_list$NCells <- length(s_list$channels)
    s_list$meanfiringrate <- s_list$meanfiringrate[names(s_list$meanfiringrate) %in% eltd_keep]
    
    s_list$firing_rates$frate <- s_list$firing_rates$frate[, eltd_keep]
    s_list$firing_rates$frate.avg <- apply(s_list$firing_rates$frate, 1, mean)
  }
  ##
  return(s_list)
}