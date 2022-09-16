##' Download Targets
##' @return data.frame in long format with days as rows, and time, site_id, variable, and observed as columns
download_targets <- function(){
  readr::read_csv("https://data.ecoforecast.org/neon4cast-targets/aquatics/aquatics-targets.csv.gz", guess_max = 1e6)
}

##' Download Site metadata
##' @return metadata dataframe
download_site_meta <- function(){
  site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") 
  site_data %>% filter(as.integer(aquatics) == 1)
}


##' append historical meteorological data into target file
##' @param target targets dataframe
##' @return updated targets dataframe with added weather data
merge_met_past <- function(target){
  
  ## connect to data
  df_past <- neon4cast::noaa_stage3()
  
  ## filter for site and variable
  sites <- unique(target$site_id)
  noaa_past <- df_past %>%
    dplyr::filter(site_id %in% sites,
                  variable == "air_temperature") %>%
    dplyr::collect()
  
  ## Aggregate (to day) and convert units of drivers
  noaa_past_mean <- noaa_past %>% 
    mutate(date = lubridate::as_date(time)) %>% 
    group_by(date, site_id) %>% 
    summarize(air_temperature = mean(predicted, na.rm = TRUE), .groups = "drop") %>% 
    rename(datetime = date) %>% 
    mutate(air_temperature = air_temperature - 273.15)
  
  ## Merge in past NOAA data into the targets file, matching by date.
  target <- target %>% 
    select(datetime, site_id, variable, observed) %>%
    filter(variable %in% c("temperature", "oxygen")) %>%
    tidyr::pivot_wider(names_from = "variable", values_from = "observed")
  
  target <- left_join(target, noaa_past_mean, by = c("datetime","site_id"))
  
}

##' Download NOAA GEFS weather forecast
##' @param forecast_date start date of forecast
##' @return dataframe
download_met_forecast <- function(forecast_date){
  ## connect to data
  df_future <- neon4cast::noaa_stage2()
  
  noaa_date <- forecast_date - lubridate::days(1)  #Need to use yesterday's NOAA forecast because today's is not available yet
  
  ## filter available forecasts by date and variable
  met_future <- df_future %>%
    dplyr::filter(cycle == 0,
                  start_date == as.character(noaa_date),
                  time >= lubridate::as_datetime(forecast_date), 
                  variable == "air_temperature") %>%
    dplyr::collect()
  
  ## aggregate to daily
  met_future <- met_future %>% 
    mutate(time = lubridate::as_date(time)) %>% 
    group_by(time, site_id, ensemble) %>%
    summarize(air_temperature = mean(predicted), .groups = "drop") %>% 
    mutate(air_temperature = air_temperature - 273.15) %>%
    select(time, site_id, air_temperature, ensemble)
  
  return(met_future)
}