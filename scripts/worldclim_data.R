
# Install the packages to ensure they are up to date
install.packages(c("raster", "geodata"))

# Load the packages
pacman::p_load(char = c("raster", "geodata"))

# Countries for which to get climatic data
countries <- c("DEU", "PRT", "ITA", "BEL", "GRC", "MLT")


countries_data <- countries |> 
  # Iterate through countries
  lapply(function(country) {
    # Report country
    print(country)
    # Get administrative boundaries for the country
    gadm_country <- gadm(country, path = here::here("temp_files"), level = 0)
    
    # Define climatic variables to get
    clim_vars <- c("tavg"
                   , "prec")
    
    country_data <- clim_vars |> 
      # Iterate through climatic variables
      lapply(function(var) {
        # report variable
        print(var)
        # Download climate data from WorldClim version 2.1
        worldclim_country(country = country
                          , var = var
                          , path = here::here("temp_files")
                          , res = 10
                          , version = 2.1) |> 
          # Extract values
          values() |> 
          # Trasnsform into a tibble
          as_tibble() |> 
          # Set the column names by month
          set_colnames(month.name) |> 
          # Month and values as columns
          pivot_longer(everything(), names_to = "month") |> 
          # Calculate mean
          summarise(mean(value, na.rm = T)) |> 
          # Name as variable
          set_colnames(var)
      })
    print("merging")
    # Merge variable values of the country in a data frame
    country_data <- country_data |> 
      reduce(merge
             , all = T
             , sort = F) |> 
      # Set up more informative variable names
      set_colnames(c("mean_temp_C"
                     , "precipitation_mm")) |> 
      # Add country name
      mutate(country = country_codes() |> 
               filter(ISO3 == country) |> 
               pull(NAME_ISO))
    
  }) |> 
  # Merge data of countries in a single data frame
  reduce(merge
         , all = T
         , sort = F) |> 
  as.data.frame()

# Separate temperature and precipitation data in independent data frames
temperature <- countries_data |>
  dplyr::select(-precipitation_mm)

precipitation <- countries_data |> 
  dplyr::select(-mean_temp_C)
  
# Store data frames
save(list = c("temperature", "precipitation"), file = here::here("worldclim_tables.Rdata"))

rm("countries", "clim_vars", "countries_data", "temperature", "precipitation")
