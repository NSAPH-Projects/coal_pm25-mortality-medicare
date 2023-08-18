# renv::init()
# renv::install("lhenneman/disperseR@dev")

library( fst)
library( data.table)
library( parallel)
library( dplyr)
library( foreign)
library( sf)
library( raster)
library( ggplot2)

dir_data <- 'data/data' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

## ==================================================== ##
##  Read in covariates data
## ==================================================== ##
files.mort <- list.files(#"/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
  '/n/dominici_nsaph_l3/data/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/',
  pattern = "\\.fst",
  full.names = TRUE)


# For prototyping, just read in 2 years
# f <- f[1:2]

# which variables to read in
myvars <- c("year", "zip", "sex", "race", "age", "dual", "entry_age_break", "statecode",
            "followup_year", "followup_year_plus_one", "dead", "pm25_ensemble",
            "mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome", "medianhousevalue",
            "poverty", "education", "popdensity", "pct_owner_occ",
            "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax")
covar_vars <- c( myvars[12:26], 'region')



# define regions
NORTHEAST <- c("NY", "MA", "PA", "RI", "NH", "ME", "VT", "CT", "NJ")  
SOUTH <- c("DC", "VA", "NC", "WV", "KY", "SC", "GA", "FL", "AL", "TN", "MS", "AR", "MD", "DE", "OK", "TX", "LA")
MIDWEST <- c("OH", "IN", "MI", "IA", "MO", "WI", "MN", "SD", "ND", "IL", "KS", "NE")
WEST <- c("MT", "CO", "WY", "ID", "UT", "NV", "CA", "OR", "WA", "AZ", "NM")

# apply read and aggregate functions to all years
national_merged_list <- 
  lapply( files.mort,
          function(f){
            print( f)
            
            # read in the file
            file.in <- read.fst( f, columns = myvars, as.data.table = TRUE)
            
            # create zip codes as 5- character
            file.in[, zip := sprintf("%05d", zip)]
            
            # define region
            file.in[, region := 
                      ifelse( statecode %in% NORTHEAST, "NORTHEAST",
                              ifelse( statecode %in% SOUTH, "SOUTH",
                                      ifelse( statecode %in% MIDWEST, "MIDWEST",
                                              ifelse( statecode %in% WEST, "WEST",
                                                      NA))))]
            
            # take only complete cases
            colnames_check <- names( file.in)[ !names( file.in) == 'pm25_ensemble']
            national_merged <- file.in[ complete.cases( file.in[,..colnames_check]), ]
            
            # grab the year
            year.in <- unique( national_merged$year)
            
            # covariates
            covariates <- 
              national_merged[, lapply( .SD, min),
                              .SDcols =  covar_vars,
                              by = .( zip, year, statecode)]
            covariates[, `:=` ( year_fac = as.factor( year),
                                region = as.factor( region))]
            
            # aggregate_data
            # Generate count data for each individual characteristics and follow-up year
            national_merged[, time_count := followup_year_plus_one - followup_year]
            dead_personyear <- 
              national_merged[, lapply( .SD, sum),
                              .SDcols =  c( 'dead', 'time_count'),
                              by = .( zip, year, sex, race, dual, entry_age_break, followup_year)]
            confounders <- 
              national_merged[, lapply( .SD, min),
                              .SDcols =  covar_vars,
                              by = .( zip, year, sex, race, dual, entry_age_break, followup_year)]
            
            # merge the datasets
            aggregate_data <- merge(dead_personyear,
                                    confounders,
                                    by = c("zip", "year", "sex", "race", "dual", "entry_age_break", "followup_year"))
            
            # take summary statistics of key individual characteristics
            key_chars_sex <- 
              national_merged[, .N, by = 'sex']
            key_chars_race <- 
              national_merged[, .N, by = 'race']
            key_chars_dual <- 
              national_merged[, .N, by = 'dual']
            key_chars_age <- 
              national_merged[, .( age = mean( age), 
                                   age.sd = sd( age),
                                   age.min = min( age),
                                   age.max = max( age),
                                   N = .N)]
            
            # set names
            setnames( key_chars_sex, 'sex', 'char')
            setnames( key_chars_race, 'race', 'char')
            setnames( key_chars_dual, 'dual', 'char')
            
            # re-organize
            key_chars_sex[, char := paste0( 'sex.', char)]
            key_chars_race[, char := paste0( 'race.', char)]
            key_chars_dual[, char := paste0( 'dual.', char)]
            key_chars_sex.m <- dcast( key_chars_sex,
                                      . ~ char, value.var = 'N')[, `.` := NULL]
            key_chars_race.m <- dcast( key_chars_race,
                                       . ~ char, value.var = 'N')[, `.` := NULL]
            key_chars_dual.m <- dcast( key_chars_dual,
                                       . ~ char, value.var = 'N')[, `.` := NULL]
            key_chars_all <- 
              data.table( key_chars_sex.m, key_chars_race.m, 
                          key_chars_dual.m, key_chars_age)
            key_chars_all[, year := year.in]
            
            # create out list of years
            print( dim( national_merged))
            out <- 
              list( covariates = covariates,
                    aggregate_data = aggregate_data,
                    key_chars = key_chars_all)
            return( out)
          }) 

# extract all the covariates
covariates <- 
  lapply( national_merged_list,
          '[[', 'covariates') %>% 
  rbindlist

# extract all the aggregate data
aggregate_data <- 
  lapply( national_merged_list,
          '[[', 'aggregate_data') %>% 
  rbindlist

dim(aggregate_data)
# [1] 47878752       25



## ==================================================== ##
# read in new confounders that include ozone, no2, etc
## ==================================================== ##
confounders.f <- 
  list.files( file.path(dir_data, "cache_data", "confounders"), full.names = TRUE)

confounders_new <- 
  lapply( confounders.f, fread, drop = 'V1') %>% rbindlist
confounders_new[, zip := formatC( ZIP, width = 5, flag = '0')]

confounders_new[zip == '99363']
covariates[zip == '99363']

# merge with aggregate data
aggregate_data.additional_species <- 
  merge( aggregate_data,
         confounders_new[, .( zip, year, ozone.current_year,
                              no2.current_year, ozone_summer.current_year)],
         by = c( 'zip', 'year'), all.x = TRUE)

# merge with covariates
covariate.additional_species <- 
  merge( covariates,
         confounders_new[, .( zip, year, ozone.current_year,
                              no2.current_year, ozone_summer.current_year)],
         by = c( 'zip', 'year'), all.x = TRUE)



dim( aggregate_data.additional_species)
dim( covariate.additional_species)
dim( aggregate_data)
dim( covariates)

## ==================================================== ##
# save the data
## ==================================================== ##

write.fst(covariate.additional_species, file.path(dir_data, "cache_data", "covariates.fst"))
write.fst(aggregate_data.additional_species, file.path(dir_data, "cache_data", "aggregate_data.fst"))

covariates <- 
  read.fst( file.path(dir_data, "cache_data", "covariates.fst"), as.data.table = TRUE)
aggregate_data <- 
  read.fst( file.path(dir_data, "cache_data", "aggregate_data.fst"), as.data.table = TRUE)

## ==================================================== ##
# summarise key characteristics
## ==================================================== ##
# summarise key characteristics
key_characteristics <- 
  lapply( national_merged_list,
          '[[', 'key_chars') %>% 
  rbindlist

# sum person-years
sum( key_characteristics$N)

# percent female
key_characteristics[, sum( sex.2) / sum( N)]

# percent white
key_characteristics[, sum( race.1) / sum( N)]

# percent black
key_characteristics[, sum( race.2) / sum( N)]

# age
key_characteristics[, sum( age * N / sum( N))]

# dual
key_characteristics[, sum( dual.1) / sum( N)]

## ==================================================== ##
##  Read in the annual fst unit pm25 impacts
## ==================================================== ##
# set directory structure
# disperseR.base <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/'
# disperseR::create_dirs( disperseR.base)

# define pm25 exposure directory
# exp25_dir <- paste0( exp_dir, '25_new')

# exp25_dir2 <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/main/output/zips_model.lm.cv_single_poly'
exp25_dir2 <- '/n/dominici_nsaph_l3/Lab/projects/analytic/coal_exposure_pm25/zips_model.lm.cv_single_poly'

# files for raw, un-corrected hyads
exp_dir.raw <- 'data/data/cache_data/hyads_raw_exp'

zips.files.tot.yr <- list.files( exp25_dir2,
                                 pattern = 'zips_.*total_\\d{4}\\.fst',
                                 full.names = TRUE)
zips.files.tot_raw.yr <- list.files( exp_dir.raw,
                                     pattern = 'zips_.*total_\\d{4}\\.fst',
                                     full.names = TRUE)

# name each file by its year
names( zips.files.tot.yr) <- 1999:2020
names( zips.files.tot_raw.yr) <- 1999:2020


hyads_zips_tot <- lapply( 1999:2020,
                          function( y){
                            file.in <- zips.files.tot.yr[ as( y, 'character')]
                            in.dt <- read.fst( file.in, as.data.table = TRUE)
                            in.dt[, year := y]
                          }) %>% rbindlist()
hyads_zips_tot_raw <- lapply( 1999:2020,
                              function( y){
                                file.in <- zips.files.tot_raw.yr[ as( y, 'character')]
                                in.dt <- read.fst( file.in, as.data.table = TRUE)
                                in.dt[, year := y]
                              }) %>% rbindlist()

# scale raw hyads by its mean for scale issues
hyads_zips_tot_raw[, hyads := hyads / mean( hyads)]

# naming convention
setnames( hyads_zips_tot, c( 'vals.out', 'ZIP'),
          c( 'Y1', 'zip'))
setnames( hyads_zips_tot_raw, c( 'hyads', 'ZIP'),
          c( 'Y1_raw', 'zip'))

# merge the two
hyads_zips_tot_merge <- 
  merge( hyads_zips_tot,
         hyads_zips_tot_raw[,.( zip, year, Y1_raw)],
         by = c( 'zip', 'year'))

fwrite( hyads_zips_tot_merge,
        file.path(dir_data, "platform_data", 'pm25_annual_zip.csv'))
fwrite( hyads_zips_tot_merge,
        file.path("~/Desktop", 'pm25_annual_zip.csv'))


## ==================================================== ##
## scale coal pm by observation scaling factors by region-year
#   train the scaling factors
## ==================================================== ##
species.data_dir <- 'data/data/cache_data/aqs_annual'
years <- 1999:2020

## function: read AQS files
AQS_dataset <- function( pollutant.codes = c( 88307,  88321, 84316),
                         .species.data_dir = species.data_dir){
  aqs.files <- list.files( .species.data_dir,
                           pattern = 'annual_conc_by_monitor_.*\\.csv', full.names = TRUE)
  aqs.in <- rbindlist( lapply( aqs.files, fread))
  
  # read in monitor information
  aqs.mons <- fread( 'data/data/cache_data/aqs_annual/aqs_sites.csv')
  
  # merge site info and obs
  merge.names <- names( aqs.in)[which( names( aqs.in) %in% names( aqs.mons))]
  aqs.in.info <- merge( aqs.in, aqs.mons, by = merge.names, all.x = TRUE)
  
  # subset to input parameter.codes
  aqs.p <- aqs.in.info[grepl( paste0( pollutant.codes, collapse = '|'), `Parameter Code`),]
  
  # convert to spatial object
  crs.wgs84 <-  sf::st_crs( "+proj=longlat +datum=WGS84 +no_defs")
  crs.nad83 <-  sf::st_crs( "+proj=longlat +datum=NAD83 +no_defs")
  crs.use   <-  sf::st_crs( "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m")
  aqs.wgs84.sf <- sf::st_as_sf( aqs.p[ Datum == 'WGS84'], coords = c( 'Longitude', 'Latitude'), crs = crs.wgs84)
  aqs.nad83.sf <- sf::st_as_sf( aqs.p[ Datum == 'NAD83'], coords = c( 'Longitude', 'Latitude'), crs = crs.nad83)
  aqs.wgs84.sf <- sf::st_transform( aqs.wgs84.sf, crs.use)
  aqs.nad83.sf <- sf::st_transform( aqs.nad83.sf, crs.use)
  aqs.sf <- data.table( rbind( aqs.wgs84.sf, aqs.nad83.sf))
  
  # unique measurement ID
  aqs.sf$mID <- paste( aqs.sf$`State Code`, aqs.sf$`County Code`, 
                       aqs.sf$`Site Num`, aqs.sf$`Parameter Code`, sep = '.')
  aqs.sf$ID <- paste0( as.integer( aqs.sf$`State Code`), 
                       as.integer( aqs.sf$`County Code`), 
                       as.integer( aqs.sf$`Site Num`), 
                       as.integer( aqs.sf$`Parameter Code`))
  
  # AQS monitor location setting
  aqs.sf[ is.na( `Location Setting`) | `Location Setting` == '', `Location Setting` := 'MISSING']
  
  out <- sf::st_as_sf( aqs.sf)
  return( out)
}

## Function: evaluation metrics
eval.fn <- function( Yhat, Yact, mod.name){
  num.diff <- sum( Yhat - Yact, na.rm = T)
  abs.diff <- sum( abs( Yhat - Yact), na.rm = T)
  denom <- sum( Yact, na.rm = T)
  metrics <- data.table( year = mod.name,
                         N = length( Yhat),
                         NMB = num.diff / denom,
                         NME = abs.diff / denom,
                         MB   = num.diff / length( Yhat),
                         RMSE = sqrt( sum( ( Yhat - Yact) ^ 2, na.rm = T) / length( Yhat)),
                         R.p = cor( Yhat, Yact, use = 'complete.obs') ^ 2,
                         R.s = cor( Yhat, Yact, use = 'complete.obs', method = 'spearman') ^ 2)
  return( metrics)
}

# read in sulfate from improve monitors
AQS_sulfate <- AQS_dataset( pollutant.codes = c( 88403))

# give states a bin ##REMINDER--check Washington, DC!
state_bins <- data.table( state = state.name,
                          STATE = state.abb)
state_bins[ state %in% c( 'Mississippi', 'Georgia', 'Alabama', 
                          'Florida', 'South Carolina', 
                          'Tennessee', 'North Carolina'),
            statebin_facility := 'Southeast']
state_bins[ state %in% c( 'West Virginia', 'Delaware', 'Kentucky',
                          'Maryland', 'Virginia', 'New Jersey', 
                          'Pennsylvania'),
            statebin_facility := 'Mid-Atlantic']
state_bins[ state %in% c( 'Arkansas', 'Louisiana', 'Oklahoma', 'Texas'),
            statebin_facility := 'Southcentral']
state_bins[ state %in% c( 'New York', 'New Hampshire', 'Massachusetts', 'Maine',
                          'Connecticut', 'Vermont', 'Rhode Island'),
            statebin_facility := 'Northeast']
state_bins[ state %in% c( 'Nevada', 'Utah', 'Arizona', 'New Mexico',
                          'Colorado', 'Wyoming', 'California', 'Idaho',
                          'Washington', 'Oregon', 'Montana'),
            statebin_facility := 'West']
state_bins[ state %in% c( 'Wisconsin','Michigan', 'Illinois',
                          'Indiana', 'Ohio'),
            statebin_facility := 'Midwest east']
state_bins[ state %in% c( 'North Dakota', 'South Dakota', 
                          'Minnesota', 'Kansas', 'Iowa', 
                          'Nebraska', 'Missouri'),
            statebin_facility := 'Midwest west']

# merge with observationdataset
AQS_sulfate_reg <- 
  merge( AQS_sulfate, state_bins,
         by.x = 'State Name', by.y = 'state')
AQS_sulfate_reg <- AQS_sulfate_reg[grep( 'IMPROVE', AQS_sulfate_reg$`Method Name`),]

# gridded hyads file location
hyads_file_loc <- '/n/dominici_nsaph_l3/Lab/projects/analytic/coal_exposure_pm25/grids_model.lm.cv_single_poly'

#coordinate reference system projection string for spatial data
p4s <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m'

# get the names of the gridded HyADS output files
grid.files.yr <- list.files( hyads_file_loc,
                             pattern = 'grids_pm25_total_\\d{4}\\.fst',
                             full.names = TRUE)

# read select files
grid.dat <- lapply( grid.files.yr,
                    function( f){
                      year.f <- gsub( '^.*_|\\.fst', '', f)
                      
                      in.f <- read.fst( f, as.data.table = T)
                      in.f[, year := year.f]
                    }) %>% rbindlist

# calculate average by year
grid.dat[, .( mean = mean( vals.out, na.rm = TRUE),
              X_25 = quantile( vals.out, .25, na.rm = TRUE),
              X_75 = quantile( vals.out, .75, na.rm = TRUE)), by = year]

# dcast to get year columns, convert to sfåˆ
grid.dat.c <- dcast( grid.dat, x + y ~ year, value.var = 'vals.out')
grid.dat.sf <- 
  rasterFromXYZ( grid.dat.c, crs = p4s) %>%
  rasterToPolygons() %>%
  st_as_sf() 

## link the IMPROVE & HyADS datasets
years <- 1999:2020
eval_out <- 
  lapply( years,
          function( y) {
            print( y)
            # pull coal pm for this year
            hyads.yr.name <- paste0( 'X', y)
            hyads.yr <- grid.dat.sf[, hyads.yr.name]
            
            # pulll the observations
            obs.yr <- AQS_sulfate_reg[ AQS_sulfate_reg$Year == y,
                                       c( 'Arithmetic Mean', 'statebin_facility')]
            
            # join obs in available grid cells
            hyads_obs_join <- 
              st_join( st_transform( obs.yr, p4s), 
                       hyads.yr, 
                       join = st_within) %>%
              as.data.table() %>% 
              na.omit() %>%
              setnames( c( 'Arithmetic Mean', hyads.yr.name),
                        c( 'obs', 'hyads'))
            
            # add year label
            hyads_obs_join[, year := y]
            # hyads_obs_join <- hyads_obs_join[hyads > quantile(hyads, .2, na.rm = TRUE)]
            
            # do the evaluation
            hyads_obs_eval <- 
              hyads_obs_join[, eval.fn( hyads, obs, y), 
                             by = statebin_facility]
            
            return( list( hyads_obs_eval = hyads_obs_eval,
                          hyads_obs_join = hyads_obs_join))
            
          }
  ) 

eval_out.m <- 
  lapply( eval_out, '[[', 'hyads_obs_eval') %>%
  rbindlist( ) %>%
  melt( id.vars = c( 'year', 'statebin_facility'))

# create table of hyads & improve by region
# plot average change over time by region
scatter.dt <- 
  lapply( eval_out, '[[', 'hyads_obs_join') %>%
  rbindlist( )
scatter.dt_reg <- scatter.dt[, .( hyads = mean( hyads),
                                  obs = mean( obs)),
                             by = .( statebin_facility, year)] %>%
  melt( id.vars = c( 'statebin_facility', 'year'))

# relative to 1999
scatter.dt_reg <- scatter.dt_reg[ year == 1999][, year := NULL] %>%
  setnames( 'value', 'value99') %>%
  merge( scatter.dt_reg, by = c( 'statebin_facility', 'variable'))
scatter.dt_reg[ , value.rel := value / value99]

# melt to create single plot
scatter.dt_reg.m <- 
  melt( scatter.dt_reg,
        id.vars = c( 'statebin_facility', 'variable', 'year'),
        measure.vars = c( 'value', 'value.rel'),
        variable.name = 'model')

# change the names
scatter.dt_reg.m[ model == 'value', 
                  model.name := "Average~'['*µg~m^{-3}*']'"]
scatter.dt_reg.m[ model == 'value.rel', 
                  model.name := "Average~relative~to~1999"]
scatter.dt_reg.m[, statebin_facility := gsub( ' ', '~', statebin_facility)]

## Figure S8: IMPROVE sulfate vs. HyADS
gg_region_adj <- 
  ggplot( scatter.dt_reg.m,
          aes( x = year, y = value, color = variable, group = variable)) + 
  geom_hline( yintercept = 0) +
  geom_line( size = 2) +
  scale_color_brewer( palette = 'Dark2',
                      labels = c( 'hyads' = expression( Coal~PM[2.5]),
                                  'obs' = expression( Observed~SO[4]~PM[2.5]))) +
  # labs( y = expression(Coal[SO2]~PM[2.5]~and~observed~SO[4]~PM[2.5]~relative~to~1999)) +
  facet_grid( model.name ~ statebin_facility, 
              scales = 'free_y',
              switch = 'y',
              space = 'free_y',
              labeller = label_parsed) + 
  theme_minimal() + 
  theme( axis.text = element_text( size = 14),
         axis.text.x = element_text( angle = 30),
         axis.title = element_blank(),
         legend.background = element_rect( fill = 'white', color = 'white'),
         legend.position = c( .9, .9),
         legend.text = element_text( size = 14),
         legend.text.align = 0,
         legend.title = element_blank(),
         plot.background = element_rect( color = NULL),
         plot.margin = margin(.5,.5,.5,.5, "cm"),
         strip.placement = "outside",
         strip.text = element_text( size = 18),
         strip.text.x = element_text( size = 18))
# gg_region_adj

ggsave( 'figures/hyads_sulfate_adjustments_rel.png',
        gg_region_adj, height = 6, width = 12, units = 'in', scale = 1.25)

# save the adjustment factors
scatter.dt_reg <- dcast( scatter.dt_reg,
                         statebin_facility + year ~ variable,
                         value.var = 'value.rel')
scatter.dt_reg[, obs_adjust := obs / hyads]

# save the output
write.csv( scatter.dt_reg,
           file.path( dir_data, 'coalpm25_to_obs_adjustment_factors.csv'))



## ==================================================== ##
## do the scaling scale coal pm by observation scaling factors by region-year
## ==================================================== ##
zip_state <- covariates[, .( zip, statecode)] %>% unique()

# read the adjustment factors
hyads_scale_factors <- 
  fread( file.path( dir_data, 'coalpm25_to_obs_adjustment_factors.csv'), drop = 'V1')

# define state bins
state_bins <- data.table( state = c( state.name, 'Washington, DC'),
                          STATE = c( state.abb, 'DC')) %>%
  merge( zip_state, by.x = 'STATE', by.y = 'statecode')

state_bins[ state %in% c( 'Mississippi', 'Georgia', 'Alabama', 
                          'Florida', 'South Carolina', 
                          'Tennessee', 'North Carolina'),
            statebin_facility := 'Southeast']
state_bins[ state %in% c( 'West Virginia', 'Delaware', 'Kentucky',
                          'Maryland', 'Virginia', 'New Jersey', 
                          'Pennsylvania', 'Washington, DC'),
            statebin_facility := 'Mid-Atlantic']
state_bins[ state %in% c( 'Arkansas', 'Louisiana', 'Oklahoma', 'Texas'),
            statebin_facility := 'Southcentral']
state_bins[ state %in% c( 'New York', 'New Hampshire', 'Massachusetts', 'Maine',
                          'Connecticut', 'Vermont', 'Rhode Island'),
            statebin_facility := 'Northeast']
state_bins[ state %in% c( 'Nevada', 'Utah', 'Arizona', 'New Mexico',
                          'Colorado', 'Wyoming', 'California', 'Idaho',
                          'Washington', 'Oregon', 'Montana'),
            statebin_facility := 'West']
state_bins[ state %in% c( 'Wisconsin','Michigan', 'Illinois',
                          'Indiana', 'Ohio'),
            statebin_facility := 'Midwest east']
state_bins[ state %in% c( 'North Dakota', 'South Dakota', 
                          'Minnesota', 'Kansas', 'Iowa', 
                          'Nebraska', 'Missouri'),
            statebin_facility := 'Midwest west']

# merge hyads zips with states
hyads_zips_tot_state <- 
  merge( hyads_zips_tot_merge, state_bins,
         by = 'zip', all = TRUE) %>%
  merge( hyads_scale_factors,
         by = c( 'year', 'statebin_facility'),
         all = TRUE)

hyads_zips_tot_state[, Y1.adj := Y1 * obs_adjust]
hyads_zips_tot_state <- hyads_zips_tot_state[!is.na( year)]
## ==================================================== ##
## save the data
## ==================================================== ##
write.fst( hyads_zips_tot_state, file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'))
hyads_zips_tot_state <- 
  read.fst(  file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
             as.data.table = TRUE)

## ==================================================== ##
## data summaries
## ==================================================== ##
summary( hyads_zips_tot_state[year == 1999])
summary( hyads_zips_tot_state[year == 2020])
summary( hyads_zips_tot_state[, Y1.adj])
summary( hyads_zips_tot_state[, Y1])
mean( hyads_zips_tot_state[, Y1.adj], na.rm = T)
mean( hyads_zips_tot_state[, Y1], na.rm = T)
sd( hyads_zips_tot_state[, Y1], na.rm = T)

## ==================================================== ##
## check correlation of each variable with coal PM
## ==================================================== ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj'), as.data.table = TRUE)
covariates_use <- merge(covariates, dat_annual, by = c("zip", "year")) #, all.x = TRUE)
aggregate_data_use <- merge(aggregate_data, dat_annual, by = c("zip", "year")) #, all.x = TRUE)

# calculate portion of PM that is not coal PM
covariates_use[, pm_not_coalPM := pm25_ensemble - Y1]
aggregate_data_use[, pm_not_coalPM := pm25_ensemble - Y1]

# calculate correlations
cor( covariates_use[, .( Y1, pm_not_coalPM, pm25_ensemble, mean_bmi, smoke_rate,
                         hispanic, pct_blk, medhouseholdincome,
                         medianhousevalue, poverty, education, popdensity,
                         pct_owner_occ, summer_tmmx, winter_tmmx, summer_rmax,
                         winter_rmax, ozone.current_year,
                         no2.current_year, ozone_summer.current_year)],
     use = 'pairwise.complete.obs')


aggregate_data_use[, lapply( .SD, function( x){ list( mean = sum( x * time_count, na.rm = TRUE) / sum( time_count),
                                                      sd = sd( x, na.rm = TRUE),
                                                      min = min( x, na.rm = TRUE),
                                                      max = max( x, na.rm = TRUE))}),
                   .SDcols = c( 'Y1', 'pm_not_coalPM', 'pm25_ensemble','ozone.current_year',
                                'no2.current_year', 'ozone_summer.current_year')]















