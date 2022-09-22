#devtools::install_github( 'lhenneman/disperseR@dev')
library( disperseR)
library( fst)
library( data.table)
library( parallel)
library( xgboost)
library( dplyr)
library( foreign)

dir_data <- 'data/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

## ==================================================== ##
##  Read in covariates data
## ==================================================== ##
files.mort <- list.files("/nfs/nsaph_ci3/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2",
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




# save(covariates, file = paste0(dir_data, "covariates.RData"))
# save(aggregate_data, file = paste0(dir_data, "aggregate_data.RData"))

write.fst(covariates, file.path(dir_data, "cache_data", "covariates.fst"))
write.fst(aggregate_data, file.path(dir_data, "cache_data", "aggregate_data.fst"))

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
disperseR.base <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/'
disperseR::create_dirs( disperseR.base)

# define pm25 exposure directory
exp25_dir <- paste0( exp_dir, '25_new')

exp25_dir2 <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/main/output/zips_model.lm.cv_single_poly'
zips.files.tot.yr <- list.files( exp25_dir2,
                                 pattern = 'zips_.*total_\\d{4}\\.fst',
                                 full.names = TRUE)
names( zips.files.tot.yr) <- 1999:2020
hyads_zips_tot <- lapply( 1999:2020,
                          function( y){
                            file.in <- zips.files.tot.yr[ as( y, 'character')]
                            in.dt <- read.fst( file.in, as.data.table = TRUE)
                            in.dt[, year := y]
                          }) %>% rbindlist()
setnames( hyads_zips_tot, c( 'vals.out', 'ZIP'),
          c( 'Y1', 'zip'))

fwrite( hyads_zips_tot,
        file.path(dir_data, "platform_data", 'pm25_annual_zip.csv'))
fwrite( hyads_zips_tot,
        file.path("~/Desktop", 'pm25_annual_zip.csv'))

## ==================================================== ##
## scale coal pm by observation scaling factors by region-year
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
  merge( hyads_zips_tot, state_bins,
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
summary( hyads_zips_tot_state[year == 1999])
summary( hyads_zips_tot_state[, Y1.adj])
summary( hyads_zips_tot_state[, Y1])
mean( hyads_zips_tot_state[, Y1.adj], na.rm = T)
mean( hyads_zips_tot_state[, Y1], na.rm = T)
sd( hyads_zips_tot_state[, Y1], na.rm = T)
## ==================================================== ##
## population-weighted hyads
## ==================================================== ##





