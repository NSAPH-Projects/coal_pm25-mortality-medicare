library( fst)
library( data.table)
library( gnm)
library( magrittr)
library( pryr)
library( sf)
library( pbmcapply)

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/' 
dir_out <- 'results/'

aggregate_data <- read.fst( paste0(dir_data, "aggregate_data.fst"),
                            as.data.table = TRUE)

# Caculate the number of death and person-time in each zip code and year
dead_personyear <- 
  aggregate_data[, lapply( .SD, sum),
                 .SDcols =  c( 'dead', 'time_count'),
                 by = .( zip, year)]

# retrieve ZIP code PM concs
aggregate_pm25 <- 
  aggregate_data[, .( zip, year, pm25_ensemble)] %>% unique

rm( aggregate_data)

## ====================================================== ##
# project deaths into future years - assign average from 2014-2016
## ====================================================== ##
# create datasets needed to fill
zips_year <- 
  expand.grid( 
    zip = unique( dead_personyear[ year %in% 2014:2016]$zip),
    year = 2017:2020)

#calculate average over 3 years for baseline
dat_year_avgs <- 
  dead_personyear[ year %in% 2014:2016][, .( deaths.fill = mean( dead, na.rm = T),
                                             denom.fill = mean( time_count, na.rm = T)),
                                        by = .( zip)]

# merge all the fill data
dat_year_fill <- merge( dead_personyear, zips_year, by = c( 'zip', 'year'), all = T) %>%
  merge( dat_year_avgs, by = c( 'zip'), all = T)%>%
  merge( aggregate_pm25, by = c( 'zip', 'year'), all = T)

# fill in the missing months with deaths
dat_year_fill[ year <= 2016 | (year > 2016 & !is.na( dead)), 
               `:=` ( deaths.fill = dead,
                      denom.fill = time_count)]

# merge with pm25 ensemble data

# dat_year_fill <- dat_year_fill[zip %in% zips_all]

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( 'data/cache_data/hyads_pm25_annual.fst',
                        columns = c('zip','year', 'Y1'), 
                        as.data.table = TRUE)

# merge with medicare/confounders
dat_year_fill <- merge(dat_year_fill, dat_annual, 
                       by = c("zip", "year")) #, all.x = TRUE)


## ====================================================== ##
# load poisson model results
## ====================================================== ##
num_uniq_zip <- length( unique( dat_year_fill$zip))

# read in the coefficients
poisson_coefs <- fread( paste0(dir_out, "poisson_model_coefs.csv"))

# read the confidence interval dataset
# loads the object bootstrap_CI
load( paste0(dir_out, "loglinear_coefs_boots.RData"))

# calculate 95% CI
num_uniq_zip <- length( unique( dat_year_fill$zip))
boot_beta_pm <- 
  1.96 * sd( bootstrap_CI$pm25_ensemble) * 
  sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)
boot_beta_hy <- 
  1.96 * sd( bootstrap_CI$hyads) * 
  sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)
boot_beta_hy.adj <- 
  1.96 * sd( bootstrap_CI$hyads.adj) * 
  sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)

# create list of betas
betas_pm <- 
  poisson_coefs[ model == 'pm25_ensemble']$Estimate +
  c( -boot_beta_pm, 0, +boot_beta_pm)
betas_hy <- 
  poisson_coefs[ model == 'hyads']$Estimate +
  c( -boot_beta_hy, 0, +boot_beta_hy)
betas_hy.adj <- 
  poisson_coefs[ model == 'hyads.adj']$Estimate +
  c( -boot_beta_hy.adj, 0, +boot_beta_hy.adj)
betas_kr <- log( c( 1.035, 1.056, 1.078)) / 10

exp(1 * betas_pm)
exp(1 * betas_hy)
exp(1 * betas_hy.adj)
exp(1 * betas_kr)

## ==================================================== ##
##  Read in the ZIP code spatial inputs (State variable)
## ==================================================== ##
# get the common functions
source( 'code/utilities/common_functions.R')

direct.dat <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/main/input/zcta_500k'
# direct.dat <- '~/Dropbox/Harvard/Manuscripts/Energy_Transitions/Data_and_code/data/gis'

zips <- zip_sf_reader( direct.dat, 'RCE') %>%
  data.table()
setnames( zips, 'ZIP', 'zip')

# convert zip, month, year, and state to factors
dat_year_fill <- merge( dat_year_fill, 
                        zips[, .( zip, STATE)], by = 'zip', all.x = T)

write.fst( dat_year_fill,
           'data/dat_year_fill.fst')
## ====================================================== ##
# sum deaths by state & year as input for adjoint comparison
## ====================================================== ##
sum_deaths_year <- dat_year_fill[, .( deaths.fill = sum( deaths.fill, na.rm = T),
                                      denom.fill = sum( denom.fill, na.rm = T)),
                                 by = .( year, STATE)]
write.fst( sum_deaths_year,
           'data/total_deaths_by_state_year.fst')
#
## ====================================================== ##
# sum all deaths
## ====================================================== ##
sum( dat_year_fill$deaths.fill)
## ====================================================== ##
# calculate deaths avoided for all coal pm2.5
## ====================================================== ##
deaths_by_year_all.hy <- 
  lapply( betas_hy,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_hy == b))]
          }) %>% rbindlist

deaths_by_year_all.pm <- 
  lapply( betas_pm,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_pm == b))]
          }) %>% rbindlist

deaths_by_year_all.kr <- 
  lapply( betas_kr,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_kr == b))]
          }) %>% rbindlist

# add model labels
deaths_by_year_all.hy[, model := 'hyads']
deaths_by_year_all.pm[, model := 'pm25_ensemble']
deaths_by_year_all.kr[, model := 'pm25_krewski']

# combine into single data.table
deaths_by_year_all <- 
  rbind( deaths_by_year_all.hy,
         deaths_by_year_all.pm,
         deaths_by_year_all.kr)

# sum by state
deaths_by_year_state_all <- 
  deaths_by_year_all[, .( 
    deaths_coef   = sum( deaths_coef, na.rm = T)),
    by = .( year, STATE, model, beta)]

deaths_by_year_all.c <- 
  dcast( deaths_by_year_state_all,
         year + STATE + model ~ beta,
         value.var = 'deaths_coef')
setnames( deaths_by_year_all.c, unique( deaths_by_year_state_all$beta),
          paste0( 'deaths_coef_', unique( deaths_by_year_state_all$beta)))

# save it
write.fst( deaths_by_year_all.c, paste0(dir_out, "deaths_by_year_totallHyADS.fst"))

## ====================================================== ##
# calculate deaths avoided for total pm2.5
## ====================================================== ##
deaths_by_year_all.pm25_ensemble <- 
  lapply( betas_pm,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * pm25_ensemble) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_pm == b))]
          }) %>% rbindlist

# sum by state
deaths_by_year_state_all.pm25_ensemble <- 
  deaths_by_year_all.pm25_ensemble[, .( 
    deaths_coef   = sum( deaths_coef, na.rm = T)),
    by = .( year, STATE, beta)]

deaths_by_year_state_all.pm25_ensemble.c <- 
  dcast( deaths_by_year_state_all.pm25_ensemble,
         year + STATE ~ beta,
         value.var = 'deaths_coef')
setnames( deaths_by_year_state_all.pm25_ensemble.c, 
          unique( deaths_by_year_state_all.pm25_ensemble$beta),
          paste0( 'deaths_coef_', unique( deaths_by_year_state_all.pm25_ensemble$beta)))

# save it
write.fst( deaths_by_year_state_all.pm25_ensemble.c, paste0(dir_out, "deaths_by_year_totall_pm25_ensemble.fst"))




