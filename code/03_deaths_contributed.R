#renv::install("lhenneman/disperseR@dev")

library( fst)
library( data.table)
library( gnm)
library( magrittr)
library( pryr)
library( sf)
library( pbapply)

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/data' 
dir_out <- 'results/'

aggregate_data <- 
  read.fst( 
    file.path(dir_data, "cache_data", "aggregate_data.fst"),
    as.data.table = TRUE)

# Caculate the number of death and person-time in each zip code and year
dead_personyear <- 
  aggregate_data[, lapply( .SD, sum),
                 .SDcols =  c( 'dead', 'time_count'),
                 by = .( zip, year)]

# retrieve ZIP code PM concs
aggregate_pm25 <- 
  aggregate_data[, .( zip, year, pm25_ensemble)] %>% unique

# match data used to fit models
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), as.data.table = TRUE)

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
dat_annual <- read_fst( 'data/data/cache_data/hyads_pm25_annual.fst',
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), 
                        as.data.table = TRUE)

# merge with medicare/confounders
dat_year_fill <- merge(dat_year_fill, dat_annual, 
                       by = c("zip", "year")) #, all.x = TRUE)
summary( data.table( dat_year_fill)[year == 1999])
summary( data.table( dat_year_fill)[year == 2020])

# calculate population-weighted hyads
data.table( dat_year_fill)[, .( Y1.pw = sum( Y1 * denom.fill) / 
                                  sum( denom.fill)),
                           by = year]

## ====================================================== ##
# load poisson model results
## ====================================================== ##
num_uniq_zip <- length( unique( dat_year_fill$zip))

# read in the coefficients
poisson_coefs <- fread( paste0(dir_out, "poisson_model_coefs.csv"), drop = 'V1')

# read the confidence interval dataset
# loads the object bootstrap_CI
load( paste0(dir_out, "loglinear_coefs_boots.RData"))

# calculate 95% CI
num_uniq_zip <- length( unique( dat_year_fill$zip))

beta_ci <- 
  lapply( names( bootstrap_CI)[ names( bootstrap_CI) != 'boots_id'],
          function( x){
            print( x)
            
            # calculate CI width
            ci <- 1.96 * sd( unlist( bootstrap_CI[,..x])) * 
              sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)
            
            # three points: mean estimate, ± CI
            out <- poisson_coefs[ model == x]$Estimate +
              c( -ci, 0, +ci) %>% data.table()
            
            # add labels
            out[, `:=` ( model = ..x,
                         ci = c( 'low', 'mid', 'high'))]
            setnames( out, '.', 'value')
            return( out)
          }) %>% rbindlist()

# add on the Krewski et al estimates
betas_kr <- data.table( model = 'krewski', 
                        low = log( 1.035) / 10, 
                        mid = log( 1.056) / 10, 
                        high = log( 1.078) / 10)

beta_ci.dt <- 
  dcast( beta_ci, model ~ ci, value.var = 'value') %>%
  rbind( betas_kr)

## ==================================================== ##
##  Check out some of the betas
## ==================================================== ##
# hyads, adj. by rPM; by NO2; by rPM & NO2
beta_ci.dt[ , lapply( .SD, function( x) round( exp(x), 4)),
                 .SDcols = c( 'mid', 'low', 'high'),
                 by = model]

# check RR for hyads & raw hyads per standard deviation
sd_hyads <- sd( dat_annual$Y1)
sd_hyads_raw <- sd( dat_annual$Y1_raw)

exp( beta_ci.dt[ model %in% c( 'hyads'), .( low, mid, high)] * sd_hyads)
exp( beta_ci.dt[ model %in% c( 'hyads_raw'), .( low, mid, high)] * sd_hyads_raw)


## ==================================================== ##
##  Read in the ZIP code spatial inputs (State variable)
## ==================================================== ##
# create zip code shape file reader function
# read zcta shapefile and crosswalk
zip_sf_reader <- function( d = direct.dat){
  zcta_shapefile <- file.path( d, 'cb_2017_us_zcta510_500k.shp')
  
  # zcta-ZIP crosswalk file downloaded from 'http://mcdc2.missouri.edu/data/corrlst/'
  cw <- disperseR::crosswalk
  zips <- st_read(zcta_shapefile)
  setnames( zips, 'ZCTA5CE10', 'ZCTA')
  zips <- merge( zips, cw, by = "ZCTA", all = F, allow.cartesian = TRUE)
  
  return( zips)
}

direct.dat <- 'data/data/zcta_500k'

zips <- zip_sf_reader( direct.dat) %>%
  data.table()
setnames( zips, 'ZIP', 'zip')

# convert zip, month, year, and state to factors
dat_year_fill <- merge( dat_year_fill, 
                        zips[, .( zip, STATE)], by = 'zip', all.x = T)

write.fst( dat_year_fill,
           'data/data/cache_data/dat_year_fill.fst')
## ====================================================== ##
# sum deaths by state & year as input for adjoint comparison
## ====================================================== ##
sum_deaths_year <- dat_year_fill[, .( deaths.fill = sum( deaths.fill, na.rm = T),
                                      denom.fill = sum( denom.fill, na.rm = T)),
                                 by = .( year, STATE)]
write.fst( sum_deaths_year,
           'data/data/cache_data/total_deaths_by_state_year.fst')
#
## ====================================================== ##
# fix deaths/rates and hyads at 1999 levels
## ====================================================== ##
dat_1999 <- 
  dat_year_fill[year == 1999, .( zip, deaths.fill, denom.fill, Y1)]

setnames( dat_1999,
          c( 'deaths.fill', 'denom.fill', 'Y1'),
          paste0( c( 'deaths.fill', 'denom.fill', 'Y1'), '_99'))

dat_year_fill <- 
  merge( dat_year_fill,
         dat_1999,
         by = c( 'zip'),
         all = TRUE)

## ====================================================== ##
# sum all deaths
## ====================================================== ##
sum( dat_year_fill$deaths.fill)

## ====================================================== ##
# calculate deaths avoided for all coal pm2.5
## ====================================================== ##
death_calculator <- 
  function( model_n, high_mid_low = c( 'low', 'mid', 'high'), betas = beta_ci.dt,
            fulldata = copy( dat_year_fill),
            deaths_use = 'deaths.fill',
            denom_use = 'denom.fill',
            exposure_use = 'Y1',
            by.vars = c( 'zip', 'year', 'STATE')){
    print( model_n)
    
    # change names in data table to match selections
    setnames( fulldata,
              c( deaths_use, denom_use, exposure_use),
              c( 'deaths_use', 'denom_use', 'exposure_use'))
    
    # run the function for each beta
    deaths_out <- 
      lapply( high_mid_low,
              function( hml){
                # select correct beta
                beta_use <- betas[ model == model_n, ..hml] %>% unlist
                
                # calculate deaths
                deaths_by_zip.b <- 
                  fulldata[, .( 
                    deaths_coef = deaths_use / denom_use * ( exp( beta_use * exposure_use) - 1) * denom_use#,
                  ), by = by.vars]
                
                beta_names <- c( 'low' = '1', 'mid' = '2', 'high' = '3')
                deaths_by_zip.b[, `:=` (beta = beta_names[hml],
                                        model = model_n)]
                
                return( deaths_by_zip.b)
              }) %>% rbindlist
    return( deaths_out)
  }

# seperate betas that are hyads.adj or raw hyads
betas_hyads.adj <- c( 'hyads.adj')
betas_hyads_raw <- c( 'hyads_raw')

deaths_by_year_all <- 
  lapply( beta_ci.dt$model[!(beta_ci.dt$model %in% betas_hyads.adj | 
                               beta_ci.dt$model %in% betas_hyads_raw)],
          death_calculator) %>% rbindlist
deaths_by_year_all.adj <- 
  lapply( beta_ci.dt$model[beta_ci.dt$model %in% betas_hyads.adj],
          death_calculator,
          exposure_use = 'Y1.adj') %>% rbindlist
deaths_by_year_all_raw <- 
  lapply( beta_ci.dt$model[beta_ci.dt$model %in% betas_hyads_raw],
          death_calculator,
          exposure_use = 'Y1_raw') %>% rbindlist

deaths_by_year_all <- 
  list( deaths_by_year_all,
        deaths_by_year_all.adj,
        deaths_by_year_all_raw) %>%
  rbindlist( )

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

# check sums
c( 
  low = deaths_by_year_all.c[ model == 'hyads', sum( deaths_coef_1)],
  mid = deaths_by_year_all.c[ model == 'hyads', sum( deaths_coef_2)],
  high = deaths_by_year_all.c[ model == 'hyads', sum( deaths_coef_3)])
c( 
  low = deaths_by_year_all.c[ model == 'hyads_raw', sum( deaths_coef_1)],
  mid = deaths_by_year_all.c[ model == 'hyads_raw', sum( deaths_coef_2)],
  high = deaths_by_year_all.c[ model == 'hyads_raw', sum( deaths_coef_3)])

# save it
write.fst( deaths_by_year_all.c, paste0(dir_out, "deaths_by_year_total_coal_pm25.fst"))

## ====================================================== ##
# calculate deaths avoided for total pm2.5
## ====================================================== ##
deaths_by_year_all.pm25_ensemble <- 
  death_calculator( "pm25_ensemble",
                    exposure_use = 'pm25_ensemble') 

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
write.fst( deaths_by_year_state_all.pm25_ensemble.c, paste0(dir_out, "deaths_by_year_total_pm25_ensemble.fst"))

## ====================================================== ##
# calculate deaths avoided for hyads; constant death rate or constant hyads
## ====================================================== ##
deaths_all.hyads_constdeath <- 
  death_calculator( "hyads",
                    exposure_use = 'Y1',
                    deaths_use = 'deaths.fill_99',
                    denom_use = 'denom.fill_99') 
deaths_all.hyads_consthy <- 
  death_calculator( "hyads",
                    exposure_use = 'Y1_99',
                    deaths_use = 'deaths.fill',
                    denom_use = 'denom.fill') 

# sum by year
deaths_by_year_all.hyads_constdeath <- 
  deaths_all.hyads_constdeath[, .( 
    deaths_coef   = sum( deaths_coef, na.rm = T)),
    by = .( year, beta)]
deaths_by_year_all.hyads_consthy <- 
  deaths_all.hyads_consthy[, .( 
    deaths_coef   = sum( deaths_coef, na.rm = T)),
    by = .( year, beta)]

# add model label
deaths_by_year_all.hyads_constdeath[, case := 'Constant deaths & death rates']
deaths_by_year_all.hyads_consthy[, case := 'Constant coal PM2.5']

# rbind & dcast
deaths_by_year_sensitivities <- 
  rbind( deaths_by_year_all.hyads_constdeath,
         deaths_by_year_all.hyads_consthy)

deaths_by_year_sensitivities.c <- 
  dcast( deaths_by_year_sensitivities,
         year + case ~ beta,
         value.var = 'deaths_coef')
setnames( deaths_by_year_sensitivities.c, 
          unique( deaths_by_year_sensitivities$beta),
          paste0( 'deaths_coef_', unique( deaths_by_year_sensitivities$beta)))

# calculate 1999 values
deaths_by_year_sensitivities.c_99 <- 
  deaths_by_year_sensitivities.c[year == 1999][, year := NULL]
setnames( deaths_by_year_sensitivities.c_99,
          c( 'deaths_coef_1', 'deaths_coef_2', 'deaths_coef_3'),
          c( 'deaths_coef_1_99', 'deaths_coef_2_99', 'deaths_coef_3_99'))
deaths_by_year_sensitivities.c <- 
  merge( deaths_by_year_sensitivities.c,
         deaths_by_year_sensitivities.c_99,
         by = c( 'case'))
deaths_by_year_sensitivities.c[, `:=`( 
  deaths_coef_1_diff = (deaths_coef_1 - deaths_coef_1_99) / deaths_coef_1_99,
  deaths_coef_2_diff = (deaths_coef_2 - deaths_coef_2_99) / deaths_coef_2_99,
  deaths_coef_3_diff = (deaths_coef_3 - deaths_coef_3_99) / deaths_coef_3_99
)]

# save it
fwrite( deaths_by_year_sensitivities.c, 
        paste0(dir_out, "deaths_by_year_sensitivities.csv"))

## ====================================================== ##
# get filenames for unit-level coal PM2.5
## ====================================================== ##
# define pm25 exposure directory
exp_dir <- '/n/dominici_nsaph_l3/Lab/projects/analytic/coal_exposure_pm25/zips'
exp25_dir <- '/n/dominici_nsaph_l3/Lab/projects/analytic/coal_exposure_pm25/zips_model.lm.cv_single_poly'

# get files for unit-specific hyads
zip.units.yr <- list.files( exp25_dir,
                            pattern = 'zips_.*byunit.*\\d{4}\\.fst',
                            full.names = TRUE)
zip.units.yr_raw <- list.files( exp_dir,
                            pattern = 'zips_.*byunit.*\\d{4}\\.fst',
                            full.names = TRUE)

# put together data.table for easy reference
zip.units.yr.dt <- data.table( f = zip.units.yr,
                               y = 1999:2020)
zip.units.yr.dt_raw <- data.table( f = zip.units.yr_raw,
                               y = 1999:2020)
## ====================================================== ##
# do the risk assessment
## ====================================================== ##
# write the function
# cuts defines number of units run at a time
# cuts = 2 seems to work okay for 100 samples per zip
# cuts = 10 seems to work okay for 3 samples per zip
risk_assessmenter <- function( n, fyms, 
                               model_n, 
                               scale_by = NULL,
                               cuts.n = 2){
  print( n)
  # message( paste( 'Mem used is', mem_used()))
  fym <- fyms[n]
  f <- fym$f
  y <- fym$y
  
  # read in the unit file
  hyads.in <- read_fst( f, as.data.table = T)
  setnames( hyads.in, 'ZIP', 'zip')
  unames <- names( hyads.in)[ !(names( hyads.in) %in% 'zip')]
  
  # do the risk assessment for each unit
  if( !is.null( cuts.n)){
    ucuts <- split( unames, ceiling( seq_along( unames)/ cuts.n))
  } else
    ucuts <- list( unames)
  
  # fill in the hyads input
  hyads.in[, `:=` ( year = y) ] #factor( as( y, 'character')))]
  
  # merge with the annual data
  hyads.in.medi <- merge( dat_year_fill, hyads.in, by = c( 'zip', 'year'))
  hyads.in.medi.m <- melt( hyads.in.medi, id.vars = names( dat_year_fill))
  
  # trim down
  hyads.in.medi.m <- hyads.in.medi.m[ deaths.fill != 0 & 
                                        !is.na( deaths.fill)]
  
  # scale if called for
  if( !is.null( scale_by))
    hyads.in.medi.m[, value := value / scale_by]
  
  # print( hyads.in.medi.m)
  # do the risk assessment for each unit
  deaths_by_unit <- 
    pblapply( 
      ucuts,
      function( u, hyads_medi){
        hyads.in.medi.u <- hyads_medi[variable %in% u]
        
        # calculate deaths by each of the betas provided
        deaths_by_zip_unit_year <- 
          death_calculator( model_n = model_n,
                            betas = beta_ci.dt,
                            exposure_use = 'value',
                            fulldata = copy( hyads.in.medi.u),
                            by.vars = c( 'zip', 'STATE', 'variable')) 
        
        
        # sum each unit's deaths
        out.2 <- deaths_by_zip_unit_year[, .( 
          deaths_coef   = sum( deaths_coef, na.rm = T)),
          by = .( variable, beta, STATE)]
        
        return( out.2)
      }, 
      hyads_medi = hyads.in.medi.m[,.( zip, STATE, deaths.fill, denom.fill, value, variable)]) %>% rbindlist
  
  # dcast by which beta
  deaths_by_unit.out <- 
    dcast( deaths_by_unit,
           variable + STATE ~ beta,
           value.var = 'deaths_coef')
  setnames( deaths_by_unit.out, unique( deaths_by_unit$beta),
            paste0( 'deaths_coef_', unique( deaths_by_unit$beta)))
  
  # give year and month names
  deaths_by_unit.out[, `:=`( year = y)]
  
  
  return( deaths_by_unit.out)
}

# ============================================== #
# run the function
# ============================================== #
# 100gb, 10 core
# hyads, zip-code betas
deaths_by_year <- lapply( 1, #:nrow( zip.units.yr.dt), 
                          risk_assessmenter,
                          fyms = zip.units.yr.dt,
                          model_n = 'hyads',
                          cuts.n = 75) %>% rbindlist
write.fst(deaths_by_year, paste0(dir_out, "deaths_by_year_hyRR.fst"))

deaths_by_year_pm <- lapply( 1:nrow( zip.units.yr.dt), 
                             risk_assessmenter,
                             fyms = zip.units.yr.dt,
                             model_n = 'pm25_ensemble',
                             cuts.n = 75) %>% rbindlist
write.fst(deaths_by_year_pm, paste0(dir_out, "deaths_by_year_pmRR.fst"))

deaths_by_year_kr <- lapply( 1:nrow( zip.units.yr.dt), 
                             risk_assessmenter,
                             fyms = zip.units.yr.dt,
                             model_n = 'hyads',
                             cuts.n = 75) %>% rbindlist
write.fst(deaths_by_year_kr, paste0(dir_out, "deaths_by_year_krRR.fst"))

deaths_by_year_raw <- lapply( 1:nrow( zip.units.yr.dt_raw), 
                          risk_assessmenter,
                          fyms = zip.units.yr.dt_raw,
                          model_n = 'hyads_raw',
                          scale_by = 1318896214, # comes from mean(hyads_zips_tot_raw$hyads in 00_organize_data.R)
                          cuts.n = 75) %>% rbindlist
write.fst(deaths_by_year_raw, paste0(dir_out, "deaths_by_year_hy_rawRR.fst"))


