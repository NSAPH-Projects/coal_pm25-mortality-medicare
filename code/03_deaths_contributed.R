#renv::install("lhenneman/disperseR@dev")

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
                        columns = c('zip','year', 'Y1', 'Y1.adj'), 
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
  poisson_coefs[ model == 'Y1']$Estimate +
  c( -boot_beta_hy, 0, +boot_beta_hy)
betas_hy.adj <- 
  poisson_coefs[ model == 'hyads.adj']$Estimate +
  c( -boot_beta_hy.adj, 0, +boot_beta_hy.adj)
betas_kr <- log( c( 1.035, 1.056, 1.078)) / 10

exp(1 * betas_pm)
exp(1 * betas_hy)
exp(1 * betas_hy.adj)
exp(1 * betas_kr)

( exp(1 * betas_hy) - 1 ) / ( exp(1 * betas_pm) - 1)

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
           'data/cache_data/dat_year_fill.fst')
## ====================================================== ##
# sum deaths by state & year as input for adjoint comparison
## ====================================================== ##
sum_deaths_year <- dat_year_fill[, .( deaths.fill = sum( deaths.fill, na.rm = T),
                                      denom.fill = sum( denom.fill, na.rm = T)),
                                 by = .( year, STATE)]
write.fst( sum_deaths_year,
           'data/cache_data/total_deaths_by_state_year.fst')
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
#hyads HR
deaths_by_year_all.hy <- 
  lapply( betas_hy,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_hy == b))]
          }) %>% rbindlist

# PM HR (from Xiao Wu Science Advances)
deaths_by_year_all.pm <- 
  lapply( betas_pm,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_pm == b))]
          }) %>% rbindlist

# Krueski HR
deaths_by_year_all.kr <- 
  lapply( betas_kr,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_kr == b))]
          }) %>% rbindlist

# sulfate-adjusted coal PM2.5
deaths_by_year_all.hy_adj <- 
  lapply( betas_hy.adj,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1.adj) - 1) * denom.fill#,
              ), by = .( zip, year, STATE)]
            deaths_by_zip.b[, beta := as.character( which( betas_hy.adj == b))]
          }) %>% rbindlist



# add model labels
deaths_by_year_all.hy[, model := 'hyads']
deaths_by_year_all.pm[, model := 'pm25_ensemble']
deaths_by_year_all.kr[, model := 'pm25_krewski']
deaths_by_year_all.hy_adj[, model := 'hyads.adj']

# combine into single data.table
deaths_by_year_all <- 
  rbind( deaths_by_year_all.hy,
         deaths_by_year_all.pm,
         deaths_by_year_all.kr,
         deaths_by_year_all.hy_adj)

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
deaths_by_year_all.c[ model == 'hyads', sum( deaths_coef_2)]
deaths_by_year_all.c[ model == 'hyads.adj', sum( deaths_coef_1)]
deaths_by_year_all.c[ model == 'hyads.adj', sum( deaths_coef_2)]
deaths_by_year_all.c[ model == 'hyads.adj', sum( deaths_coef_3)]

# save it
write.fst( deaths_by_year_all.c, paste0(dir_out, "deaths_by_year_total_coal_pm25.fst"))

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
write.fst( deaths_by_year_state_all.pm25_ensemble.c, paste0(dir_out, "deaths_by_year_total_pm25_ensemble.fst"))

## ====================================================== ##
# calculate deaths avoided for hyads; constant death rate or constant hyads
## ====================================================== ##
deaths_all.hyads_constdeath <- 
  lapply( betas_hy,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill_99 / denom.fill_99 * ( exp( b * Y1) - 1) * denom.fill_99#,
              ), by = .( year)]
            deaths_by_zip.b[, beta := as.character( which( betas_hy == b))]
          }) %>% rbindlist

deaths_all.hyads_consthy <- 
  lapply( betas_hy,
          function( b){
            deaths_by_zip.b <- 
              dat_year_fill[, .( 
                deaths_coef = deaths.fill / denom.fill * ( exp( b * Y1_99) - 1) * denom.fill#,
              ), by = .( year)]
            deaths_by_zip.b[, beta := as.character( which( betas_hy == b))]
          }) %>% rbindlist

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
# set directory structure
disperseR.base <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/'
disperseR::create_dirs( disperseR.base)

# define pm25 exposure directory
exp25_dir <- '/nfs/home/H/henneman/shared_space/ci3_nsaph/LucasH/disperseR/main/output/zips_model.lm.cv_single_poly'

# get files for unit-specific hyads
zip.units.yr <- list.files( exp25_dir,
                            pattern = 'zips_.*byunit.*\\d{4}\\.fst',
                            full.names = TRUE)

# put together data.table for easy reference
zip.units.yr.dt <- data.table( f = zip.units.yr,
                               y = 1999:2020)
## ====================================================== ##
# do the risk assessment
## ====================================================== ##
# write the function
# cuts defines number of units run at a time
# cuts = 2 seems to work okay for 100 samples per zip
# cuts = 10 seems to work okay for 3 samples per zip
risk_assessmenter <- function( n, fyms, betas, cuts.n = 2){
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
  # print( hyads.in.medi.m)
  # do the risk assessment for each unit
  deaths_by_unit <- 
    pbmclapply( 
      ucuts,
      function( u, hyads_medi, betas_use){
        print( u)
        hyads.in.medi.u <- hyads_medi[variable %in% u]
        
        # calculate deaths by each of the betas provided
        deaths_by_zip_unit_year <- 
          lapply( betas_use,
                  function( b){
                    deaths_by_zip_unit.b <- 
                      hyads.in.medi.u[, .( 
                        deaths_coef = deaths.fill / denom.fill * ( exp( b * value) - 1) * denom.fill#,
                        # deaths_coef = deaths.fill / denom.fill * (1 - (1 / exp( b)^value)) * denom.fill#,
                      ), by = .( zip, STATE, variable)]
                    deaths_by_zip_unit.b[, beta := as.character( which( betas_use == b))]
                  }) %>% rbindlist
        
        # sum each unit's deaths
        out.2 <- deaths_by_zip_unit_year[, .( 
          deaths_coef   = sum( deaths_coef, na.rm = T)),
          by = .( variable, beta, STATE)]
        
        return( out.2)
      }, 
      hyads_medi = hyads.in.medi.m[,.( zip, STATE, deaths.fill, denom.fill, value, variable)], 
      betas_use = betas, 
      mc.cores = NSAPHutils::get_cpus()) %>% rbindlist
  
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
deaths_by_year <- lapply( 1:nrow( zip.units.yr.dt), 
                          risk_assessmenter,
                          fyms = zip.units.yr.dt,
                          betas = betas_hy,
                          cuts.n = 3) %>% rbindlist
write.fst(deaths_by_year, paste0(dir_out, "deaths_by_year_hyRR.fst"))

deaths_by_year_pm <- lapply( 1:nrow( zip.units.yr.dt), 
                             risk_assessmenter,
                             fyms = zip.units.yr.dt,
                             betas = betas_pm,
                             cuts.n = 3) %>% rbindlist
write.fst(deaths_by_year_pm, paste0(dir_out, "deaths_by_year_pmRR.fst"))

deaths_by_year_kr <- lapply( 1:nrow( zip.units.yr.dt), 
                             risk_assessmenter,
                             fyms = zip.units.yr.dt,
                             betas = betas_kr,
                             cuts.n = 3) %>% rbindlist
write.fst(deaths_by_year_kr, paste0(dir_out, "deaths_by_year_krRR.fst"))


