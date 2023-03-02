library( fst)
library( data.table)
library( gnm)
library( magrittr)

# Xiao's code
# https://github.com/NSAPH/National_Causal/blob/master/statistical_models.R

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/data' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

aggregate_data <- read.fst( file.path(dir_data, "cache_data", "aggregate_data.fst"),
                            as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), as.data.table = TRUE)
dat_annual_use <- merge(aggregate_data, dat_annual, by = c("zip", "year")) #, all.x = TRUE)

# 35 locations are lost in hyads because they are points
# most do not span the entire time series
# unique( dat_annual_use[is.na( Y1)]$zip)
# dat_annual_use[zip == '01477']
# dat_annual_use[zip == '01517']
# dat_annual_use[zip == '02565']
# dat_annual_use[zip == '22201']

# clean house
rm(list = c("aggregate_data", "dat_annual"))

# calculate portion of PM that is not coal PM
dat_annual_use[, pm_not_coalPM := pm25_ensemble - Y1]

## ======================================================= ##
# run xiao's model
# https://github.com/wxwx1993/National_Causal/blob/master/statistical_models.R
## ======================================================= ##
# Cox-equvalent conditional Poisson Regression
## PM-only model
gnm_raw <- 
  gnm(dead ~ pm25_ensemble +
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson <- coef(gnm_raw)
rm( gnm_raw)
gc()

## HyADS-only model
gnm_raw_hy <- 
  gnm(dead ~ Y1 +
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy <- coef(gnm_raw_hy)
rm( gnm_raw_hy)
gc()

# model just based on hyads adjusted to sulfate
gnm_raw_hy.adj <- 
  gnm(dead ~ Y1.adj +
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) +
        as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy.adj <- coef( gnm_raw_hy.adj)
rm( gnm_raw_hy.adj)
gc()

## HyADS model adjusted by [total PM - coal PM], 2003 and earlier
gnm_raw_hy_early <- 
  gnm(dead ~ Y1 + 
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use[ year <= 2003],
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_early <- coef(gnm_raw_hy_early)
rm( gnm_raw_hy_early)
gc()

## HyADS model adjusted by [total PM - coal PM], 2004-2007
gnm_raw_hy_mid <- 
  gnm(dead ~ Y1 + 
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use[ year >= 2004 & year <= 2007],
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_mid <- coef(gnm_raw_hy_mid)
rm( gnm_raw_hy_mid)
gc()

## HyADS model adjusted by [total PM - coal PM], 2008 and later
gnm_raw_hy_late <- 
  gnm(dead ~ Y1 + 
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use[ year >= 2008],
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_late <- coef(gnm_raw_hy_late)
rm( gnm_raw_hy_late)
gc()

## HyADS-only model, only at ZIPs with non-missing PM
# gnm_raw_hy_pmmissing <- 
#   gnm(dead ~ Y1 +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use[!is.na( pm25_ensemble)],
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_pmmissing <- coef(gnm_raw_hy_pmmissing)
# rm( gnm_raw_hy_pmmissing)
# gc()

## ====================================================
#  adjusted models - total PM
## ====================================================
## HyADS model adjusted by total PM
gnm_raw_hy_pm <- 
  gnm(dead ~ Y1 + pm25_ensemble +
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_pm <- coef(gnm_raw_hy_pm)
rm( gnm_raw_hy_pm)
gc()

# model just based on hyads adjusted to sulfate with pm adjustment
# gnm_raw_hy.adj_pm <- 
#   gnm(dead ~ Y1.adj + pm25_ensemble +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) +
#         as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use,
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy.adj_pm <- coef( gnm_raw_hy.adj_pm)
# rm( gnm_raw_hy.adj_pm)
# gc()

## ====================================================
#  adjusted models - [total PM - coal PM]
## ====================================================
## HyADS model adjusted by [total PM - coal PM]
gnm_raw_hy_dpm <- 
  gnm(dead ~ Y1 + pm_not_coalPM +
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_dpm <- coef(gnm_raw_hy_dpm)
rm( gnm_raw_hy_dpm)
gc()

## HyADS model adjusted by [total PM - coal PM], 2003 and earlier
# gnm_raw_hy_dpm_early <- 
#   gnm(dead ~ Y1 + pm_not_coalPM +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use[ year <= 2003],
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_dpm_early <- coef(gnm_raw_hy_dpm_early)
# rm( gnm_raw_hy_dpm_early)
# gc()

## HyADS model adjusted by [total PM - coal PM], 2004-2007
# gnm_raw_hy_dpm_mid <- 
#   gnm(dead ~ Y1 + pm_not_coalPM +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use[ year >= 2004 & year <= 2007],
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_dpm_mid <- coef(gnm_raw_hy_dpm_mid)
# rm( gnm_raw_hy_dpm_mid)
# gc()

## HyADS model adjusted by [total PM - coal PM], 2008 and later
# gnm_raw_hy_dpm_late <- 
#   gnm(dead ~ Y1 + pm_not_coalPM +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use[ year >= 2008],
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_dpm_late <- coef(gnm_raw_hy_dpm_late)
# rm( gnm_raw_hy_dpm_late)
# gc()

# model just based on hyads adjusted to sulfate with [total PM - coal PM] adjustment
# gnm_raw_hy.adj_dpm <- 
#   gnm(dead ~ Y1.adj + pm_not_coalPM +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) +
#         as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use,
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy.adj_dpm <- coef( gnm_raw_hy.adj_dpm)
# rm( gnm_raw_hy.adj_dpm)
# gc()

## ====================================================
#  adjusted models - [total PM - coal PM] and other pollutants
## ====================================================
## HyADS model adjusted by [total PM - coal PM] and other pollutants
# gnm_raw_hy_pollutants <- 
#   gnm(dead ~ Y1 + pm_not_coalPM + ozone.current_year +
#         no2.current_year + ozone_summer.current_year +
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use,
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_pollutants <- coef(gnm_raw_hy_pollutants)
# rm( gnm_raw_hy_pollutants)
# gc()

## HyADS model adjusted by [total PM - coal PM] and other pollutants
gnm_raw_hy_dpm_no2 <- 
  gnm(dead ~ Y1 + pm_not_coalPM + no2.current_year + 
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_dpm_no2 <- coef( gnm_raw_hy_dpm_no2)
rm( gnm_raw_hy_dpm_no2)
gc()

## HyADS model adjusted by [total PM - coal PM] and other pollutants
gnm_raw_hy_no2 <- 
  gnm(dead ~ Y1 + no2.current_year + 
        mean_bmi + smoke_rate + hispanic + pct_blk +
        medhouseholdincome + medianhousevalue +
        poverty + education + popdensity + pct_owner_occ +
        summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
        as.factor(year) + as.factor(region) +
        offset(log(time_count)),
      eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
      data = dat_annual_use,
      family = poisson(link = "log"))

# save only the summaries, clear large object
Poisson_hy_no2 <- coef( gnm_raw_hy_no2)
rm( gnm_raw_hy_no2)
gc()

## raw HyADS model adjusted by [total PM - coal PM] and other pollutants
# gnm_raw_hy_raw_dpm_no2 <- 
#   gnm(dead ~ Y1_raw + pm_not_coalPM + no2.current_year + 
#         mean_bmi + smoke_rate + hispanic + pct_blk +
#         medhouseholdincome + medianhousevalue +
#         poverty + education + popdensity + pct_owner_occ +
#         summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#         as.factor(year) + as.factor(region) +
#         offset(log(time_count)),
#       eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
#       data = dat_annual_use,
#       family = poisson(link = "log"))

# save only the summaries, clear large object
# Poisson_hy_raw_dpm_no2 <- coef(gnm_raw_hy_raw_dpm_no2)
# rm( gnm_raw_hy_raw_dpm_no2)
# gc()

## ====================================================
#  view/manipulate
## ====================================================
# view the coefficients
exp(1 * Poisson['pm25_ensemble'])
exp(1 * Poisson_hy['Y1'])
exp(1 * c( Poisson_hy_pm['Y1'], Poisson_hy_pm['pm25_ensemble']))
exp(1 * c( Poisson_hy_dpm['Y1'], Poisson_hy_dpm['pm_not_coalPM']))
exp(1 * c( Poisson_hy_dpm_early['Y1'], Poisson_hy_dpm_early['pm_not_coalPM']))
exp(1 * Poisson_hy_pmmissing['Y1'])
exp(1 * Poisson_hy.adj['Y1.adj'])
exp(1 * c( Poisson_hy.adj_pm['Y1.adj'], Poisson_hy.adj_pm['pm25_ensemble']))
exp(1 * c( Poisson_hy_pollutants['Y1'], Poisson_hy_pollutants['pm_not_coalPM'], 
           Poisson_hy_pollutants['no2.current_year'], Poisson_hy_pollutants['ozone_summer.current_year']))



# save the coefficients as data.table
poisson_coefs_vec <- 
  c( pm25_ensemble   = unname( Poisson['pm25_ensemble']),
     hyads =              unname( Poisson_hy['Y1']),
     hyads_early =        unname( Poisson_hy_early['Y1']),
     hyads_mid   =        unname( Poisson_hy_mid['Y1']),
     hyads_late =         unname( Poisson_hy_late['Y1']),
     hyads_pm =           unname( Poisson_hy_pm['Y1']),
     hyads_dpm =          unname( Poisson_hy_dpm['Y1']),
     hyads_dpm_early =    unname( Poisson_hy_dpm_early['Y1']),
     hyads_dpm_mid   =    unname( Poisson_hy_dpm_mid['Y1']),
     hyads_dpm_late =     unname( Poisson_hy_dpm_late['Y1']),
     hyads_pollutants =   unname( Poisson_hy_pollutants['Y1']),
     Poisson_hy_dpm_no2= unname( Poisson_hy_dpm_no2['Y1']),
     hyads_raw_dpm_no2=   unname( Poisson_hy_raw_dpm_no2['Y1_raw']),
     hyads.adj =       unname( Poisson_hy.adj['Y1.adj']),
     hyads.adj_pm =    unname( Poisson_hy.adj_pm['Y1.adj']),
     hyads.adj_dpm =   unname( Poisson_hy.adj_dpm['Y1.adj']))

poisson_coefs <- 
  data.table( Estimate = poisson_coefs_vec,
              model = names( poisson_coefs_vec))

# save the coefficients
write.csv( poisson_coefs, file = paste0(dir_out, "poisson_model_coefs.csv"))

# save the summaries
# save(Poisson, file = paste0(dir_out, "Poisson_pm25.RData"))
# save(Poisson_hy, file = paste0(dir_out, "Poisson_hyads.RData"))

## ==================================================== ##
##  table of zip code characteristics
## ==================================================== ##
# load the covariate data
covariates <- read.fst( paste0(dir_data, "covariates.fst"),
                        as.data.table = TRUE)

# merge with hyads
dat_annual_covars <- 
  merge(covariates, dat_annual, 
        by = c("zip", "year")) 

#  choose covariates to summarize
covars <- c( 'mean_bmi', 'smoke_rate', 'hispanic', 'pct_blk',
             'medhouseholdincome', 'medianhousevalue',
             'poverty', 'education', 'popdensity', 'pct_owner_occ',
             'summer_tmmx', 'winter_tmmx', 'summer_rmax', 'winter_rmax',
             'Y1'
)

# create a table with some statistics
dat_annual_covars[,..covars] %>% 
  psych::describe(na.rm = TRUE) %>% 
  select(mean, sd, median, min, max) %>% 
  knitr::kable(digits = 2)

# calculate stats for individual level characteristics
# define three race variables
all_count <- sum( dat_annual_use$time_count)
dat_annual_female <- 
  dat_annual_use[  sex == 2, sum( time_count)] / all_count
dat_annual_white <- 
  dat_annual_use[  race == 1, sum( time_count)] / all_count
dat_annual_black <- 
  dat_annual_use[  race == 2, sum( time_count)] / all_count
dat_annual_dual <- 
  dat_annual_use[  dual == 1, sum( time_count)] / all_count



## ==================================================== ##
## check out biggest impacting facilities
## ==================================================== ##



