library( fst)
library( data.table)
library( gnm)
library( magrittr)

# Xiao's code
# https://github.com/NSAPH/National_Causal/blob/master/statistical_models.R

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

aggregate_data <- read.fst( file.path(dir_data, "cache_data", "aggregate_data.fst"),
                            as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj'), as.data.table = TRUE)
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

## ======================================================= ##
# run xiao's model
# https://github.com/wxwx1993/National_Causal/blob/master/statistical_models.R
## ======================================================= ##
# Cox-equvalent conditional Poisson Regression
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

# model just based on hyads adjusted to sulfate
gnm_raw_hy_adj <- 
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
Poisson_hy.adj <- coef( gnm_raw_hy_adj)
rm( gnm_raw_hy_adj)

# view the coefficients
exp(10 * Poisson['pm25_ensemble'])
exp(10 * Poisson_hy['Y1'])
exp(1 * Poisson_hy.adj['Y1.adj'])

exp( 10 * 0.02454625)

# save the coefficients as data.table
poisson_coefs <- rbind( Poisson['pm25_ensemble'],
                        Poisson_hy['Y1'],
                        Poisson_hy.adj['Y1.adj']) %>%
  as.data.table
setnames( poisson_coefs, 'pm25_ensemble', 'Estimate')
poisson_coefs[, model := c( 'pm25_ensemble', 'hyads', 'hyads.adj')]

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





