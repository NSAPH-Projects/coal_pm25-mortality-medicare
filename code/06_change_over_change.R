library( fst)
library( data.table)
library( magrittr)
library( ggplot2)

# Xiao's code
# https://github.com/NSAPH/National_Causal/blob/master/statistical_models.R

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/data' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

aggregate_data <- read.fst( file.path(dir_data, "cache_data", "aggregate_data.fst"),
                            as.data.table = TRUE)
covariates <- 
  read.fst( file.path(dir_data, "cache_data", "covariates.fst"), as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj'), as.data.table = TRUE)

## ======================================================= ##
# create death & rate sums by zip-year
## ======================================================= ##
dat_annual_zip_death <- 
  aggregate_data[, .( death_rate = sum( dead) / sum( time_count),
                      dead = sum( dead),
                      time_count = sum( time_count)),
                 by = .( year, zip)]

dat_annual_zip_year <- merge(covariates, 
                             dat_annual_zip_death, by = c("zip", "year")) #, all.x = TRUE)

dat_annual_hy <- 
  merge( dat_annual_zip_year,
         dat_annual,
         by = c( 'zip', 'year'))

# clean house
rm(list = c("aggregate_data", "covariates", "dat_annual"))

## ======================================================= ##
# calculcate differences from 2000
## ======================================================= ##
# extract year 2000 data
dat_annual_use2000 <- 
  dat_annual_hy[year == 2000, .( zip, death_rate, pm25_ensemble, Y1)]
setnames( dat_annual_use2000,
          c( 'death_rate', 'pm25_ensemble', 'Y1'),
          c( 'death_rate2000', 'pm25_ensemble2000', 'Y12000'))

# merge year 2000 data with original dataset
dat_annual_use_delt <- 
  merge( dat_annual_hy,
         dat_annual_use2000,
         by = 'zip')

# calculate differences in exposure & rates from 2000
dat_annual_use_delt[, `:=` ( death_rate_delta = death_rate - death_rate2000,
                             log_death_rate_delta = log( death_rate) - log( death_rate2000),
                             pm25_ensemble_delta = pm25_ensemble - pm25_ensemble2000,
                             Y1_delta = Y1 - Y12000)]

## ======================================================= ##
# a slightly different approach—averages from 2000 - 2004 and 2012-2016
## ======================================================= ##
dat_annual_early <- 
  dat_annual_hy[year %in% 2000:2004, 
                .( death_rate_early = sum( dead) / sum( time_count), 
                   death_rate_mean_early = mean( death_rate), 
                   pm25_ensemble_early = mean( pm25_ensemble), 
                   Y1_early = mean( Y1)),
                by = 'zip']
dat_annual_late <- 
  dat_annual_hy[year %in% 2012:2016, 
                .( death_rate_late = sum( dead) / sum( time_count), 
                   death_rate_mean_late = mean( death_rate), 
                   pm25_ensemble_late = mean( pm25_ensemble), 
                   Y1_late = mean( Y1)),
                by = 'zip']

dat_annual_early[, `:=` (Y1quantile = cut(Y1_early, quantile(Y1_early, probs = 0:10/10),
                                          labels = FALSE, include.lowest = TRUE),
                         PMquantile = cut(pm25_ensemble_early, quantile(pm25_ensemble_early, 
                                                                        probs = 0:10/10, na.rm = TRUE),
                                          labels = FALSE, include.lowest = TRUE))]

dat_annual_avgs <- 
  merge( dat_annual_early,
         dat_annual_late, by = 'zip')

# calculate differences in exposure & rates from 2000
dat_annual_avgs[, `:=` ( death_rate_delta = death_rate_late - death_rate_early,
                         death_rate_mean_delta = death_rate_mean_late - death_rate_mean_early,
                         pm25_ensemble_delta = pm25_ensemble_late - pm25_ensemble_early,
                         Y1_delta = Y1_late - Y1_early)]

# merge year 2000 data with original dataset
dat_annual_avgs_delt <- 
  merge( dat_annual_hy[year == 2000],
         dat_annual_avgs,
         by = 'zip')

# calculate average change at different quantiles of early exposure
dat_annual_avgs_quant <- 
  dat_annual_avgs[, .( delta_death = mean( death_rate_late - death_rate_early),
                       delta_Y1 = mean( Y1_late - Y1_early)
  ), by = Y1quantile]
ggplot( dat_annual_avgs_quant,
        aes( x = delta_Y1, y = delta_death)) + 
  geom_point()

dat_annual_avgs_quant_PM <- 
  dat_annual_avgs[, .( delta_death = mean( death_rate_late - death_rate_early),
                       delta_PM = mean( pm25_ensemble_late - pm25_ensemble_early, na.rm = T)
  ), by = PMquantile]
ggplot( dat_annual_avgs_quant_PM,
        aes( x = delta_PM, y = delta_death)) + 
  geom_point()

#evidence of this being driven by space and not time?
# if average of ∆death rate is ~-0.007, would it be possible to 
# see the influence of a RR of 1.5% from coal PM?
# average baseline death rate is 0.05 death/enrolled
# average change in death rate is 0.007 death/enrolled
# so, we're looking at changes of 0.007 / 0.05 = 14%, an order of 
# magnitude larger than the change expected from coal PM
## ======================================================= ##
# plots!
## ======================================================= ##

ggplot( dat_annual_use_delt[year %in% c( 2005, 2010, 2016)],
        aes( x = Y1_delta,
             y = death_rate_delta)) + 
  geom_point( alpha = .1) + 
  geom_smooth() +
  facet_wrap( . ~ year) + 
  theme_bw()

ggplot( dat_annual_avgs,
        aes( x = Y1_delta,
             y = death_rate_delta)) + 
  geom_point( alpha = .1) + 
  geom_smooth() +
  theme_bw()

## ======================================================= ##
# linear models of rates!
## ======================================================= ##
lm_fn <- 
  function( formula_in, x.name){
    # print( formula_in)
    
    lm1 <- lm( formula_in) #, weights = weights_in)
    slope <- coef( lm1)[x.name]
    se <- summary( lm1)$coefficients[x.name, 'Std. Error']
    return( list( slope = slope, se = se))
  }

# slopes over time
slopes_over_time <- 
  dat_annual_use_delt[ year > 2000 ,#& time_count > 1000, 
                       .( Y1_delt_slope = 
                            lm_fn( death_rate_delta ~ Y1_delta +
                                     mean_bmi + smoke_rate + hispanic + pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(region), x.name = 'Y1_delta')$slope,
                          Y1_delt_se =
                            lm_fn( death_rate_delta ~ Y1_delta +
                                     mean_bmi + smoke_rate + hispanic + pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(region), x.name = 'Y1_delta')$se,
                          PM_delt_slope =
                            lm_fn( death_rate_delta ~ pm25_ensemble_delta +
                                     mean_bmi + smoke_rate + hispanic + pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(region), x.name = 'pm25_ensemble_delta')$slope,
                          PM_delt_se =
                            lm_fn( death_rate_delta ~ pm25_ensemble_delta +
                                     mean_bmi + smoke_rate + hispanic + pct_blk +
                                     medhouseholdincome + medianhousevalue +
                                     poverty + education + popdensity + pct_owner_occ +
                                     summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                     as.factor(region), x.name = 'pm25_ensemble_delta')$se
                       ),
                       by = year]


slopes_over_time_avgs <- 
  dat_annual_avgs_delt[ ,#time_count > 1000, 
                        .( Y1_delt_slope = 
                             lm_fn( death_rate_delta ~ Y1_delta +
                                      mean_bmi + smoke_rate + hispanic + pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(region), x.name = 'Y1_delta')$slope,
                           Y1_delt_se =
                             lm_fn( death_rate_delta ~ Y1_delta +
                                      mean_bmi + smoke_rate + hispanic + pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(region), x.name = 'Y1_delta')$se,
                           PM_delt_slope =
                             lm_fn( death_rate_delta ~ pm25_ensemble_delta +
                                      mean_bmi + smoke_rate + hispanic + pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(region), x.name = 'pm25_ensemble_delta')$slope,
                           PM_delt_se =
                             lm_fn( death_rate_delta ~ pm25_ensemble_delta +
                                      mean_bmi + smoke_rate + hispanic + pct_blk +
                                      medhouseholdincome + medianhousevalue +
                                      poverty + education + popdensity + pct_owner_occ +
                                      summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
                                      as.factor(region), x.name = 'pm25_ensemble_delta')$se
                        )]

summary( lm( death_rate_mean_delta ~ Y1_delta +
               mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(region), data = dat_annual_avgs_delt, 
             weights = time_count))
summary( lm( death_rate_delta ~ pm25_ensemble_delta+
               mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(region), data = dat_annual_avgs_delt, weights = time_count))

slopes_over_time.m <- 
  melt( slopes_over_time,
        measure.vars = list( 'slope' = c( 'Y1_delt_slope', 'PM_delt_slope'),
                             'se' = c( 'Y1_delt_se', 'PM_delt_se')),
        id.vars = 'year')

ggplot( slopes_over_time.m,
        aes( x = as.factor( year), y = slope, color = variable,
             ymin = slope - se, ymax = slope + se)) + 
  geom_errorbar( position = position_dodge( width = .3), width = 0) + 
  geom_point( position = position_dodge( width = .3)) + 
  scale_color_brewer( palette = 'Dark2', 
                      labels = c( '1' = 'coal PM2.5',
                                  '2' = 'PM2.5'))



