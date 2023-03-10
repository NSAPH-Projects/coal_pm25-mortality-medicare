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
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), as.data.table = TRUE)

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
cols_2000 <- names( dat_annual_hy)[ !( names( dat_annual_hy) %in% c( 'year', 'statecode', 'year_fac'))]
cols_2000diff <- names( dat_annual_hy)[ !( names( dat_annual_hy) %in% c( 'zip', 'year', 'statecode', 'year_fac'))]
dat_annual_use2000 <- 
  dat_annual_hy[year == 2000, ..cols_2000]
setnames( dat_annual_use2000,
          cols_2000diff,
          paste0( cols_2000diff, '_2000'))

# merge year 2000 data with original dataset
dat_annual_use_delt <- 
  merge( dat_annual_hy[ year %in% 1999:2016],
         dat_annual_use2000,
         by = 'zip')

# calculate differences in exposure & rates from 2000
dat_annual_use_delt[, `:=` ( death_rate_delta = death_rate - death_rate_2000,
                             pm25_ensemble_delta = pm25_ensemble - pm25_ensemble_2000,
                             Y1_delta = Y1 - Y1_2000,
                             Y1_raw_delta = Y1_raw - Y1_raw_2000,
                             mean_bmi_delta = mean_bmi - mean_bmi_2000,
                             smoke_rate_delta = smoke_rate - smoke_rate_2000,
                             hispanic_delta = hispanic - hispanic_2000,
                             pct_blk_delta = pct_blk - pct_blk_2000,
                             medhouseholdincome_delta = medhouseholdincome - medhouseholdincome_2000,
                             medianhousevalue_delta = medianhousevalue - medianhousevalue_2000,
                             poverty_delta  =  poverty - poverty_2000,
                             education_delta = education - education_2000,
                             popdensity_delta = popdensity - popdensity_2000,
                             pct_owner_occ_delta = pct_owner_occ - pct_owner_occ_2000,
                             summer_tmmx_delta = summer_tmmx - summer_tmmx_2000,
                             winter_tmmx_delta = winter_tmmx - winter_tmmx_2000,
                             summer_rmax_delta = summer_rmax - summer_rmax_2000,
                             winter_rmax_delta = winter_rmax - winter_rmax_2000,
                             ozone.current_year_delta = ozone.current_year - ozone.current_year_2000,
                             no2.current_year_delta = no2.current_year - no2.current_year_2000,
                             ozone_summer.current_year_delta = ozone_summer.current_year - ozone_summer.current_year_2000,
                             dead_delta = dead - dead_2000,
                             time_count_delta = time_count - time_count_2000)]

# center and scale
scale_cols_all <- grep( '_2000|_delt', names( dat_annual_use_delt), value = TRUE)
scale_cols <- scale_cols_all[ !( scale_cols_all %in% c( 'pm25_ensemble_2000', 'region_2000', 'death_rate_2000',
                                                        'dead_2000', 'Y1_2000', 'Y1_raw_2000', 'Y1.adj_2000',
                                                        'pm25_ensemble_delta', 'region_delta', 'death_rate_delta',
                                                        'dead_delta', 'Y1_delta', 'Y1.adj_delta'))]
scale_cols_names <- paste0( scale_cols, '_scale')

dat_annual_use_delt[, (scale_cols_names) := lapply( .SD, function( x) (x - mean( x, na.rm = T)) / sd( x, na.rm = T)),
                    .SDcols = scale_cols, 
                    by = year]


dat_annual_use_delt[ year > 2008, `:=` (Y1quantile = cut( Y1_delta, 
                                                          quantile(Y1_delta, probs = 0:20/20),
                                                          labels = FALSE, include.lowest = TRUE),
                                        Y1_rawquantile = cut( Y1_raw_delta, 
                                                              quantile(Y1_raw_delta, probs = 0:20/20),
                                                              labels = FALSE, include.lowest = TRUE),
                                        PMquantile = cut( pm25_ensemble_delta, 
                                                          quantile( pm25_ensemble_delta, 
                                                                    probs = 0:20/20, na.rm = TRUE),
                                                          labels = FALSE, include.lowest = TRUE)),
                     by = year]
## ======================================================= ##
# linear models of rates!
## ======================================================= ##
lm_fn <- 
  function( formula_in, x.name){
    
    lm1 <- lm( formula_in) #, weights = weights_in)
    print( summary( lm1))
    
    coef.names <- grepl( x.name, names( coef( lm1)))
    intercept <- coef( lm1)[ 1]
    slope <- coef( lm1)[ coef.names]
    se <- summary( lm1)$coefficients[ coef.names, 'Std. Error']
    return( list( slope = slope, se = se, intercept = intercept))
  }

# slopes over time
slopes_over_time_Y1 <- 
  dat_annual_use_delt[ year > 2008, 
                       .( Y1_delt_slope = 
                            lm_fn( death_rate_delta ~ -1 + as.factor( Y1quantile) + time_count_2000_scale +
                                     mean_bmi_2000_scale + smoke_rate_2000_scale + hispanic_2000_scale + pct_blk_2000_scale +
                                     medhouseholdincome_2000_scale + medianhousevalue_2000_scale +
                                     poverty_2000_scale + education_2000_scale + popdensity_2000_scale + pct_owner_occ_2000_scale +
                                     summer_tmmx_2000_scale + winter_tmmx_2000_scale + summer_rmax_2000_scale + winter_rmax_2000_scale +
                                     time_count_delta_scale +
                                     mean_bmi_delta_scale + smoke_rate_delta_scale + hispanic_delta_scale + pct_blk_delta_scale +
                                     medhouseholdincome_delta_scale + medianhousevalue_delta_scale +
                                     poverty_delta_scale + education_delta_scale + popdensity_delta_scale + pct_owner_occ_delta_scale +
                                     summer_tmmx_delta_scale + winter_tmmx_delta_scale + summer_rmax_delta_scale + winter_rmax_delta_scale +
                                     as.factor(region)
                                   , x.name = 'Y1quantile')$slope,
                          pollutant = 'Coal PM2.5' ),
                       by = year]
# slopes over time
# slopes_over_time_Y1_raw <- 
#   dat_annual_use_delt[ year > 2008, 
#                        .( Y1_delt_slope = 
#                             lm_fn( death_rate_delta ~ -1 + as.factor( Y1_rawquantile) + time_count_2000_scale +
#                                      mean_bmi_2000_scale + smoke_rate_2000_scale + hispanic_2000_scale + pct_blk_2000_scale +
#                                      medhouseholdincome_2000_scale + medianhousevalue_2000_scale +
#                                      poverty_2000_scale + education_2000_scale + popdensity_2000_scale + pct_owner_occ_2000_scale +
#                                      summer_tmmx_2000_scale + winter_tmmx_2000_scale + summer_rmax_2000_scale + winter_rmax_2000_scale +
#                                      time_count_delta_scale +
#                                      mean_bmi_delta_scale + smoke_rate_delta_scale + hispanic_delta_scale + pct_blk_delta_scale +
#                                      medhouseholdincome_delta_scale + medianhousevalue_delta_scale +
#                                      poverty_delta_scale + education_delta_scale + popdensity_delta_scale + pct_owner_occ_delta_scale +
#                                      summer_tmmx_delta_scale + winter_tmmx_delta_scale + summer_rmax_delta_scale + winter_rmax_delta_scale +
#                                      as.factor(region)
#                                    , x.name = 'Y1_rawquantile')$slope,
#                           pollutant = 'Raw coal PM2.5' ),
#                        by = year]
# slopes_over_time_PM <- 
#   dat_annual_use_delt[ year > 2008, 
#                        .( Y1_delt_slope = 
#                             lm_fn( death_rate_delta ~ -1 + as.factor( PMquantile) + time_count_2000_scale +
#                                      mean_bmi_2000_scale + smoke_rate_2000_scale + hispanic_2000_scale + pct_blk_2000_scale +
#                                      medhouseholdincome_2000_scale + medianhousevalue_2000_scale +
#                                      poverty_2000_scale + education_2000_scale + popdensity_2000_scale + pct_owner_occ_2000_scale +
#                                      summer_tmmx_2000_scale + winter_tmmx_2000_scale + summer_rmax_2000_scale + winter_rmax_2000_scale +
#                                      time_count_delta_scale +
#                                      mean_bmi_delta_scale + smoke_rate_delta_scale + hispanic_delta_scale + pct_blk_delta_scale +
#                                      medhouseholdincome_delta_scale + medianhousevalue_delta_scale +
#                                      poverty_delta_scale + education_delta_scale + popdensity_delta_scale + pct_owner_occ_delta_scale +
#                                      summer_tmmx_delta_scale + winter_tmmx_delta_scale + summer_rmax_delta_scale + winter_rmax_delta_scale +
#                                      as.factor(region)
#                                    , x.name = 'PMquantile')$slope,
#                           pollutant = 'PM2.5'),
#                        by = year]        

slopes_over_time <- 
  rbindlist( 
    list( slopes_over_time_Y1#,
          # slopes_over_time_Y1_raw,
          # slopes_over_time_PM
    ) )

slopes_over_time[, quant := rep( 1:20, 8)]


gg_first_diffs <- 
  ggplot( slopes_over_time,
          aes( x = quant, y = Y1_delt_slope * 10000)) + 
  geom_hline( yintercept = 0) +
  geom_smooth( method = 'lm') +
  labs( y = 'Change in mortality rate per 10,000 since 2000',
        x = 'Exposure change quantile (larger reductions to the left)') +
  geom_point() + 
  facet_grid( . ~ year) + 
  expand_limits( y = 0) +
  theme_minimal() +
  theme( axis.text = element_text( size = 12),
         axis.text.x = element_blank(),
         axis.title = element_text( size = 14),
         strip.text = element_text( size = 14))
ggsave( gg_first_diffs,
        filename = 'figures/first_diffs_since_2000.pdf',
        height = 5, width = 10, unit = 'in', device = cairo_pdf)


## ======================================================= ##
# plots of death rates in areas with highest ventiles
## ======================================================= ##
large_zips <- unique( dat_annual_use_delt[ Y1quantile == 1 & time_count_2000 > 1000]$zip)
length( large_zips)

dat_annual_use_delt[, state_name := state.name[which( state.abb == statecode)]]

labeller_state <- 
  function( abb){
   df.out <-  lapply( abb[,1], function(a) {
      state.name[which( state.abb == a)] %>% data.frame
    }) %>% rbindlist()
   colnames( df.out) <- 'statecode'
   return( df.out)
  }

gg_zips_rates <- 
  ggplot( dat_annual_use_delt[zip %in% large_zips],
          aes( y = death_rate, x = year, group = zip#, linewidth = time_count_2000
               )) + 
  geom_line( alpha = .1) + 
  scale_size_continuous( 
    name = 'Number of Medicare enrollees',
    range = c( 0.05, 2)) +
  facet_wrap( . ~ statecode, labeller = labeller_state) + 
  labs( y = 'Annual death rate', x = 'Year',
        linewidth = 'Medicare beneficiaries') +
  expand_limits( y = 0) +
  guides( linewidth = guide_legend(title.position="top")) +
  theme_minimal() + 
  theme( axis.title = element_text( size = 14),
         axis.text = element_text( size = 12),
         legend.direction = 'horizontal',
         legend.position = c( .6, .08),
         legend.title = element_text( size = 14),
         legend.text = element_text( size = 12),
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text = element_text( size = 14))
ggsave( gg_zips_rates,
        filename = 'figures/zip_death_rates.pdf',
        height = 8, width = 10, unit = 'in', device = cairo_pdf)


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

dat_annual_early[, `:=` (Y1quantile = cut(Y1_early, quantile(Y1_early, probs = 0:20/20),
                                          labels = FALSE, include.lowest = TRUE),
                         PMquantile = cut(pm25_ensemble_early, quantile(pm25_ensemble_early, 
                                                                        probs = 0:20/20, na.rm = TRUE),
                                          labels = FALSE, include.lowest = TRUE))]

dat_annual_avgs <- 
  merge( dat_annual_early,
         dat_annual_late, by = 'zip') %>%
  merge( )

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


