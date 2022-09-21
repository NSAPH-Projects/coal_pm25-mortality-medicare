library( data.table)
library( raster)
library( sf)
library( fst)
library( ggplot2)
library( magrittr)
## need medicare deaths by state (and total 48 states US)
# need for the Cox model constant hyads & PM

## =========================================================== ##
## Load adjoint runs, perform rankings, merge
## =========================================================== ##

read.adj <- function( base.dir = 'data/adjoint_results'){
  
  # ID directory, list files
  adj.dir <- file.path( base.dir)
  adj.files <- data.table( file = list.files( adj.dir,
                                              pattern = '_adjoint_',
                                              full.names = T))
  
  # ID states and years
  adj.files[, `:=` ( state = gsub( '.*_adjoint_|.csv', 
                                   '', 
                                   file),
                     year = gsub( '.*_all_units_|_adjoint.*.csv', 
                                  '', 
                                  file))]
  adj.files[, state := gsub( 'USA', 'US', state)]
  
  # read in file that includes names, assign these names to read files
  adj.withnames <- fread( file.path( base.dir, 'final_merge_nei_ampd_all_units_2006_adjoint_USA.csv'))
  adj.names <- names( adj.withnames)
  
  # Read in, listify
  adj.in <- lapply( 1:nrow( adj.files),
                    function( n, dt, adj.names.in) { 
                      x <- dt[n]
                      adj.in.n <- fread( x$file)
                      names( adj.in.n) <- adj.names.in
                      
                      adj.in.n[, `:=` (state = x$state,
                                       year = as( x$year,
                                                  'numeric'))]
                      return( adj.in.n)
                    },
                    adj.files,
                    adj.names)
  adj.in.dt <- rbindlist( adj.in)
  
  # minor adjustments, rank by SOx exposure
  adj.in.dt[, uID := gsub( ' |"', '', ID)]
  adj.in.dt[, uID := gsub('_|-|\\*', '.', uID)]
  
  # minor adjustments, rank by SOx exposure
  adj.in.dt <- adj.in.dt[ SOx > 0]
  
  # set names to recognizable things
  setnames( adj.in.dt,
            c( 'PM2.5 pop exposure due to NOx ems',
               'PM2.5 pop exposure due to SOx ems',
               'Total PM2.5 pop exposure'),
            c( 'nox.pm25.popwt', 'sox.pm25.popwt', 'tot.pm25.popwt'))
  
  return( adj.in.dt)
}

## ===============================================
#  Read in the Irene geoschem data
## ===============================================
irene_adj <- read.adj( )[state == 'US']

## ================================================== ##
# read population data and merge with adjoint
## ================================================== ##
us_states.pop.dt.m <- fread( 'data/us_population_by_state_2006_2011.csv',
                             drop = 'V1')

## merge populations with adjoint
irene_adj_pop <- merge( irene_adj, us_states.pop.dt.m, by = c( 'state', 'year'))

# replace populations from Irene's email
irene_adj_pop[ year == 2006, pop := 290345000] #298345000]
irene_adj_pop[ year == 2011, pop := 307602241]

## ================================================== ##
# merge with unit information
## ================================================== ##
# read in facility info
info.list <- c( 'Facility.Name', 'Facility.ID..ORISPL.', 'Unit.ID', 'State')
facility_info <- fread( '~/shared_space/ci3_nsaph/LucasH/RFM_health/inputs/Coal_facilities/AMPD_Unit.csv',
                        select = info.list) %>% unique
setnames( facility_info, 
          c( 'Facility.ID..ORISPL.', 'Unit.ID'), 
          c( 'FacID', 'unitID'))

facility_info[, `:=` ( #state = state.name[ match( State, state.abb)],
  uID = gsub( '-|\\#|\\*', '.', paste( FacID, unitID, sep = '.')))]

# merge with adjoin data
irene_adj_pop_fac <- merge( irene_adj_pop, facility_info, by = 'uID')


## ===============================================
#  Do the deaths by units calculation
## ===============================================
tot_deaths_by_state_year <- read.fst( 'data/cache_data/total_deaths_by_state_year.fst',
                                      as.data.table = T)
setnames( tot_deaths_by_state_year, 'STATE', 'state')

# sum deaths across all US
tot_deaths_us <- tot_deaths_by_state_year[ state %in% state.abb[!(state.abb %in% c('HI', 'AK'))],
                                           .( deaths.fill = sum( deaths.fill),
                                              denom.fill = sum( denom.fill)), by = year]
tot_deaths_us[, state := 'US']
tot_deaths_by_state_year <- rbind( tot_deaths_by_state_year, tot_deaths_us)

# link total deaths and irene's data
irene_adj_pop_fac_tot <- merge( irene_adj_pop_fac, tot_deaths_by_state_year,
                                by = c( 'year', 'state'))

## ====================================================== ##
# load poisson model results
## ====================================================== ##
dir_out <- 'results/'

dat_year_fill <- read.fst( 'data/cache_data/dat_year_fill.fst')
num_uniq_zip <- length( unique( dat_year_fill$zip))

# read in the coefficients
poisson_coefs <- fread( paste0(dir_out, "poisson_model_coefs.csv"))

# read the confidence interval dataset
# loads the object bootstrap_CI
load( paste0(dir_out, "loglinear_coefs_boots.RData"))

# calculate 95% CI
boot_beta_pm <- 
  1.96 * sd( bootstrap_CI$pm25_ensemble) * 
  sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)
boot_beta_hy <- 
  1.96 * sd( bootstrap_CI$hyads) * 
  sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)

# create list of betas
betas_pm <- 
  poisson_coefs[ model == 'pm25_ensemble']$Estimate +
  c( -boot_beta_pm, 0, +boot_beta_pm)
betas_hy <- 
  poisson_coefs[ model == 'hyads']$Estimate +
  c( -boot_beta_hy, 0, +boot_beta_hy)
betas_kr <- log( c( 1.035, 1.056, 1.078)) / 10


## ====================================================== ##
# calculate deaths by unit
## ====================================================== ##
irene_adj_pop_fac_tot[, popwgt_sox := sox.pm25.popwt / pop] #(.99 * 298000000)] #pop]
irene_adj_pop_fac_tot[, `:=` ( deaths_adj_1 = deaths.fill / denom.fill * ( exp( betas_hy[1] * popwgt_sox) - 1) * denom.fill, 
                               deaths_adj_2 = deaths.fill / denom.fill * ( exp( betas_hy[2] * popwgt_sox) - 1) * denom.fill, 
                               deaths_adj_3 = deaths.fill / denom.fill * ( exp( betas_hy[3] * popwgt_sox) - 1) * denom.fill)]

sum( irene_adj_pop_fac_tot[state == 'US' & year == 2006]$deaths_adj_2)
sum( irene_adj_pop_fac_tot[state == 'US' & year == 2011]$deaths_adj_2)

# sum by facility
irene_adj_fac <- irene_adj_pop_fac_tot[, .( deaths_adj_1 = sum( deaths_adj_1),
                                            deaths_adj_2 = sum( deaths_adj_2),
                                            deaths_adj_3 = sum( deaths_adj_3)),
                                       by = .( FacID, year)]

## ===============================================
#  Read in hyads death data for comparison
## ===============================================
# load rdata with tables of deaths attributed to various geographies/facilities
load( 'data/cache_data/annual_deaths_by.RData')

# merge hyads death predictions with adjoint
h_adj <- merge( irene_adj_fac, 
                deaths_by_year_fac_statebinf[model == 'hyads'], 
                by = c( 'FacID', 'year'))

# calculate regressions
h_adj_lm <- h_adj[, .( N = .N,
                       int = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['(Intercept)','Estimate'],
                       int.se = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['(Intercept)','Std. Error'],
                       slo = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['deaths_adj_2','Estimate'],
                       slo.se = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['deaths_adj_2','Std. Error'],
                       r2 = cor( deaths_adj_2, deaths_coef_2) ^ 2,
                       nmb = ( 1 / sum( deaths_adj_2)) * sum( deaths_coef_2 - deaths_adj_2),
                       nme = ( 1 / sum( deaths_adj_2)) * sum( abs( deaths_coef_2 - deaths_adj_2)),
                       rmse = sqrt( ( 1 / .N) * sum( ( ( deaths_coef_2 - mean( deaths_coef_2)) - ( deaths_adj_2 - mean( deaths_adj_2))) ^ 2) ) ),
                  by = .( statebin_facility, year)]
h_adj_us <- h_adj[, .( N = .N,
                       int = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['(Intercept)','Estimate'],
                       int.se = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['(Intercept)','Std. Error'],
                       slo = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['deaths_adj_2','Estimate'],
                       slo.se = coef( summary( lm( deaths_coef_2 ~ deaths_adj_2)))['deaths_adj_2','Std. Error'],
                       r2 = cor( deaths_adj_2, deaths_coef_2)^2,
                       nmb = ( 1 / sum( deaths_adj_2)) * sum( deaths_coef_2 - deaths_adj_2),
                       nme = ( 1 / sum( deaths_adj_2)) * sum( abs( deaths_coef_2 - deaths_adj_2)),
                       rmse = sqrt( ( 1 / .N) * sum( ( ( deaths_coef_2 - mean( deaths_coef_2)) - ( deaths_adj_2 - mean( deaths_adj_2))) ^ 2) ) ),
                  by = .( year)]

h_adj_lm_out <- h_adj_lm[, .( N = N,
                              Intercept = paste( sprintf( '%.1f', int), '±', sprintf( '%.0f', int.se)),
                              Slope = paste( sprintf( '%.2f', slo), '±', sprintf( '%.2f', slo.se)),
                              R2 = sprintf( '%.2f', r2),
                              NMB = scales::percent( nmb),
                              NME = scales::percent( nme),
                              RMSE = sprintf( '%.0f', rmse)),
                         by = .( statebin_facility, year)]
h_adj_lm_out_us <- h_adj_us[, .( N = N,
                              Intercept = paste( sprintf( '%.1f', int), '±', sprintf( '%.0f', int.se)),
                              Slope = paste( sprintf( '%.2f', slo), '±', sprintf( '%.2f', slo.se)),
                              R2 = sprintf( '%.2f', r2),
                              NMB = scales::percent( nmb),
                              NME = scales::percent( nme),
                              RMSE = sprintf( '%.0f', rmse)),
                         by = .( year)]

# get the annual results
sum_cols <- c( 'deaths_adj_1', 'deaths_adj_2', 'deaths_adj_3',
               'deaths_coef_1', 'deaths_coef_2', 'deaths_coef_3')
h_adj[, lapply( .SD, sum),
      .SDcols = sum_cols,
      by = .( year)]

## ===============================================
#  facility deaths comparisons
## ===============================================
geos_chem_eval_plot <- 
  ggplot( h_adj[model == 'hyads'], 
          aes( y = deaths_coef_2, x = deaths_adj_2)) + 
  geom_abline( slope = 1, intercept = 0, size = .25) + 
  geom_errorbar( aes( ymin = deaths_coef_1, ymax = deaths_coef_3), size = 1) +
  geom_errorbarh( aes( xmin = deaths_adj_1, xmax = deaths_adj_3), size = 1) +
  geom_smooth( formula = y ~ x, method = 'lm', se = T, fullrange = T, size = .5) +
  labs( y = expression( HyADS~Coal~PM[2.5]), x = 'GEOS-Chem Adjoint') +
  coord_cartesian( xlim = c( .01, max( h_adj$deaths_coef_3, h_adj$deaths_adj_3)),
                   ylim = c( .01, max( h_adj$deaths_coef_3, h_adj$deaths_adj_3))) +
  # scale_x_log10( ) + 
  # scale_y_log10( ) + 
  facet_grid( year ~ statebin_facility,
              switch = 'y') + 
  theme_bw() + 
  theme( axis.text = element_text( size = 12),
         axis.text.x = element_text( angle = 30),
         axis.title = element_text( size = 18),
         strip.background = element_blank(),
         strip.placement = 'outside',
         strip.text = element_text( size = 18))

ggsave( 'figures/geos_chem_eval_region_deaths.png',
        geos_chem_eval_plot, height = 4, width = 10, 
        unit = 'in', scale = 1.75)

## ===============================================
#  ranks by GEOS-Chem & hyads
## ===============================================
perc.rank <- function(x) trunc(rank(x))/length(x)

# how much deaths associated with 50% emissions?
# merge facility names with facilities dataset
load( 'data/units_coal_1997_2021.rda')
units_all <- units_updated %>% as.data.table
units_all[, FacID := gsub( '-.*$', '', ID) %>% as.integer()]
facs_all <- units_all[, .( SOx = sum( SOx)),
                      by = .( year, FacID)]

deaths_by_fac_statebinf_emiss <- 
  merge( facs_all,
         h_adj[model == 'hyads'],
         by = c( 'FacID', 'year')) 

deaths_by_fac_statebinf_emiss[, `:=` ( SOx_pct = perc.rank( SOx),
                                       deaths_pct_hyads = perc.rank( deaths_coef_2),
                                       deaths_pct_adj = perc.rank( deaths_adj_2))]#,
deaths_by_fac_statebinf_emiss.m <- 
  melt( deaths_by_fac_statebinf_emiss,
        id.vars = c( 'year', 'FacID', 'SOx_pct', 'statebin_facility'),
        measure.vars = c( 'deaths_pct_hyads', 'deaths_pct_adj'))

# by = statebin_facility]
setkey( deaths_by_fac_statebinf_emiss.m, SOx_pct)


# calculate vertical distance from 
deaths_by_fac_statebinf_emiss.m[, `:=` ( min.plot.y = min(value, SOx_pct ),
                                         max.plot.y = max(value, SOx_pct )),
                                by = .( FacID, year, variable, statebin_facility)]

pct_emiss_death.gg <-
  ggplot( deaths_by_fac_statebinf_emiss.m,
          aes( x = SOx_pct, y = value,
               ymin = min.plot.y,
               ymax = max.plot.y,
               color = variable)) + 
  labs( x = expression( SO[2]~emissions~percentile),
        y = 'Excess deaths percentile') +
  geom_abline( slope = 1, intercept = 0,
               color = 'grey80') +
  geom_linerange( alpha = .25) +
  geom_point( alpha = .25, size = 1) + 
  geom_smooth( se = FALSE, alpha = .5) +
  scale_x_continuous( labels = scales::percent_format()) +
  scale_y_continuous( labels = scales::percent_format()) +
  scale_color_brewer( palette = 'Dark2',
                      labels = c( 'deaths_pct_hyads' = expression( HyADS~Coal~PM[2.5]),
                                  'deaths_pct_adj' = expression( GEOS-Chem~Adjoint~PM[2.5]~Sensitivities))) +
  facet_grid( year ~ statebin_facility,
              switch = 'y') + 
  theme_minimal() +
  theme( axis.text = element_text( size = 12),
         axis.text.x = element_text( angle = 30),
         axis.title = element_text( size = 16),
         legend.position = c( .8, -.09),
         legend.direction = 'horizontal',
         legend.title = element_blank(),
         legend.text = element_text( size = 12),
         panel.grid.minor = element_blank(),
         strip.placement = 'outside',
         strip.text = element_text( size = 16))

# save it
ggsave( pct_emiss_death.gg,
        filename = 'figures/deaths_emiss_adj_pct.png',
        height = 7, width = 15, unit = 'in')

## ===============================================
#  evaluate hyads & geos-Chem
## ===============================================









