library( gnm)
library( pbapply)
library( dplyr)
library( data.table)
library( fst)

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/data' 
dir_out <- 'results/' 

aggregate_data <- 
  read.fst( 
    file.path(dir_data, "cache_data", "aggregate_data.fst"),
    as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), as.data.table = TRUE)
dat_annual_use <- merge(aggregate_data, dat_annual, by = c("zip", "year")) #, all.x = TRUE)

# clean house
rm( dat_annual, aggregate_data)

# calculate portion of PM that is not coal PM
dat_annual_use[, pm_not_coalPM := pm25_ensemble - Y1]
## =================================================== ##
# Bootstrap on zip-code cluster to obtain robust CIs account for spatial correlation
## =================================================== ##

# split the data
# dat_annual_use.list <- split( dat_annual_use, list(dat_annual_use$zip))
num_uniq_zip <- length( unique( dat_annual_use$zip))

# bootstrap function to fun in parallel
bootstrap_CI.fn <- 
  function( boots_id,
            dat_all = dat_annual_use,
            num_zips = num_uniq_zip){
    set.seed(boots_id)
    # print(boots_id)
    zip_sample <- 
      sample( unique( dat_all$zip),
              floor( 2 * sqrt( num_uniq_zip)), replace=T) 
    dat_annual_boots <- dat_all[zip %in% zip_sample]
    
    ## =========================================
    # models not adjusted for PM
    ## =========================================
    gnm_raw <-
      gnm( dead ~ pm25_ensemble +
             mean_bmi + smoke_rate + hispanic + pct_blk +
             medhouseholdincome + medianhousevalue +
             poverty + education + popdensity + pct_owner_occ +
             summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
             as.factor(year) + as.factor(region) +
             offset(log(time_count)),
           eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
           data = dat_annual_boots,
           family = poisson(link = "log"))
    Poisson <- coef(gnm_raw)
    rm( gnm_raw)
    
    gnm_raw_hy <- 
      gnm( dead ~ Y1 +
             mean_bmi + smoke_rate + hispanic + pct_blk +
             medhouseholdincome + medianhousevalue +
             poverty + education + popdensity + pct_owner_occ +
             summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
             as.factor(year) + as.factor(region) +
             offset(log(time_count)),
           eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
           data = dat_annual_boots,
           family = poisson(link = "log"))
    Poisson_hy <- coef(gnm_raw_hy)
    rm( gnm_raw_hy)
    
    ## HyADS model, 2003 and earlier
    gnm_raw_hy_early <- 
      gnm(dead ~ Y1 + 
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots[ year <= 2003],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_early <- coef(gnm_raw_hy_early)
    rm( gnm_raw_hy_early)
    
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
          data = dat_annual_boots[ year >= 2004 & year <= 2007],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_mid <- coef(gnm_raw_hy_mid)
    rm( gnm_raw_hy_mid)
    
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
          data = dat_annual_boots[ year >= 2008],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_late <- coef(gnm_raw_hy_late)
    rm( gnm_raw_hy_late)
    
    # model just based on hyads adjusted to sulfate
    gnm_raw_hy.adj <- 
      gnm( dead ~ Y1.adj +
             mean_bmi + smoke_rate + hispanic + pct_blk +
             medhouseholdincome + medianhousevalue +
             poverty + education + popdensity + pct_owner_occ +
             summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
             as.factor(year) +
             as.factor(region) +
             offset(log(time_count)),
           eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
           data = dat_annual_boots,
           family = poisson(link = "log"))
    Poisson_hy.adj <- coef( gnm_raw_hy.adj)
    rm( gnm_raw_hy.adj)
    
    ## ====================================================
    #  adjusted models - total pm
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
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_pm <- coef(gnm_raw_hy_pm)
    rm( gnm_raw_hy_pm)
    
    # model just based on hyads adjusted to sulfate with pm adjustment
    gnm_raw_hy.adj_pm <- 
      gnm(dead ~ Y1.adj + pm25_ensemble +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) +
            as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy.adj_pm <- coef( gnm_raw_hy.adj_pm)
    rm( gnm_raw_hy.adj_pm)
    
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
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_dpm <- coef(gnm_raw_hy_dpm)
    rm( gnm_raw_hy_dpm)
    
    ## HyADS model adjusted by [total PM - coal PM], 2003 and earlier
    gnm_raw_hy_dpm_early <- 
      gnm(dead ~ Y1 + pm_not_coalPM +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots[ year <= 2003],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_dpm_early <- coef(gnm_raw_hy_dpm_early)
    rm( gnm_raw_hy_dpm_early)
    
    ## HyADS model adjusted by [total PM - coal PM], 2004-2007
    gnm_raw_hy_dpm_mid <- 
      gnm(dead ~ Y1 + pm_not_coalPM +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots[ year >= 2004 & year <= 2007],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_dpm_mid <- coef(gnm_raw_hy_dpm_mid)
    rm( gnm_raw_hy_dpm_mid)
    
    ## HyADS model adjusted by [total PM - coal PM], 2008 and later
    gnm_raw_hy_dpm_late <- 
      gnm(dead ~ Y1 + pm_not_coalPM +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots[ year >= 2008],
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_dpm_late <- coef(gnm_raw_hy_dpm_late)
    rm( gnm_raw_hy_dpm_late)
    
    # model just based on hyads adjusted to sulfate with pm adjustment
    gnm_raw_hy.adj_dpm <- 
      gnm(dead ~ Y1.adj + pm_not_coalPM +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) +
            as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy.adj_dpm <- coef( gnm_raw_hy.adj_dpm)
    rm( gnm_raw_hy.adj_dpm)
    
    ## HyADS model adjusted by [total PM - coal PM] and other pollutants
    gnm_raw_hy_pollutants <- 
      gnm(dead ~ Y1 + pm_not_coalPM + ozone.current_year +
            no2.current_year + ozone_summer.current_year +
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_pollutants <- coef(gnm_raw_hy_pollutants)
    rm( gnm_raw_hy_pollutants)
    
    ## HyADS model adjusted by [total PM - coal PM] and other pollutants
    gnm_raw_hy_pollutants_no_o3 <- 
      gnm(dead ~ Y1 + pm_not_coalPM + no2.current_year + 
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_pollutants_no_o3 <- coef(gnm_raw_hy_pollutants_no_o3)
    rm( gnm_raw_hy_pollutants_no_o3)

    ## raw HyADS model adjusted by [total PM - coal PM] and other pollutants
    gnm_raw_hy_raw_pollutants_no_o3 <- 
      gnm(dead ~ Y1_raw + pm_not_coalPM + no2.current_year + 
            mean_bmi + smoke_rate + hispanic + pct_blk +
            medhouseholdincome + medianhousevalue +
            poverty + education + popdensity + pct_owner_occ +
            summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
            as.factor(year) + as.factor(region) +
            offset(log(time_count)),
          eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
          data = dat_annual_boots,
          family = poisson(link = "log"))
    
    # save only the summaries, clear large object
    Poisson_hy_raw_pollutants_no_o3 <- coef(gnm_raw_hy_raw_pollutants_no_o3)
    rm( gnm_raw_hy_raw_pollutants_no_o3)
    gc()
    
    # create output data table
    out <- 
      data.table( boots_id = boots_id,
                  pm25_ensemble =      Poisson['pm25_ensemble'],
                  hyads =              Poisson_hy['Y1'],
                  hyads_early =        Poisson_hy_early['Y1'],
                  hyads_mid =          Poisson_hy_mid['Y1'],
                  hyads_late =         Poisson_hy_late['Y1'],
                  hyads_pm =           Poisson_hy_pm['Y1'],
                  hyads_dpm =          Poisson_hy_dpm['Y1'],
                  hyads_dpm_early =    Poisson_hy_dpm_early['Y1'],
                  hyads_dpm_mid =      Poisson_hy_dpm_mid['Y1'],
                  hyads_dpm_late =     Poisson_hy_dpm_late['Y1'],
                  hyads_pollutants =   Poisson_hy_pollutants['Y1'],
                  hyads_pollutants_no_o3=unname( Poisson_hy_pollutants_no_o3['Y1']),
                  hyads_raw_pollutants_no_o3=unname( Poisson_hy_raw_pollutants_no_o3['Y1_raw']),
                  hyads.adj =          Poisson_hy.adj['Y1.adj'],
                  hyads.adj_pm =       Poisson_hy.adj_pm['Y1.adj'],
                  hyads.adj_dpm =      Poisson_hy.adj_dpm['Y1.adj'])
    
    gc()
    return( out)    
  }

# run the bootstrap code
# 150 GB RAM
bootstrap_CI <- 
  pblapply( # pbmclapply
    1:500,
    bootstrap_CI.fn) %>% rbindlist # mc.cores =  get_cpus()

# save the results
save( bootstrap_CI,
      file=paste0(dir_out, "loglinear_coefs_boots.RData"))


