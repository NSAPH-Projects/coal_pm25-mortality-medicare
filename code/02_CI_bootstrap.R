library( gnm)
library( pbmcapply)
library( dplyr)
library( data.table)
library( fst)

## ======================================================= ##
#  load Xiao's data
## ======================================================= ##
dir_data <- 'data/' 
dir_out <- 'results/' 

aggregate_data <- 
  read.fst( 
    file.path(dir_data, "cache_data", "aggregate_data.fst"),
    as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj'), as.data.table = TRUE)
dat_annual_use <- merge(aggregate_data, dat_annual, by = c("zip", "year")) #, all.x = TRUE)

# clean house
rm( dat_annual, aggregate_data)

## =================================================== ##
# define get_cpu function
## =================================================== ##
get_cpus <- function() {
  ncpus <- system("grep -i ^requestcpus $_CONDOR_JOB_AD", intern = T)
  ncpus <- as.numeric(gsub("\\D+", "", ncpus))
  return(ncpus)
}

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
    zip_sample <- 
      sample( unique( dat_all$zip),
              floor( 2 * sqrt( num_uniq_zip)), replace=T) 
    dat_annual_boots <- dat_all[zip %in% zip_sample]
    
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
    
    # model just based on hyads adjusted to sulfate
    gnm_raw_hy_adj <- 
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
    Poisson_hy.adj <- coef( gnm_raw_hy_adj)
    rm( gnm_raw_hy_adj)
    
    out <- 
      data.table( boots_id = boots_id,
                  pm25_ensemble = Poisson['pm25_ensemble'],
                  hyads = Poisson_hy['Y1'],
                  hyads.adj = Poisson_hy.adj['Y1.adj'])
    
    return( out)    
  }

# run the bootstrap code
# 60GB, 5 Core
bootstrap_CI <- 
  pblapply( # pbmclapply
    1:500,
    bootstrap_CI.fn) %>% rbindlist # mc.cores =  get_cpus()

# save the results
save( bootstrap_CI,
      file=paste0(dir_out, "loglinear_coefs_boots.RData"))


