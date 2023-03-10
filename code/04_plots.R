library( fst)
library( data.table)
library( sf)
library( ggplot2)
library( magrittr)
library( viridis)
library( ggrepel)
library( ggpattern)

dir_data <- 'data/data' 
dir_out <- 'results/'

## ================================================== ##
#  load and manipulate deaths data
## ================================================== ##
deaths_by_year_unit.hy <- read.fst( paste0(dir_out, "deaths_by_year_hyRR.fst"), as.data.table = TRUE)
deaths_by_year_unit.pm <- read.fst( paste0(dir_out, "deaths_by_year_pmRR.fst"), as.data.table = TRUE)
deaths_by_year_unit.kr <- read.fst( paste0(dir_out, "deaths_by_year_krRR.fst"), as.data.table = TRUE)

# combine into single data.table
deaths_by_year_unit.hy[, model := 'hyads']
deaths_by_year_unit.pm[, model := 'pm25_ensemble']
deaths_by_year_unit.kr[, model := 'pm25_krewski']
deaths_by_year_unit <- 
  rbind( deaths_by_year_unit.hy, 
         deaths_by_year_unit.pm,
         deaths_by_year_unit.kr)

# read in deaths from total hyads exposure
deaths_by_state_year_all <- read.fst(paste0(dir_out, "deaths_by_year_total_coal_pm25.fst"), as.data.table = TRUE)
deaths_by_state_year_all_pm25 <- read.fst(paste0(dir_out, "deaths_by_year_total_pm25_ensemble.fst"), as.data.table = TRUE)

## ======================================================= ##
# load hyads data
## ======================================================= ##
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj', 'Y1_raw'), as.data.table = TRUE)

## ================================================== ##
# give states a bin ##REMINDER--check Washington, DC!
## ================================================== ##
state_bins <- data.table( state = state.name,
                          STATE = state.abb)
state_bins[ state %in% c( 'Mississippi', 'Georgia', 'Alabama', 
                          'Florida', 'South Carolina', 
                          'Tennessee', 'North Carolina'),
            statebin_facility := 'Southeast']
state_bins[ state %in% c( 'West Virginia', 'Delaware', 'Kentucky',
                          'Maryland', 'Virginia', 'New Jersey', 
                          'Pennsylvania'),
            statebin_facility := 'Mid-Atlantic']
state_bins[ state %in% c( 'Arkansas', 'Louisiana', 'Oklahoma', 'Texas'),
            statebin_facility := 'Southcentral']
state_bins[ state %in% c( 'New York', 'New Hampshire', 'Massachusetts', 'Maine',
                          'Connecticut', 'Vermont', 'Rhode Island'),
            statebin_facility := 'Northeast']
state_bins[ state %in% c( 'Nevada', 'Utah', 'Arizona', 'New Mexico',
                          'Colorado', 'Wyoming', 'California', 'Idaho',
                          'Washington', 'Oregon', 'Montana'),
            statebin_facility := 'West']
state_bins[ state %in% c( 'Wisconsin','Michigan', 'Illinois',
                          'Indiana', 'Ohio'),
            statebin_facility := 'Midwest east']
state_bins[ state %in% c( 'North Dakota', 'South Dakota', 
                          'Minnesota', 'Kansas', 'Iowa', 
                          'Nebraska', 'Missouri'),
            statebin_facility := 'Midwest west']

# set up a comperable dataset for the receptors
state_bins_zips <- state_bins[, .(STATE, statebin_facility)]
setnames( state_bins_zips, 'statebin_facility', 'statebin_zip')

# add region for DC
state_bins_zips <- rbind( state_bins_zips, 
                          data.table( STATE = 'DC', 
                                      statebin_zip = 'Mid-Atlantic'))

# clean up
state_bins[ ,STATE := NULL]

## ================================================== ##
# merge with unit information
## ================================================== ##
# read in facility info
info.list <- c( 'Facility.Name', 'Facility.ID..ORISPL.', 'Unit.ID', 'State',
                'Facility.Latitude',  'Facility.Longitude')
facility_info <- fread( file.path( dir_data, 'AMPD_Unit.csv'),
                        select = info.list) %>% unique
setnames( facility_info, 
          c( 'Facility.ID..ORISPL.', 'Unit.ID'), 
          c( 'FacID', 'unitID'))

# tricky text problems
facility_info[, `:=` ( state = state.name[ match( State, state.abb)],
                       uID = gsub( '-|\\#|\\*', '.', paste( FacID, unitID, sep = '.')))]
facility_info <- facility_info[!(duplicated( uID))]

# one facility giving name issues
facility_info[FacID == '1004', Facility.Name := 'Edwardsport']

# group some of the states
unit.states <- merge( facility_info, state_bins, by = 'state')[!duplicated( uID)]
# setnames( unit.states, 'statebin', 'statebin_unit')

# remove X at start of variable to give unit names
deaths_by_year_unit[, uID := gsub( '^X', '', variable)]

# merge with emissions and states
deaths_and_units <- merge( deaths_by_year_unit, unit.states, by = c( 'uID')) #, all.x = TRUE)
deaths_and_units <- merge( deaths_and_units, state_bins_zips, by = c( 'STATE'))

# clarify state names
setnames( deaths_and_units,
          c( 'STATE', 'State'),
          c( 'state_zip', 'state_facility'))

## ================================================== ##
# sum deaths by various groups
## ================================================== ##
# define columns we will sum for
sum_cols <- c( 'deaths_coef_1', 'deaths_coef_2', 'deaths_coef_3')

# deaths by unit-year-statebin_facility 
deaths_by_unit_year_statebin <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, uID, year, statebin_facility)]

# deaths by unit-year-statebin_facility 
deaths_by_fac_year_statebinz <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, year, Facility.Name, 
                           FacID, state_facility, statebin_facility, 
                           Facility.Latitude,  Facility.Longitude, state_zip)]
deaths_by_fac_year_statebinz.pltfrm <- 
  deaths_by_fac_year_statebinz[model == 'hyads'][,model := NULL]
fwrite( deaths_by_fac_year_statebinz.pltfrm,
        file.path(dir_data, "platform_data", 'pm25_facility_state.csv'))

# deaths by year
deaths_by_year <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, year)]
deaths_by_year.present <- 
  deaths_by_year[, .( deaths = 
                        paste( 
                          format(signif( deaths_coef_2, 2), big.mark=","),
                          ' (',  
                          format(signif( deaths_coef_1, 2), big.mark=","),
                          '-', 
                          format(signif( deaths_coef_3, 2), big.mark=","), ')',
                          sep = ''
                        )),
                 by = .( model, year)] %>%
  dcast( year ~ model, 
         value.var = 'deaths',
         variable.name = '')
fwrite( deaths_by_year.present,
        file.path(dir_data, "platform_data", 'deaths_by_year.csv'))

# deaths by year calculated with all HyADS
deaths_by_year_all <- 
  deaths_by_state_year_all[, lapply( .SD, sum),
                           .SDcols = sum_cols,
                           by = .( model, year)]

# deaths by year calculated with all pm25
deaths_by_year_all.pm25 <- 
  deaths_by_state_year_all_pm25[, lapply( .SD, sum),
                                .SDcols = sum_cols,
                                by = .( year)]

# total deaths
deaths_total <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model)]

# total deaths calculated with all HyADS
deaths_total_all <- 
  deaths_by_state_year_all[, lapply( .SD, sum),
                           .SDcols = sum_cols,
                           by = .( model)]

# death by statebin_zip and statebin_facility
deaths_by_statebinf_statebinz <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, statebin_facility, statebin_zip)]

# deaths by facility-year-statebin
deaths_by_year_fac_statebinf <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, year, FacID, statebin_facility)]

# deaths by facility and facility statebin
deaths_by_fac_statebinf <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, FacID, statebin_facility)]

# deaths by statebin
deaths_by_statebinf <- 
  deaths_and_units[, lapply( .SD, sum),
                   .SDcols = sum_cols,
                   by = .( model, statebin_facility)]

## ================================================== ##
# plot total deaths by year
## ================================================== ##
# create labels for each model type
modnames <- c(   'hyads' = expression( Coal~PM[2.5]~(this~study)),
                 'pm25_ensemble'  = expression( PM[2.5]~(Wu~et~al.*','~2020)),
                 'pm25_krewski' = expression( PM[2.5]~(Krewski~et~al.*','~2008)))

deaths_year.gg <-
  ggplot( deaths_by_year,
          aes( x = year,
               # fill = model
          )) + 
  geom_rect( xmin = 2017, xmax = 2020,
             ymin = 0, ymax = 10, fill = 'grey95', color = NA) +
  geom_ribbon_pattern( aes( ymax = deaths_coef_3 / 10000,
                            ymin = deaths_coef_1 / 10000,
                            pattern_angle = model),
                       size = .2,
                       fill = NA,
                       color = 'black',
                       pattern_spacing = .025,
                       # pattern_angle = 0,
                       pattern_color = NA,
                       pattern_fill = 'grey50',
                       alpha = .8) + 
  scale_pattern_angle_discrete( labels = modnames,
                                name = 'Risk Ratio source') + 
  labs( y = 'Annual Excess Deaths (10,000)') +
  # coord_cartesian( ylim = c( 0, NULL)) +
  theme_bw() + 
  theme( axis.text = element_text( size = 12),
         axis.title.y = element_text( size = 16),
         axis.title.x = element_blank(),
         legend.position = c(.65,.8),
         legend.text = element_text( size = 12),
         legend.text.align = 0,
         legend.title = element_text( size = 14),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))
deaths_year.gg

ggsave( deaths_year.gg,
        filename = 'figures/deaths_per_year.png',
        height = 4, width = 8, unit = 'in')

## ================================================== ##
# plot total deaths by year sensitivities
## ================================================== ##
# create labels for each model type
modnames.s <- c( 'Constant deaths & death rates' = expression( Change~from~coal~PM[2.5]),
                 'Constant coal PM2.5'  = expression( Change~from~deaths~and~death~rates))

deaths_by_year_sensitivities.c <- fread( paste0(dir_out, "deaths_by_year_sensitivities.csv"))


deaths_year.gg.sens <-
  ggplot( deaths_by_year_sensitivities.c,
          aes( x = year,
               linetype = case)) + 
  geom_rect( xmin = 2017, xmax = 2020,
             ymin = -10, ymax = 10, fill = 'grey95', color = NA) +
  geom_hline( yintercept = 0) +
  geom_line( aes( y = deaths_coef_2_diff),
             linewidth = 1.25,
             alpha = .8) + 
  scale_y_continuous( labels = scales::percent) +
  scale_linetype_manual( labels = modnames.s,
                         values = c( 'solid', 'dotted')) + 
  labs( y = 'Annual Excess Deaths Change from 1999 (%)') +
  # coord_cartesian( ylim = c( 0, NULL)) +
  theme_bw() + 
  theme( axis.text = element_text( size = 12),
         axis.title.y = element_text( size = 12),
         axis.title.x = element_blank(),
         legend.key.width = unit( 1, 'in'),
         legend.position = c(.3,.13),
         legend.text = element_text( size = 12),
         legend.text.align = 0,
         legend.title = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))
deaths_year.gg.sens

ggsave( deaths_year.gg.sens,
        filename = 'figures/deaths_per_year_sensitivity.png',
        height = 4, width = 8, unit = 'in')

## ================================================== ##
# sum deaths by various groups
## ================================================== ##
# add the number of deaths to each statebin name
deaths_by_statebinf[, statebin_lab := 
                      eval( expression( paste0( statebin_facility, '\n',
                                                format( round( deaths_coef_2, -2), big.mark = ','), 
                                                '\n(', 
                                                gsub( ' ', '', format( round( deaths_coef_1, -2), big.mark = ',')), 
                                                '-', 
                                                gsub( ' ', '', format( round( deaths_coef_3, -2), big.mark = ',')), ')')))]

# merge with statebin totals
deaths_by_fac_year_statebin_lab <- 
  merge( deaths_by_year_fac_statebinf, 
         deaths_by_statebinf[, .( statebin_facility, statebin_lab, model)],
         by = c( 'statebin_facility', 'model'))

# do some ordering
setorder( deaths_by_statebinf, -deaths_coef_2)
deaths_by_statebinf[, statebin_lab := factor( statebin_lab,
                                              levels = c( statebin_lab))]

## ================================================== ##
# facility names => update!
## ================================================== ##

# merge facility names with facilities dataset
facility_info_noyear <- unique( facility_info[, .( Facility.Name, FacID, State, state)])
deaths_by_fac_info <- merge( deaths_by_fac_statebinf, facility_info_noyear, by = 'FacID')
deaths_by_fac_info <- merge( deaths_by_fac_info, deaths_by_statebinf[, .( model, statebin_facility, statebin_lab)],
                             by = c( 'statebin_facility', 'model'))
deaths_by_fac_info <- deaths_by_fac_info[model == 'hyads']

# number of labels per region
nlabs_region <- data.table( statebin_facility = c( "Mid-Atlantic", "Midwest east",
                                                   "Southeast", "Southcentral", 
                                                   "Midwest west", "West", "Northeast"),
                            n = c( 6, 7, 5, 3, 2, 1, 1))
deaths_by_fac_info <- merge( deaths_by_fac_info, nlabs_region, by = 'statebin_facility')
setkey( deaths_by_fac_info, deaths_coef_2)

# create labels for  most-death facilities & remove "Power station"
deaths_by_fac_info.lab <- deaths_by_fac_info[ , tail( .SD, unique( n)), by = .( statebin_facility)]
deaths_by_fac_info.lab[, lab := eval( expression( paste0( Facility.Name, ', ', State)))]#, '\n',
# format( round( deaths_coef_2), big.mark = ','), 
# ' (', 
# gsub( ' ', '', format( round( deaths_coef_1), big.mark = ',')), 
# '-', 
# gsub( ' ', '', format( round( deaths_coef_3), big.mark = ',')), ')')))]
deaths_by_fac_info.lab[, lab := gsub( ' Power Station', '', lab)]

# five year bins
deaths_by_fac_year_statebin_lab[, yearbin := cut( year, c( 1999, 2005, 2011, 2017, 2020), 
                                                  include.lowest = T, right = F,
                                                  ordered_result = T,
                                                  labels = c( '1999-2004', '2005-2010', '2011-2016', '2017-2020')
)]

# group some states together
deaths_by_fac_year_statebin_lab[, lab := '']
setorder( deaths_by_fac_statebinf, -deaths_coef_2)
deaths_by_fac_year_statebin_lab[, statebin_lab := factor( statebin_lab, levels = levels( deaths_by_statebinf$statebin_lab))]
deaths_by_fac_year_statebin_lab[, FacID := factor( FacID, levels = rev( unique( deaths_by_fac_statebinf$FacID)))]
deaths_by_fac_info.lab[, FacID := factor( FacID, levels = rev( unique( deaths_by_fac_statebinf$FacID)))]
deaths_by_fac_info.lab[, statebin_lab := factor( statebin_lab, levels = levels( deaths_by_statebinf$statebin_lab))]

## ================================================== ##
# save some of the data 
## ================================================== ##
save( deaths_by_fac_year_statebin_lab, deaths_by_fac_info.lab,
      deaths_by_fac_info, deaths_by_unit_year_statebin, 
      deaths_by_year, deaths_by_year_fac_statebinf, 
      deaths_by_fac_statebinf, deaths_by_statebinf, state_bins,
      file = 'data/data/cache_data/annual_deaths_by.RData')

## ================================================== ##
# Medicare deaths/death rates 
## ================================================== ##
# load in all deaths (including projected) by year
sum_deaths_year_state <- read.fst( 'data/data/cache_data/total_deaths_by_state_year.fst', as.data.table = TRUE)
sum_deaths_year <- sum_deaths_year_state[,.( deaths.fill = sum( deaths.fill),
                                             denom.fill = sum( denom.fill)),
                                         by = year]
sum_deaths_year[, dates := as.Date( paste( year, '01', '01', sep = '-'))]

# seperate into observed and projected
sum_deaths_year.pre <- sum_deaths_year[year <= 2016]
sum_deaths_year.post <- sum_deaths_year[year >= 2016]

ggdeaths <- ggplot( mapping = aes( x = dates, y = deaths.fill / 1e6)) + 
  geom_line( data = sum_deaths_year.pre) + 
  geom_line( data = sum_deaths_year.post, linetype = 'dashed') + 
  labs( y = 'Deaths (millions)') + 
  scale_y_continuous( label = scales::comma) + 
  expand_limits( y = 0) + 
  theme_bw() + 
  theme( axis.text = element_text( size = 18),
         axis.title = element_text( size = 20),
         axis.title.x = element_blank())
ggdeathrate <- ggplot( mapping = aes( x = dates, 
                                      y = deaths.fill / denom.fill * 1e3)) + 
  geom_line( data = sum_deaths_year.pre) + 
  geom_line( data = sum_deaths_year.post, linetype = 'dashed') + 
  labs( y = 'Death rate per 1,000') + 
  scale_y_continuous( label = scales::comma) + 
  expand_limits( y = 0) + 
  theme_bw() + 
  theme( axis.text = element_text( size = 18),
         axis.title = element_text( size = 20),
         axis.title.x = element_blank())

gg_combine <- 
  cowplot::plot_grid( ggdeaths, 
                      ggdeathrate,
                      labels = NULL, ncol = 2, 
                      align = 'v', axis = 'lr') 

ggsave( gg_combine, file = 'figures/annual_deaths_combine.png', width = 11, height = 4, scale = 1)


## ================================================== ##
# plot it! 
## ================================================== ##
gg_bardeaths <-
  ggplot( deaths_by_fac_year_statebin_lab[model == 'hyads'], #[sample(1:7000, 500)],
          aes( y = deaths_coef_2, x = FacID, label = lab)) + 
  labs( y = 'Excess Deaths') +
  geom_col( aes( fill = yearbin), 
            color = NA, width = 1, 
            position = position_stack( reverse = T)) + 
  scale_fill_viridis( option = 'B', 
                      discrete = T,
                      direction = -1,
                      begin = .2,
                      end = .9,
                      guide = guide_legend( ncol = 1)) +
  geom_label_repel( data = deaths_by_fac_info.lab,#[statebin == 'Mid-Atlantic'],
                    size = 7,
                    seed = 123,
                    force = 1,
                    # angle = 90,
                    direction  = "x",
                    label.padding = unit(0.3, "lines"),
                    vjust = .5,
                    xlim = c( NA, 15),
                    ylim = c( 500, 11000),
                    segment.size = 0.2) +
  facet_grid( statebin_lab ~ ., switch = 'y',
              scales = 'free_y', space = 'free') +
  scale_y_continuous(expand = c(0,0),
                     breaks = seq( 0, 20000, 2500)) +
  coord_flip( ylim = c( 0, 11000), clip = 'on') + #
  # coord_cartesian( ) +   # This keeps the labels from disappearing
  theme_bw() + 
  theme( axis.text = element_text( size = 22),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         axis.title = element_text( size = 22),
         axis.title.y = element_blank(),
         legend.direction = 'horizontal',
         legend.key.size = unit( 2, 'lines'),
         legend.position = c( .9, .45),
         legend.text = element_text( size = 20),
         legend.title = element_blank(),
         panel.background = element_rect( fill = NA),
         panel.border = element_blank(),
         # panel.grid.major.x = element_rect( color = 'grey90'),
         panel.grid.major.y = element_blank(),
         # plot.margin = unit(c(.1,2,.1,.1), "lines"),
         panel.spacing = unit( 0.1, 'cm'),
         plot.background = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 16, 
                                    margin = margin( 0, 0, 0, 0)),
         strip.text.y.left = element_text( angle = 0))

# gg_bardeaths
ggsave( gg_bardeaths,
        filename = 'figures/deaths_coal_pm25_unit_year.png',
        height = 15, width = 20, unit = 'in')
ggsave( gg_bardeaths,
        filename = 'figures/deaths_coal_pm25_unit_year.pdf',
        height = 15, width = 20, unit = 'in', device = cairo_pdf)

## ================================================== ##
# region plot with emissions
## ================================================== ##
# get units dataset
load( 'data/data/units_coal_1997_2021.rda')
units_all <- units_updated %>% as.data.table
units_all[, FacID := gsub( '-.*$', '', ID)]

# get USA dataset -- just lower 48 states
states48 <- state.name[!(state.name %in% c( 'Alaska', 'Hawaii'))]
usa <- rnaturalearth::ne_states( 
  country = "United States of America",
  returnclass = 'sf')
usa <- usa[ usa$name %in% states48,] #%>% as.data.table

usa_region <- merge( usa, state_bins, by.x = 'name', by.y = 'state')
usa_region.sf <- lapply( c( na.omit( unique( state_bins$statebin_facility))),
                         function( reg){
                           sfin <- usa_region[usa_region$statebin_facility == reg,'statebin_facility']
                           sfout <- st_union( sfin) %>% st_sf %>% st_cast( 'MULTIPOLYGON')
                           sfout$statebin_facility <- reg
                           
                           # class( sfout$x) <- c( "sfc_MULTIPOLYGON", 'sfc')
                           # st_crs( sfout$x) <- 4326#st_crs( usa)
                           return( sfout)
                         }) %>% rbindlist #%>% st_as_sf( sf_column_name = 'geometry') #dplyr::bind_rows()

# create labels
usa_region.sf$xlab <- c( -76, -118, -105, -75, -70, -83, -97)
usa_region.sf$ylab <- c(  32,   50,   28,  46,  38,  49,  50)

# calculate sox in million kg
units_all[, SOx_kg := SOx * 907.185 / 1e6]

# create spaital object
units_all.sf <- st_as_sf(units_all, coords = c( 'Longitude', 'Latitude'),
                         crs = st_crs( 'WGS84'))

# sum units SO2 emissions by facility
map_emissions <- 
  ggplot( ) +
  geom_sf( data = usa,
           # aes( geometry = 'geometry'),
           color = 'grey50',
           fill = 'white', size = .1) +
  geom_sf( data = usa_region.sf,
           fill = NA,
           aes( geometry = geometry),
           # aes( fill = statebin),
           color = 'black', size = .8) +
  geom_label( data = usa_region.sf,
              size = 6.5,
              aes( label = statebin_facility,
                   x = xlab, y = ylab)) +
  stat_summary_hex( aes( x = Longitude, y = Latitude, z = SOx_kg),
                    data = units_all[ Latitude > 25 & Latitude < 60,], 
                    fun = 'sum',
                    inherit.aes = F, bins = 40,
                    drop = F,
                    color = 'black', alpha = .7, size = .1
  ) +
  scale_fill_viridis( name = expression(paste(SO[2], ' emissions [million kg]')),
                      option = 'D',
                      direction = -1,
                      begin = .25,
                      limits = c( 0, 1000),
                      oob = scales::squish,
                      guide = guide_colorbar( title.position='bottom',
                                              title.hjust = 0.5,
                                              title.vjust = 0)) +
  theme_minimal() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         legend.direction = 'horizontal',
         legend.key.width = unit( .32, 'inches'),
         legend.position = c( .17, .15),
         legend.text = element_text( size = 14),#, angle = 30),
         legend.title = element_text( size = 16),
         panel.grid = element_blank(),
         plot.background = element_rect( color = 'black', fill = 'white',
                                         size = .2))
map_emissions

ggsave( 'figures/emissions_hex_plot.png',
        map_emissions,
        height = 2.9, width = 5, unit = 'in', scale = 1.6)
ggsave( 'figures/emissions_hex_plot.pdf',
        map_emissions, device = cairo_pdf,
        height = 2.9, width = 5, unit = 'in', scale = 1.6)

## ================================================== ##
# plot specific facilities change over time
## ================================================== ##
nlabs_region2 <- data.table( statebin_facility = c( "Mid-Atlantic", "Midwest east",
                                                    "Southeast", "Southcentral", 
                                                    "Midwest west", "West", "Northeast"),
                             n = c( 2, 2, 2, 2, 2, 2, 2))
deaths_by_fac_info2 <- merge( deaths_by_fac_info, nlabs_region2, by = 'statebin_facility')
setkey( deaths_by_fac_info2, deaths_coef_2)

# select facilities and rank them
deaths_by_fac_info.lab <- deaths_by_fac_info2[ , tail( .SD, unique( n.y)), by = .( statebin_facility)]
deaths_by_fac_info.lab[, deaths_rank := frank( deaths_coef_2), by = statebin_facility]

largest_facilities <- 
  deaths_by_year_fac_statebinf[ FacID %in% deaths_by_fac_info.lab$FacID &
                                  model == 'hyads'] %>%
  merge( deaths_by_fac_info.lab[, .( FacID, deaths_rank)], by = 'FacID')
largest_facilities[, lab := '']

# create labels
deaths_by_fac.lab1 <- 
  deaths_by_fac_info.lab[ FacID %in% largest_facilities$FacID,
                          .( FacID, Facility.Name, State)]

# merge with 1999 deaths
deaths_by_fac.lab <- 
  merge( deaths_by_fac.lab1, 
         largest_facilities[year == 1999],
         by = c( 'FacID'))
deaths_by_fac.lab[, lab := paste0( Facility.Name, ', ', State)]
deaths_by_fac.lab[, lab := gsub( ' Power Station| Power', '', lab)]

# look at keystone
deaths_by_year_fac_statebinf[year %in% 1999:2008 & 
                               FacID == 3136 & model == 'hyads',
                             lapply( .SD, function( f) sum( f) / .N),
                             .SDcols = sum_cols]
deaths_by_year_fac_statebinf[year %in% 2011:2020 & 
                               FacID == 3136 & model == 'hyads',
                             lapply( .SD, function( f) sum( f) / .N),
                             .SDcols = sum_cols]

change_over_time.gg <- 
  ggplot( largest_facilities,
          aes( x = year, y = deaths_coef_2, label = lab,
               color = factor( deaths_rank), group = as.factor( FacID))) + 
  labs( y = 'Annual Excess Deaths') +
  geom_line( aes( x = year, y = deaths_coef_2, group = FacID),
             color = 'grey20', alpha = 1,
             linewidth = .11,
             inherit.aes = FALSE,
             data = deaths_by_fac_year_statebin_lab[model == 'hyads']) +
  geom_line( size = 1) + 
  geom_text_repel( data = deaths_by_fac.lab,
                   direction = 'y',
                   xlim = c( NA, 1999),
                   hjust = 0,
                   size = 4,
                   x = 1999) +
  scale_color_viridis( discrete = TRUE,
                       begin = .4,
                       end = .8) + 
  facet_wrap(  . ~ statebin_facility,
               nrow = 2) + 
  scale_x_continuous( limits = c( 1980, NA),
                      breaks = seq( 2000, 2020, 10),
                      minor_breaks = seq( 2000, 2020, 5)) +
  theme_minimal() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_text( size = 14),
         axis.title.x = element_blank(),
         legend.position = 'none',
         panel.grid.major = element_line( color = 'grey90', size = .2),
         panel.grid.minor = element_blank(),
         strip.text = element_text( size = 14))
change_over_time.gg

# save it
ggsave( change_over_time.gg,
        filename = 'figures/deaths_coal_pm25_2facilities_year.pdf',
        height = 5, width = 15, unit = 'in', device = cairo_pdf)
ggsave( change_over_time.gg,
        filename = 'figures/deaths_coal_pm25_2facilities_year.png',
        height = 5, width = 15, unit = 'in')

## ================================================== ##
# percent of deaths attributable to total PM
## ================================================== ##
deaths_by_year_merge <- 
  merge( deaths_by_year[model == 'hyads', .( year, deaths_coef_2)],
         deaths_by_year_all.pm25[ , .( year, deaths_coef_2)], by = 'year')
setnames( deaths_by_year_merge, 
          c( 'deaths_coef_2.x', 'deaths_coef_2.y'),
          c( 'deaths_hyads', 'deaths_pm'))
deaths_by_year_merge[, `:=` (diff_coal = deaths_pm - deaths_hyads,
                             fraction_coal = deaths_hyads / deaths_pm)]

# melt to get 
deaths_by_year_merge.m <- 
  melt( deaths_by_year_merge[,.( year, deaths_hyads, deaths_pm)], id.vars = 'year')

# data.table for names
names.dt <- 
  c(   'deaths_hyads' = expression( Coal~PM[2.5]),
       'deaths_pm'  = expression( Total~PM[2.5]),
       'pm25_krewski' = expression( PM[2.5]~(Krewski~et~al*','~2008)))

deaths_coal_pm.gg <-
  ggplot(deaths_by_year_merge.m[year %in% 2000:2016],
         aes( x = year)) + 
  geom_col( aes( fill = variable,
                 y = value / 10000), 
            # alpha = .5,
            position = position_dodge2( padding = 0, 
                                        width = .7,
                                        reverse = TRUE)) +
  scale_fill_brewer( name = 'Polutant',
                     palette = 'Dark2',
                     labels = names.dt) + 
  geom_text( data = deaths_by_year_merge[year %in% 2000:2016],
             aes( label = scales::percent( fraction_coal,
                                           accuracy = 1)),
             y = -.5) +
  guides(fill = guide_legend( reverse = TRUE)) +
  labs( y = 'Annual Excess Deaths (10,000)') +
  expand_limits( y = -.75) +
  theme_bw() + 
  theme( axis.text = element_text( size = 12),
         axis.title.y = element_text( size = 16),
         axis.title.x = element_blank(),
         legend.position = c(.8,.8),
         legend.text = element_text( size = 12),
         legend.text.align = 0,
         legend.title = element_blank( ),
         strip.background = element_blank(),
         strip.text = element_text( size = 18))
ggsave( deaths_coal_pm.gg,
        filename = 'figures/deaths_per_year_pm_coal_pm25.png',
        height = 4, width = 9, unit = 'in')
ggsave( deaths_coal_pm.gg,
        filename = 'figures/deaths_per_year_pm_coal_pm25.pdf',
        height = 4, width = 9, unit = 'in', device = cairo_pdf)

# fraction before 2008
deaths_by_year_merge[ year %in% 2000:2008, sum( deaths_hyads) / sum( deaths_pm)]
deaths_by_year_merge[ year %in% 2012:2016, sum( deaths_hyads) / sum( deaths_pm)]
deaths_by_year_merge[ , sum( deaths_hyads) / sum( deaths_pm)]
deaths_by_year_merge[ , sum( deaths_pm)]

## ================================================= ##
#  plot total hyads
## ================================================= ##
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
setnames( zips, 'zip', 'ZIP')

## CHANGE THIS ##
# gridded hyads file location
hyads_file_loc <- '/n/dominici_nsaph_l3/Lab/projects/analytic/coal_exposure_pm25/zips_model.lm.cv_single_poly'

#coordinate reference system projection string for spatial data
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# get the names of the gridded HyADS output files
zips.files.yr <- list.files( hyads_file_loc,
                             pattern = 'zips_pm25_total_\\d{4}\\.fst',
                             full.names = TRUE)

# read select files
zips.dat <- lapply( zips.files.yr,
                    function( f){
                      year.f <- gsub( '^.*_|\\.fst', '', f)
                      
                      in.f <- read.fst( f, as.data.table = T)
                      in.f[, year := year.f]
                    }) %>% rbindlist

# calculate average by year
zips.dat[, .( mean = mean( vals.out, na.rm = TRUE),
              X_00 = min( vals.out),
              X_75 = quantile( vals.out, .75, na.rm = TRUE),
              X_25 = quantile( vals.out, .25, na.rm = TRUE),
              X_100 = max( vals.out, na.rm = TRUE)), by = year]

# merge with spatial data
zips.dat.sf <- merge( zips.dat,
                      zips,
                      by.x = 'ZIP', by.y = 'ZIP')

# get USA dataset
states48 <- c( state.name[!(state.name %in% c( 'Alaska', 'Hawaii'))],
               'District of Columbia')
# download some US data
states <- USAboundaries::us_states()
states <- states[states$name %in% c( states48),]

# do the plpotting
spat.gg <- ggplot( zips.dat.sf[ zips.dat.sf$year %in% c( 1999, 2006, 2013, 2020),],
                   aes( fill = vals.out, color = vals.out)) + 
  geom_sf( aes( geometry = geometry),
           color = 'transparent') + 
  scale_fill_gradient(high = "black",
                      low = "cornsilk",
                      breaks = c( 0, 1, 2, 3),
                      labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                      na.value = NA,
                      oob = scales::squish,
                      limits = c( 0, 2)) +
  scale_color_gradient(high = "black",
                       low = "cornsilk",
                       breaks = c( 0, 1, 2, 3),
                       labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                       na.value = NA,
                       oob = scales::squish,
                       limits = c( 0, 2),
                       guide = 'none') +
  geom_sf( data = states,
           inherit.aes = FALSE,
           fill = NA, color = 'grey90', size = .05) +
  labs( fill =  expression(paste( Coal, ' ', PM["2.5"], ', µg ', m^{"-3"})),
        color = expression(paste( Coal, ' ', PM["2.5"], ', µg ', m^{"-3"}))) +
  scale_x_continuous( expand = c( 0, 0)) +
  scale_y_continuous( expand = c( 0, 0),
                      breaks = 1:10) +
  facet_wrap( . ~ year, nrow = 1, strip.position = 'bottom') +
  theme_bw() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         legend.background = element_rect( color = 'black'),
         legend.direction = 'vertical',
         legend.position = c( .93, 2),
         legend.text = element_text( size = 14),
         legend.title = element_text( size = 14),
         panel.background = element_rect( fill = 'white'),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 20))

hist.gg <- ggplot( zips.dat, 
                   aes( y = vals.out, x = year, group = year)) + 
  geom_boxplot( size = .5) +
  scale_y_continuous( breaks = 0:10) +
  labs( y = expression( paste( Coal, ' ', PM["2.5"], ', µg', m^{"-3"}))) +
  theme_bw() + 
  theme( axis.text = element_text( size = 16),
         # axis.text.x = element_text( size = 20),
         axis.title = element_text( size = 18),
         axis.title.x = element_blank())

gg_combine <- cowplot::plot_grid(hist.gg,
                                 spat.gg,
                                 rel_heights = c( 1.1,1),
                                 labels = NULL, ncol = 1, 
                                 align = 'v', axis = 'lr') 

ggsave( 'figures/hyads_trends.pdf', gg_combine,
        width = 14, height = 5, scale = 1.2, device = cairo_pdf)
ggsave( 'figures/hyads_trends.png', gg_combine,
        width = 14, height = 5, scale = 1.2)



## ================================================= ##
#  plot hyads vs raw hyads
## ================================================= ##
gg_raw_hyads_vs_coal_pm <- 
  ggplot( dat_annual[ year %in% c( 1999, 2006, 2011, 2020)],
          aes( x = Y1_raw, y = Y1)) + 
  geom_point( size = .1, alpha = .1) +
  labs( x = "Raw HyADS [unitless]", y = expression(Coal~PM['2.5']~'[µg '~m^{'-3'}*']')) +
  facet_wrap( . ~ year, nrow = 1) + 
  theme_minimal() + 
  theme( axis.text = element_text( size = 12),
         axis.title = element_text( size = 14),
         strip.text = element_text( size = 14))



## ================================================= ##
#  plot other pollutants
## ================================================= ##
dir_data <- 'data/data' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'
dir_out <- 'results/' #/nfs/home/X/xwu/shared_space/ci3_xwu/National_Causal/data2016_temp/'

aggregate_data <- read.fst( file.path(dir_data, "cache_data", "aggregate_data.fst"),
                            as.data.table = TRUE)
covariates <- 
  read.fst( file.path(dir_data, "cache_data", "covariates.fst"), as.data.table = TRUE)

# load hyads data
dat_annual <- read_fst( file.path(dir_data, "cache_data", 'hyads_pm25_annual.fst'),
                        columns = c('zip','year', 'Y1', 'Y1.adj'), as.data.table = TRUE)

# create death & rate sums by zip-year
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

# calculate portion of PM that is not coal PM
dat_annual_hy[, pm_not_coalPM := pm25_ensemble - Y1]

# create dataset for plotting
other_pollutants.dt <- 
  dat_annual_hy[ year %in% c( 2000, 2008, 2016),
                 .( zip, year, pm_not_coalPM, pm25_ensemble, 
                    ozone_summer.current_year, no2.current_year)] %>%
  # melt( id.vars = c( 'zip', 'year'), variable.name = 'pollutant', value.name = 'concentration') %>%
  merge( zips, by = 'zip')

gg_pm <- 
  ggplot( other_pollutants.dt,
          aes( geometry = geometry, fill = pm25_ensemble)) + 
  geom_sf( size = 0, color = NA) +
  geom_sf( data = states,
           inherit.aes = FALSE,
           fill = NA, color = 'grey90', size = .05) +
  scale_fill_gradient(high = "black",
                      low = "cornsilk",
                      # breaks = c( 0, 1, 2, 3),
                      # labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                      na.value = NA,
                      oob = scales::squish,
                      limits = c( 0, 15)) +
  guides( color = NULL) +
  labs( fill = expression( atop( Total~PM["2.5"], µg~m^{"-3"}))) +
  facet_wrap( . ~ year, nrow = 1, strip.position = 'top') +
  theme_bw() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         legend.background = element_blank( ),
         legend.direction = 'vertical',
         legend.position = 'right',
         legend.text = element_text( size = 14),
         legend.title = element_text( size = 14, hjust = 0),
         panel.background = element_rect( fill = 'white'),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         strip.text = element_text( size = 16))
gg_pm.nc <- 
  ggplot( other_pollutants.dt,
          aes( geometry = geometry, fill = pm_not_coalPM)) + 
  geom_sf( size = 0, color = NA) +
  geom_sf( data = states,
           inherit.aes = FALSE,
           fill = NA, color = 'grey90', size = .05) +
  scale_fill_gradient(high = "black",
                      low = "lightblue",
                      # breaks = c( 0, 1, 2, 3),
                      # labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                      na.value = NA,
                      oob = scales::squish,
                      limits = c( 0, 15)) +
  guides( color = NULL) +
  labs( fill = expression( atop( Non-coal~PM["2.5"], µg~m^{"-3"}))) +
  facet_wrap( . ~ year, nrow = 1, strip.position = 'bottom') +
  theme_bw() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         legend.background = element_blank( ),
         legend.direction = 'vertical',
         legend.position = 'right',
         legend.text = element_text( size = 14),
         legend.title = element_text( size = 14, hjust = 0),
         panel.background = element_rect( fill = 'white'),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         strip.text = element_blank( ))

gg_no2 <- 
  ggplot( other_pollutants.dt,
          aes( geometry = geometry, fill = no2.current_year)) + 
  geom_sf( size = 0, color = NA) +
  geom_sf( data = states,
           inherit.aes = FALSE,
           fill = NA, color = 'grey90', size = .05) +
  scale_fill_gradient(high = "black",
                      low = "pink",
                      # breaks = c( 0, 1, 2, 3),
                      # labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                      na.value = NA,
                      oob = scales::squish,
                      limits = c( 0, 30)) +
  labs( fill = expression( atop( NO["2"], ppb))) +
  guides( color = NULL) +
  facet_wrap( . ~ year, nrow = 1, strip.position = 'bottom') +
  theme_bw() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         legend.background = element_blank( ),
         legend.direction = 'vertical',
         legend.position = 'right',
         legend.text = element_text( size = 14),
         legend.title = element_text( size = 14, hjust = 0),
         panel.background = element_rect( fill = 'white'),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         strip.text = element_blank( ))

gg_o3 <- 
  ggplot( other_pollutants.dt[1:1000],
          aes( geometry = geometry, fill = ozone_summer.current_year)) + 
  geom_sf( size = 0, color = NA) +
  geom_sf( data = states,
           inherit.aes = FALSE,
           fill = NA, color = 'grey90', size = .05) +
  scale_fill_gradient(high = "black",
                      low = "bisque",
                      # breaks = c( 0, 1, 2, 3),
                      # labels = c( '0.0', '1.0', '≥2.0', '3.0'),
                      na.value = NA,
                      oob = scales::squish,
                      limits = c( 0, 75)) +
  guides( color = NULL) +
  labs( fill = expression( atop( O["3"], ppb))) +
  facet_wrap( . ~ year, nrow = 1, strip.position = 'bottom') +
  theme_bw() + 
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         legend.background = element_blank( ),
         legend.direction = 'vertical',
         legend.position = 'right',
         legend.text = element_text( size = 14),
         legend.title = element_text( size = 14, hjust = 0, debug = T),
         legend.title.align = 0,
         panel.background = element_rect( fill = 'white'),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         strip.background = element_blank(),
         strip.text = element_blank())
gg_o3

gg_pollutants <- 
  cowplot::plot_grid( gg_pm,
                      gg_pm.nc,
                      gg_no2,
                      gg_o3,
                      rel_heights = c( 1.11, 1, 1, 1),
                      labels = NULL, ncol = 1, 
                      align = 'v', axis = 'lr') 
# gg_pollutants
ggsave( 'figures/pollutant_maps.png', gg_pollutants,
        width = 14, height = 10, scale = .8)

## ================================================== ##
# what percent of all medicare deaths?
## ================================================== ##
sum_deaths_year <- read.fst( 'data/cache_data/total_deaths_by_state_year.fst', as.data.table = TRUE)

medicare_deaths_by_year <- 
  sum_deaths_year[, .( deaths = sum( deaths.fill)), by = year]
medicare_deaths_all <- 
  sum_deaths_year[, .( deaths = sum( deaths.fill))] %>% unlist
medicare_deaths_all.2016 <- 
  sum_deaths_year[ year <= 2016, .( deaths = sum( deaths.fill))] %>% unlist

# deaths in all years
deaths_all <- deaths_by_year[model == 'hyads', 
                             lapply( .SD, sum),
                             .SDcols = sum_cols] 
deaths_all / medicare_deaths_all

# deaths in all years
deaths_by_year[model == 'hyads' & year <= 2005, 
               lapply( .SD, sum),
               .SDcols = sum_cols]# / length( 1999:2005)
deaths_by_year[model == 'hyads' & year == 2020, 
               lapply( .SD, sum),
               .SDcols = sum_cols] 
deaths_all.2016 <- deaths_by_year[model == 'hyads' & year <= 2016, 
                                  lapply( .SD, sum),
                                  .SDcols = sum_cols] 
deaths_all.2016 / medicare_deaths_all.2016

# deaths in all years, PM
deaths_all.wu <- deaths_by_year[model == 'pm25_ensemble', 
                                lapply( .SD, sum),
                                .SDcols = sum_cols] 
deaths_all.kr <- deaths_by_year[model == 'pm25_krewski', 
                                lapply( .SD, sum),
                                .SDcols = sum_cols] 
deaths_all.wu / medicare_deaths_all
deaths_all.kr / medicare_deaths_all


## ====================================================== ##
# plots of sensitivity analysis coefficients
## ====================================================== ##
# this is calculated in 03_deaths_contributed.R
num_uniq_zip <- 34983

# read in the coefficients
poisson_coefs <- fread( paste0(dir_out, "poisson_model_coefs.csv"))

# read the confidence interval dataset
# loads the object bootstrap_CI
load( paste0(dir_out, "loglinear_coefs_boots.RData"))

# which ones to plot?
plot_betas <- 
  c( 'pm25_ensemble', 'hyads', 'hyads_early', 'hyads_mid', 'hyads_late',
     'hyads_pm', 'hyads_dpm', 'hyads_dpm_early', 'hyads_dpm_mid',
     'hyads_dpm_late', 'hyads_pollutants', 'hyads_pollutants_no_o3',
     'hyads_raw_pollutants_no_o3', 'hyads.adj', 'hyads.adj_pm', 'hyads.adj_dpm')

# Add descriptors
poisson_coefs_info <- 
  data.table(model = plot_betas,
             descriptor1 = c( 'PM', 
                              rep( 'Coal~PM["2.5"]~with~no~additional~adjustment', 4),
                              'Coal~PM["2.5"]~adjusted~by~total~PM["2.5"]',
                              rep( 'Coal~PM["2.5"]~adjusted~by~non-coal~PM["2.5"]', 4),
                              'Coal~PM["2.5"]~adjusted~by~non-coal~PM["2.5"]*","~O["3"]*","~and~NO["2"]',
                              'Coal~PM["2.5"]~adjusted~by~non-coal~PM["2.5"]~and~NO["2"]',
                              'Raw~HyADS~adjusted~by~non-coal~PM["2.5"]~and~NO["2"]',
                              rep( 'Sulfate-corrected~coal~PM["2.5"]', 3)),
             descriptor2 = c( '2000-2016',
                              '1999-2016', '1999-2003', '2004-2007', '2008-2016',
                              '2000-2016', '2000-2016',
                              '1999-2003', '2004-2007', '2008-2016',
                              '2000-2016', '2000-2016', '2000-2016',
                              'No~additional~adjustment~(1999-2016)', 
                              'Adjusted~by~total~PM["2.5"]~(2000-2016)', 
                              'Adjusted~by~non-coal~PM["2.5"]~(2000-2016)'),
             descriptor3 = c( 'PM',
                              rep( 'coalPM', 11),
                              'raw_coalPM',
                              rep( 'coalPM.adj', 3)))

# add in the total deaths from each model
deaths_total_models <- 
  read.fst( paste0(dir_out, "deaths_by_year_total_coal_pm25.fst"), as.data.table = TRUE)

deaths_total_models_sum <- 
  deaths_total_models[, lapply( .SD, sum),
                      .SDcols = c( 'deaths_coef_1', 'deaths_coef_2', 'deaths_coef_3'),
                      by = model]

# merge all the data
poisson_coefs_plot <- 
  merge( poisson_coefs,
         poisson_coefs_info, by = 'model') %>%
  merge( deaths_total_models_sum, by = 'model')

poisson_coefs_plot[, `:=` (descriptor1 = factor( descriptor1, levels = unique( poisson_coefs_info$descriptor1)),
                           descriptor2 = factor( descriptor2, levels = rev( unique( poisson_coefs_info$descriptor2))),
                           descriptor3 = factor( descriptor3, levels = rev( unique( poisson_coefs_info$descriptor3))))]

# function to calculate CIs
CI_calc <- 
  function( estimate, mod, highlow){
    
    vec_boots <- unlist( bootstrap_CI[, ..mod])
    
    ci <- 1.96 * sd( vec_boots) *
      sqrt( 2 * sqrt( num_uniq_zip)) / sqrt(num_uniq_zip)
    
    if( highlow == 'high')
      return( estimate + ci)
    else
      return( estimate - ci)
  }

# create data table of central estimates and CI's
poisson_coefs_plot[, `:=` (ci_low  = CI_calc( Estimate, model, 'low'),
                           ci_high = CI_calc( Estimate, model, 'high')), by = model]


#create the plot
coef_ci.gg <-
  ggplot( poisson_coefs_plot[  !( model %in% c( 'pm25_ensemble', 'hyads_raw_pollutants_no_o3'))],
          aes( y = descriptor2,
               x = exp( Estimate), xmin = exp( ci_low), xmax = exp( ci_high)#,
               # shape = descriptor3, linetype = descriptor3
          )) + 
  scale_y_discrete( labels = scales::label_parse()) +
  scale_x_continuous( name = expression( Risk~ratio~per~µg~m^{"-3"})) +
  geom_vline( xintercept = 1) +
  geom_point() + 
  geom_errorbarh( height = 0) + 
  ggforce::facet_col( descriptor1 ~ ., scales = "free_y", 
                      space = "free", labeller = label_parsed) +
  theme_minimal() + 
  theme( axis.title.y = element_blank(),
         axis.text = element_text( size = 12),
         axis.title.x = element_text( size = 14),
         strip.background = element_blank(),
         strip.clip = 'off',
         strip.text = element_text( size = 16, hjust = 0))
# coef_ci.gg 


# create a bar chart
death_ci.gg <-
  ggplot( poisson_coefs_plot[  !( model %in% c( 'pm25_ensemble', 'hyads_raw_pollutants_no_o3'))],
          aes( y = descriptor2,
               x = deaths_coef_2 / 10000, xmin = deaths_coef_1  / 10000, xmax = deaths_coef_3  / 10000
          )) + 
  geom_col( fill = 'grey80') +
  geom_errorbarh( height = 0) + 
  scale_y_discrete( labels = scales::label_parse()) +
  scale_x_continuous( name = 'Total excess deaths (10,000)',
                      breaks = seq( 0, 100, 25)) +
  ggforce::facet_col( descriptor1 ~ ., scales = "free_y", 
                      space = "free", labeller = label_parsed) +
  theme_minimal() + 
  theme( axis.title.y = element_blank(),
         axis.text = element_text( size = 12),
         axis.text.y = element_blank(),
         axis.title.x = element_text( size = 14),
         strip.background = element_blank(),
         strip.text = element_text( size = 16, color = 'transparent'))
# death_ci.gg

cowplot::plot_grid( coef_ci.gg,
                    death_ci.gg,
                    rel_widths = c( 2, 2),
                    labels = NULL, #unique( poisson_coefs_plot[  !( model %in% c( 'pm25_ensemble', 'hyads_raw_pollutants_no_o3'))]$descriptor1), 
                    ncol = 2, 
                    align = 'h', axis = 'lr') 

design <- "
  A##
  AB#
  #BC
  ##C
"

lvls <- c("Points", "Bar")
poisson_coefs_plot[, V1 := NULL]
dat.melt <- 
  melt( poisson_coefs_plot[  !( model %in% c( 'pm25_ensemble', 'hyads_raw_pollutants_no_o3'))],
        id.vars = c( 'model', 'descriptor1', 'descriptor2', 'descriptor3'))
dat.melt[variable %in% c( 'ci_low', 'Estimate', 'ci_high'), plot := 'point']
dat.melt[variable %in% c( 'deaths_coef_1', 'deaths_coef_2', 'deaths_coef_3'), plot := 'bar']

dat_point <- 
  poisson_coefs_plot[, .( model, descriptor1, descriptor2, descriptor3, Estimate, ci_low, ci_high)]
dat_point[, plot := 'point']
dat_bar <- 
  poisson_coefs_plot[, .( model, descriptor1, descriptor2, descriptor3, deaths_coef_1, deaths_coef_2, deaths_coef_3)]
dat_bar[, plot := 'bar']

dat.melt_point <- 
  rbind( dat_point, dat_bar, fill = TRUE)


ggplot( ) + 
  scale_y_discrete( labels = scales::label_parse()) +
  scale_x_continuous( name = expression( Risk~ratio~per~µg~m^{"-3"})) +
  # geom_vline( xintercept = 1,
  #             data = dat.melt[plot == 'point']) +
  geom_point( aes( y = descriptor2,
                   x = exp( Estimate)),
              data = dat.melt_point[plot == 'point']) + 
  geom_errorbarh( aes( y = descriptor2,
                       xmin = exp( ci_low), xmax = exp( ci_high)),
                  height = 0,
                  data = dat.melt_point[plot == 'point']) +
  geom_col( aes( y = descriptor2,
                 x = deaths_coef_2 / 10000),
            data = dat.melt_point[plot == 'bar']) + 
  # facet_grid2(  descriptor1 ~ facet, scales = "free",
  #               labeller = label_parsed, space = 'free'
  #               # dir = 'h'
  #               )  +
  # force_panelsizes(cols = c(5, 1),
  #                  # rows = c(0.5),
  #                  respect = TRUE) +
  facet_grid( descriptor1 ~ plot, scales = "free",
              space = "free",
              labeller = label_parsed) +
  # ggforce::facet_row( descriptor1 ~ plot, scales = "free_y",
  #                     space = "free", labeller = label_parsed) +
  theme_minimal() +
  theme( axis.title.y = element_blank(),
         axis.text = element_text( size = 12),
         axis.title.x = element_text( size = 14),
         strip.text = element_text( size = 14))

## ====================================================== ##
# how much do mortality rates vary? 
#. approach - look at largest ZIP codes in different regions
## ====================================================== ##
dat_annual_zip_death

## ================================================== ##
# present some sums
## ================================================== ##
deaths_by_year[model == 'hyads' & year < 2009,
               lapply( .SD, sum),
               .SDcols = sum_cols]

dim( deaths_by_fac_statebinf[deaths_coef_2 > 5000])
dim( deaths_by_fac_statebinf[deaths_coef_2 > 1000])
deaths_by_fac_statebinf[deaths_coef_2 > 500 & statebin_facility == 'West']

# how much deaths associated with 50% emissions?
# merge facility names with facilities dataset
units_sum <- units_all[ , .( SOx = sum( SOx)), by = FacID]
units_sum[, FacID := as.integer( FacID)]

deaths_by_fac_statebinf_emiss <- 
  merge( units_sum,
         deaths_by_fac_statebinf[model == 'hyads'],
         by = 'FacID')

# rank the emissions and death
deaths_by_fac_statebinf_emiss[, `:=`( emiss.rank = frank( SOx),
                                      death.rank = frank( deaths_coef_2),
                                      emiss.frac = SOx / sum( SOx),
                                      death.frac = deaths_coef_2 / sum( deaths_coef_2))]

# fraction of deaths from facilities with lowest 50% emissions
fac_SOx_50l <- deaths_by_fac_statebinf_emiss[emiss.rank <= .N / 2, sum( deaths_coef_2)]
fac_SOx_50u <- deaths_by_fac_statebinf_emiss[emiss.rank > .N / 2, sum( deaths_coef_2)]

fac_SOx_50u / ( fac_SOx_50l + fac_SOx_50u)
( fac_SOx_50l + fac_SOx_50u) == alldeaths$deaths.med

perc.rank <- function(x) trunc(rank(x))/length(x)

deaths_by_fac_statebinf_emiss[, `:=` ( SOx_pct = perc.rank( SOx),
                                       deaths_pct = perc.rank( deaths_coef_2))]#,
# by = statebin_facility]
setkey( deaths_by_fac_statebinf_emiss, SOx_pct)


# calculate vertical distance from 
deaths_by_fac_statebinf_emiss[, `:=` ( min.plot.y = min(deaths_pct, SOx_pct ),
                                       max.plot.y = max(deaths_pct, SOx_pct )),
                              by = FacID]

pct_emiss_death.gg <- 
  ggplot( deaths_by_fac_statebinf_emiss,
          aes( x = SOx_pct, y = deaths_pct,
               ymin = min.plot.y,
               ymax = max.plot.y)) + 
  labs( x = expression( SO[2]~emissions~percentile),
        y = 'Associated Medicare deaths percentile') +
  geom_abline( slope = 1, intercept = 0,
               color = 'grey80') +
  geom_linerange() +
  geom_point() + 
  scale_x_continuous( labels = scales::percent_format()) +
  scale_y_continuous( labels = scales::percent_format()) +
  facet_wrap( . ~ statebin_facility, nrow = 2) + 
  theme_minimal() +
  theme( axis.text = element_text( size = 12),
         axis.title = element_text( size = 16),
         legend.position = 'none',
         panel.grid.minor = element_blank(),
         strip.text = element_text( size = 14))

# save it
ggsave( pct_emiss_death.gg,
        filename = 'figures/deaths_emiss_pct.png',
        height = 7, width = 15, unit = 'in')

# check out Gibson (6113) and Muskingum River (2872)
gibs_musk <- 
  deaths_by_fac_year_statebin_lab[FacID %in% c( 2872, 6113) & 
                                    model == 'hyads']
gibs_musk[ year <= 2006, mean( deaths_coef_2), by = FacID]
gibs_musk[ year >= 2007 & year <= 2015, .( deaths_coef_1 = sum( deaths_coef_1), 
                                           deaths_coef_2 = sum( deaths_coef_2),
                                           deaths_coef_3 = sum( deaths_coef_3)), by = FacID]

# check out Keystone (3136)
keystone <- 
  deaths_by_fac_year_statebin_lab[FacID == 3136 & 
                                    model == 'hyads']
keystone[ year <= 2008, .( deaths_coef_1 = mean( deaths_coef_1), 
                           deaths_coef_2 = mean( deaths_coef_2),
                           deaths_coef_3 = mean( deaths_coef_3)), by = FacID]
keystone[ year >= 2010, .( deaths_coef_1 = mean( deaths_coef_1), 
                           deaths_coef_2 = mean( deaths_coef_2),
                           deaths_coef_3 = mean( deaths_coef_3)), by = FacID]
