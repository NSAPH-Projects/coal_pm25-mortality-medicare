library( data.table)
library( ggplot2)
library( viridis)
library( readxl)
library( sf)
library( ncdf4)
library( raster)
library( fst)
library( tidyverse)
library( cowplot)


#coordinate reference system projection string for spatial data
p4s <-   '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m'

crs.use   <-  sf::st_crs( "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m")

# load major metropolitan areas
# from https://factfinder.census.gov/faces/tableservices/jsf/pages/productview.xhtml?src=bkmk
met.by.pop <- fread( '~/Dropbox/Harvard/Manuscripts/MobileBiasMS/metro_populations/PEP_2018_PEPANNCHG.US24PR/PEP_2018_PEPANNCHG.US24PR_with_ann.csv',
                     check.names = T)
met.by.pop.descript <- met.by.pop[1]
met.by.pop <- met.by.pop[-1]
met.by.pop[ ,GC.display.label := gsub( "[^[:alnum:]]", ' ', GC.display.label)]
met.by.pop <- met.by.pop[ -grep( ' PR Metro', GC.display.label)]

met.by.pop.X <- met.by.pop[ as( rank72018, 'numeric') <= 20]
MSA.codes <- unique( met.by.pop.X$GC.target.geo.id2)

# load counties spatial object
counties.all <- USAboundaries::us_counties( ) %>%
  st_transform( p4s)
counties.all$state_name <- NULL

# read in table of msa codes
county.msa <- data.table( read_excel( '~/Dropbox/Harvard/Manuscripts/MobileBiasMS/metro_populations/list1_Sep_2018.xls',
                                      skip = 2))
county.use <- county.msa[ `CBSA Code` %in% MSA.codes & !is.na( `CBSA Code`)]

# merge to get county-MSA spatial object
county.msa.sf <- sf::st_as_sf( merge( county.use, counties.all,
                                      by.y = c( 'statefp', 'countyfp'), 
                                      by.x = c( 'FIPS State Code', 'FIPS County Code')))
county.msa.sf <- sf::st_transform( county.msa.sf, crs.use)

# combine counties into single CSA's
msa.sf <- county.msa.sf %>%
  group_by( `CBSA Code`, `CBSA Title`) %>%
  summarize %>%
  st_as_sf()

## ================================================================ 
# take over HyADS
## ================================================================ 
# read in hyads
# gridded hyads file location
# gridded hyads file location
hyads_file_loc <- '~/Dropbox/Harvard/ARP/HyADS/hyads_longterm/exp_pm25_noint/grids_model.lm.cv_single_poly'

# get the names of the gridded HyADS output files
grid.files.unit.yr <- list.files( hyads_file_loc,
                                  pattern = 'grids_pm25_byunit_\\d{4}\\.fst',
                                  full.names = TRUE)

# read select files
grid.unit.dat <- lapply( grid.files.unit.yr,
                         function( f){
                           print( f)
                           year.f <- gsub( '^.*_|\\.fst', '', f) %>%
                             as( 'integer')
                           
                           in.f <- read.fst( f, as.data.table = T)
                           in.f[, year := year.f]
                           
                           in.f[ is.na( in.f)] <- 0
                           return (in.f)
                         }) %>% rbindlist( fill = T)

## ================================================================ 
# facility data
## ================================================================ 
## read the unit dataset
annual_dat <- fread( '~/Dropbox/Harvard/ARP/Data_AMPD_EIA/AMPD_Unit.csv')

units_and_facs_locs <- unique( annual_dat[, .( Facility.Latitude,
                                               Facility.Longitude,
                                               Facility.Name,
                                               State)])


## ================================================================ 
# area weight over csas
## ================================================================ 
hyads_unit_csa <- 
  lapply( unique( grid.unit.dat$year),
          function( yy){
            print( yy)
            grid.unit.dat.y.sf <- grid.unit.dat[ year == yy] %>%
              rasterFromXYZ( crs = p4s) %>%
              rasterToPolygons() %>%
              st_as_sf() %>%
              st_interpolate_aw( st_transform( msa.sf, p4s), extensive = F)
            
            # convert to data.table
            grid.unit.dat.y.dt <- as.data.table( grid.unit.dat.y.sf)
            grid.unit.dat.y.dt[, `:=` ( year = NULL,
                                        geometry = NULL,
                                        cbsa = msa.sf$`CBSA Code`,
                                        cbsa.n = msa.sf$`CBSA Title`)]
            
            grid.unit.dat.y.m <- melt( grid.unit.dat.y.dt, 
                                       measure.vars = grep( '^X', names( grid.unit.dat.y.dt),
                                                            value = T),
                                       id.vars = c( 'cbsa', 'cbsa.n'),
                                       value.name = 'hyads',
                                       variable.name = 'uID')
            
            # add year and uID
            grid.unit.dat.y.m[, `:=` (year = yy,
                                      uID = gsub( '^X', '', uID))]
            
            # merge with & sum by facility info
            unit_fac_xwalk <- fread( '~/Dropbox/GeorgeMason/Grants/2020/EPA_Empower/data/unit_fac_crosswalk.csv')
            unit_fac_xwalk_edit <- unit_fac_xwalk[ !duplicated( uID), 
                                                   .( FacID, uID, Facility.Name, State)] %>%
              unique()
            
            grid.unit.dat.y.m.f <- merge( grid.unit.dat.y.m, unit_fac_xwalk_edit,
                                          by = 'uID')
            
            # sum by facility data
            fac.dat <- grid.unit.dat.y.m.f[, .( hyads = sum( hyads)),
                                           by = .( cbsa, cbsa.n, year, FacID,
                                                   Facility.Name, State)]
            
            # do the ranking
            fac.dat[, `:=` ( rank = frankv( hyads, order = -1, ties.method = 'first')), 
                    by = cbsa]
            
            return( fac.dat)
          }) %>% rbindlist()

hyads_unit_csa[is.na( hyads), hyads := 0]
hyads_unit_csa[rank == 1 & year <= 2010 & 
                 cbsa %in% c( 12060, 47900, 14460, 16980,
                              35620, 19100, 37980, 33460)] %>%
  summary()

facs <- hyads_unit_csa[ rank <= 1 & 
                          cbsa == 47900]$FacID %>%
  unique()

hyads_unit_csa[, Facility.Name := gsub( 'Generating', 'Gen', Facility.Name)]
units_and_facs_locs[, Facility.Name := gsub( 'Generating', 'Gen', Facility.Name)]
# choose facilities with highest population-weighted exposure,
# plot area plot for all 

# goals
# individual facililties have potentially large impact
# facility impacts change over time
# different facilities are important in different regions
# (accompany each plot with a map?)
# maybe choose 3 population centers, 
#  -isolate ~5 most important facilities
#  -show these facilities' impacts on top of total impacts
# map the important facilities

## ================================================================ 
# get some spatial info about states, facilities
## ================================================================ 
# get states shapefile
states_use <- state.abb[ !( state.abb %in% c( 'HI', 'AK'))]
us <- USAboundaries::us_states()
states <- us[ us$state_abbr %in% states_use,]

# get all facilities
fac_locs_all <- st_as_sf( units_and_facs_locs[ Facility.Name %in% na.omit( hyads_unit_csa)$Facility.Name & 
                                                 !( State %in% 'PR')],
                          coords = c( 'Facility.Longitude', 'Facility.Latitude'))
st_crs( fac_locs_all) <- sf::st_crs( "+proj=longlat +datum=WGS84 +no_defs")


## ================================================================ 
# plot csa's using function
## ================================================================ 
csa_plotter <- function( csa = c( 47900, 12060),
                         x_dim = 12,
                         y_dim = 4,
                         y_upper = NA
){
  # first, get dt's of facilities ranked first at some point
  facs_first <- 
    lapply( csa,
            function( .csa){
              print( .csa)
              # identify facilities
              hyads_facs <- hyads_unit_csa[cbsa == .csa]
              facs <- hyads_facs[ rank <= 1]$FacID %>%
                unique()
              
              # naming conventions for summing
              hyads_facs[ FacID %in% facs, Fac.name := Facility.Name]
              hyads_facs[ is.na( Fac.name), Fac.name := 'All others']
              
              # sum by facilities & other facilities
              hyads_facs.f <- 
                hyads_facs[, .( hyads = sum( hyads)),
                           by = .( year, Fac.name)]
              
              # merge with location information
              hyads_facs.f.loc <- 
                merge( hyads_facs.f,
                       units_and_facs_locs[!duplicated( Facility.Name, State)], 
                       by.y = 'Facility.Name',
                       by.x = 'Fac.name', all.x = TRUE)
              hyads_facs.f.loc[ is.na( State), State := '']
              
              # naming convention, include states in names
              hyads_facs.f.loc[ Fac.name != 'All others', Fac.name.use := paste( Fac.name, State, sep = ', ')]
              hyads_facs.f.loc[ Fac.name == 'All others', Fac.name.use := Fac.name]
              
              # add in cbsa info
              hyads_facs.f.loc[, `:=` ( cbsa = .csa,
                                        cbsa.n = unique( hyads_facs$cbsa.n))]
              
              # sf object of coal facilities
              fac_locs <- 
                hyads_facs.f.loc[, .( Facility.Longitude, 
                                      Facility.Latitude, 
                                      Fac.name.use)] %>%
                unique() %>%
                na.omit() %>%
                st_as_sf( coords = c( 'Facility.Longitude', 'Facility.Latitude'))
              st_crs( fac_locs) <- sf::st_crs( "+proj=longlat +datum=WGS84 +no_defs")
              
              #Create a custom color scale
              library(RColorBrewer)
              myColors <- brewer.pal( nrow( fac_locs),"Dark2")[1:nrow( fac_locs)]
              names( myColors) <- sort( fac_locs$Fac.name.use)
              myColors <- c( myColors, 'All others' = 'grey70')
              colScale <- scale_color_manual(name = "grp",values = myColors)
              filScale <- scale_fill_manual(name = "grp",values = myColors)
              
              # fix facet levels for faceting
              hyads_facs.f.loc[ , Fac.name.use := factor( Fac.name.use, levels = names( myColors))]
              
              # select the appropriate msa
              msa.plot <- msa.sf[msa.sf$`CBSA Code` == .csa,] 
              
              # make sure everything is on same crs
              crs.use <- "+proj=longlat +datum=WGS84 +no_defs"
              fac_locs_all <- st_transform( fac_locs_all, crs.use)
              msa.plot <- st_transform( msa.plot, crs.use)
              fac_locs <- st_transform( fac_locs, crs.use)
              states <- st_transform( states, crs.use)
              
              # get an appropriate boundary box
              bound_box.o <- st_bbox( fac_locs)
              bound_box <- 
                c( bound_box.o$xmin - .5 * ( x_dim - ( bound_box.o$xmax - bound_box.o$xmin)),
                   bound_box.o$xmax + .5 * ( x_dim - ( bound_box.o$xmax - bound_box.o$xmin)),
                   bound_box.o$ymin - .5 * ( y_dim - ( bound_box.o$ymax - bound_box.o$ymin)),
                   bound_box.o$ymax + .5 * ( y_dim - ( bound_box.o$ymax - bound_box.o$ymin)))
              
              # create area time series
              gg_area <-
                ggplot( hyads_facs.f.loc,
                        aes( x = year, y = hyads, fill = Fac.name.use)) + 
                geom_area( color = NA) + 
                filScale +
                scale_y_continuous( labels = scales::number_format(accuracy = 0.1),
                                    expand = c( 0, 0),
                                    limits = c( NA, y_upper)) +
                scale_x_continuous( expand = c( 0, 0)) +
                labs( y = expression( paste( Coal, ' ', PM["2.5"], ', µg', m^{"-3"}))) +
                # facet_wrap( . ~ cbsa.n, strip.position = 'top') + 
                theme_bw() +
                theme( axis.text = element_text( size = 20),
                       axis.title.y = element_text( size = 20),
                       axis.title.x = element_blank(),
                       legend.background = element_blank(),
                       legend.position = c( .5, 1),
                       legend.justification = c( 0, 1),
                       legend.title = element_blank(),
                       legend.text = element_text( size = 14),
                       plot.margin = unit( c( 0, .25, 0, 0), 'in'),
                       strip.background = element_rect( fill = 'white'),
                       strip.text = element_text( size = 20, face = 'bold'))
              
              # plot units
              gg_units <-
                ggplot( ) + 
                geom_sf( data = states, fill = NA, size = .5) +
                geom_sf( data = fac_locs_all, color = 'grey70', alpha = .4) +
                geom_sf( data = msa.sf[msa.sf$`CBSA Code` == .csa,], fill = 'lightpink', alpha = .6) + 
                geom_sf( data = fac_locs,
                         aes( color = Fac.name.use), size = 4) + 
                colScale + 
                coord_sf( xlim = c( bound_box['xmin'], bound_box['xmax']),
                          ylim = c( bound_box['ymin'], bound_box['ymax'])) +
                theme_bw() + 
                theme( axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       legend.position = 'none',
                       panel.grid = element_blank(),
                       panel.border = element_blank(),
                       plot.margin = unit( c( 0, 0, 0, 0), 'in')
                       )
              
              # cowplot them together
              gg_combine <- cowplot::plot_grid(gg_area,
                                               gg_units, 
                                               rel_widths = c( 1.2, 1),
                                               labels = NULL, nrow = 1, 
                                               align = 'h', axis = 'tb')  +
                theme( plot.background = element_rect( fill = 'white',
                                                      color = "black", size = 1),
                       plot.margin = unit( c( .4, .1, .1, .1), 'in'))
              return( gg_combine)
            })
  
  msa.plot <- sapply( csa,
                      function( c)
                        msa.sf[ which( msa.sf$`CBSA Code` %in% c),]$`CBSA Title`)
  
  
  out <- plot_grid( plotlist = facs_first, ncol = 1,
                    align = 'v', axis = 'l', 
                    labels = msa.plot,
                    label_size = 20,
                    label_x = .5, hjust = .5)
  
  return( out)
}


## ================================================================ 
# Figures S1 & S2: plot csa's using function
## ================================================================ 
# do DC, ATL, Boston, chicago
csa_dc_atl_bos <- c( 12060, 47900, 14460, 16980)
csa_dc_nyc_dal_phi_min <- c( 35620, 19100, 37980, 33460)
gg_4 <- csa_plotter( csa_dc_atl_bos,
                     x_dim = 18,
                     y_dim = 8,
                     y_upper = 7)
ggsave( 'figures/hyads_facilities_csas1.png', gg_4,
        width = 7, height = 8, scale = 2.1)

gg_4_2 <- csa_plotter( csa_dc_nyc_dal_phi_min,
                       x_dim = 18,
                       y_dim = 8,
                       y_upper = 7)
ggsave( 'figures/hyads_facilities_csas2.png', gg_4_2,
        width = 7, height = 8, scale = 2.1)

##
csa_hou_dal_phx_stl <- c( '26420', '19100', '38060', '41180')
gg_4_3 <- csa_plotter( csa_hou_dal_phx_stl,
                       x_dim = 18,
                       y_dim = 8)
ggsave( '~/Dropbox/GeorgeMason/Research/Manuscripts/2021/HyADS_20yrs/figures/hyads_facilities_csasYX_20220920.png', gg_4_3,
        width = 7, height = 4, scale = 2.1)




