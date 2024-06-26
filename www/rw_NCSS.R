# title         : rw_NCSS.R
# purpose       : Reading and writing of NSCD data;
# reference     : NCSS Characterization Database [http://ssldata.nrcs.usda.gov/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : "Repository2010.mdb" MS Access dbs
# outputs       : R data frames for SoilProfileCollection;
# remarks 1     : LARGE dataset!

## Download the database:
# download.file("http://globalsoilmap.net/data/Repository2010.7z", destfile=paste(getwd(),"Repository2010.7z",sep="/"))

library(RODBC)
library(aqp)
library(sp)
# define a new function to merge the degree, min, sec columns:
cols2dms <- function(x,y,z,e){ifelse(is.na(e)|is.na(x), NA, as(char2dms(paste(x, "d", y, "'", z, "\"", e, sep="")), "numeric"))}
# load("Repository2010.RData")

# ------------------------------------------------------------
# Fetch tables - NCSS
# ------------------------------------------------------------

cNCSS <- odbcConnect(dsn="NCSS")
sqlTables(cNCSS)$TABLE_NAME
# get tables:
site <- sqlFetch(cNCSS, "site", stringsAsFactors=FALSE)
str(site)  # 37,333 profiles!!!
pedon <- sqlFetch(cNCSS, "pedon", stringsAsFactors=FALSE)
str(pedon)
pedon$observation_date <- as.character(as.Date(pedon$observation_date))
# has to remove Date format otherwise it gets problems with NA's
PSDA <- sqlFetch(cNCSS, "PSDA and Rock Fragments", stringsAsFactors=FALSE)   ## gravel content ("wpg2") is in mass percentage and needs to be converted to volume %
str(PSDA)
Organic <- sqlFetch(cNCSS, "Organic", stringsAsFactors=FALSE)
str(Organic)
CEC <- sqlFetch(cNCSS, "CEC and Bases", stringsAsFactors=FALSE)
str(CEC)
Carbon <- sqlFetch(cNCSS, "Carbon and Extractions", stringsAsFactors=FALSE)
str(Carbon)
BulkDens <- sqlFetch(cNCSS, "Bulk Density and Moisture", stringsAsFactors=FALSE)
str(BulkDens)  ## we need "Bulk Density, <2mm Fraction, Ovendry"
pH <- sqlFetch(cNCSS, "pH and Carbonates", stringsAsFactors=FALSE)
str(pH)
layer <- sqlFetch(cNCSS, "layer", stringsAsFactors=FALSE)
str(layer)
tax <- sqlFetch(cNCSS, "Taxonomy_New", stringsAsFactors=FALSE)
tax$last_correlated_date <- as.character(as.Date(tax$last_correlated_date))
Phosphorus <- sqlFetch(cNCSS, "Phosphorus", stringsAsFactors=FALSE)
MajorElements <- sqlFetch(cNCSS, "Major Elements", stringsAsFactors=FALSE)
str(MajorElements)
Salt <- sqlFetch(cNCSS, "Salt", stringsAsFactors=FALSE)
str(Salt)
TaxonomyE <- sqlFetch(cNCSS, "Taxonomy_Error", stringsAsFactors=FALSE)
str(TaxonomyE)
TaxonomyE$last_correlated_date <- as.character(as.Date(TaxonomyE$last_correlated_date))
procs <- sqlFetch(cNCSS, "analysis_procedure", stringsAsFactors=FALSE)
# save(list=c("Organic", "site", "pedon", "PSDA", "CEC", "Carbon", "BulkDens", "pH", "layer", "tax", "Phosphorus", "MajorElements", "Salt", "TaxonomyE"), file="cNCSS.RData")


# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

# add missing columns:
layer$DEPTH <- layer$hzn_top + (layer$hzn_bot - layer$hzn_top)/2
# mask out profiles with missing coordinates:
site <- site[!is.na(site$longitude_degrees)&!is.na(site$latitude_degrees),]
site$longitude_it <- ifelse(site$longitude_direction=="west", "W", "E")
site$latitude_it <- ifelse(site$latitude_direction=="north"|site$latitude_direction=="North", "N", "S")
site$LAT <- cols2dms(site$latitude_degrees, site$latitude_minutes, site$latitude_seconds, site$latitude_it)
site$LON <- cols2dms(site$longitude_degrees, site$longitude_minutes, site$longitude_seconds, site$longitude_it)
summary(site$LAT) # 253 NA's
summary(site$LON)
summary(site$latitude_seconds)
str(site)

# ------------------------------------------------------------
# create the horizon and site tables (takes time!)
# ------------------------------------------------------------
 
h1 <- merge(PSDA, CEC[,!(names(CEC) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h2 <- merge(Organic[,!(names(Organic) %in% c("result_source_key", "prep_code", "c_tot", "oc", "n_tot", "c_n_ra"))], Carbon, by=c("natural_key"), all=TRUE)
h3 <- merge(h1, h2[,!(names(h2) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
h4 <- merge(h3, layer, by=c("natural_key"), all=TRUE)
h5 <- merge(h4[,!(names(h4) %in% c("db_od"))], BulkDens[,!(names(BulkDens) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
horizon <- merge(h5[,!(names(h5) %in% c("ph_h2o"))], pH[,!(names(pH) %in% c("result_source_key", "prep_code"))], by=c("natural_key"), all=TRUE)
# names(horizon)
## fix some columns:
horizon$wpg2 <- ifelse(horizon$wpg2 > 100|horizon$wpg2 <0 , NA, horizon$wpg2)
## merge BulkDens and PSDA and derive GRAVEL content:
mBD <- (horizon$wpg2/100 * 2.6 + (1-horizon$wpg2/100) * horizon$db_od)
horizon$GRAVEL <- horizon$wpg2*mBD/2.6
## check visually:
# plot(y=horizon$GRAVEL, x=horizon$wpg2, xlim=c(0,100), ylab="GRAVEL (vol %)", xlab="GRAVEL (mass %)", main="NCSS (250K measurements)", pch="+", cex=.7)

pedon$site_key <- as.factor(pedon$site_key)
# there are also pedons with multiple site IDs!?
s1 <- merge(pedon, tax[,!(names(tax) %in% c("natural_key"))], by=c("pedon_key"), all=TRUE)
# summary(s1$site_key)  
s1 <- s1[!duplicated(s1$site_key),]
str(s1)
# merge the site table:
site$site_key <- as.factor(site$site_key)
site.s <- merge(site, s1, by=c("site_key"), all.x=TRUE)
site.s <- site.s[!is.na(site.s$LAT)&!is.na(site.s$LON), c("site_key","user_site_id","LAT","LON","horizontal_datum_name","observation_date","sampled_taxon_name","correlated_taxon_name","correlated_taxon_kind","correlated_class_name","SSL_class_name")]
str(site.s)

# subset tables:
hs <- subset(horizon[,c("site_key", "natural_key", "layer_sequence", "hzn_desgn", "hzn_bot", "hzn_top", "hzn_vert_subdvn", "clay_tot_psa", "sand_tot_psa", "silt_tot_psa", "pyr_col", "wpg2", "caco3", "ph_hist", "oc", "c_tot", "base_sum", "cec_sum", "cec_nh4", "ph_h2o", "ph_kcl", "db_od")], !is.na(horizon$DEPTH)&!is.na(horizon$site_key)&horizon$layer_sequence<15)
# strange - there are layer_sequence numbers up to 259?!
str(hs)
# summary(as.factor(horizon$layer_sequence))
hs$site_key <- as.factor(hs$site_key)
# remove duplicate site keys (there are many!!):
hs$layer_ID <- as.factor(paste(hs$site_key, hs$layer_sequence, sep="_"))
hs <- hs[!duplicated(hs$layer_ID),]
# remove IDs that do not exist in the site table:
sel <- hs$site_key %in% site.s$site_key
length(levels(as.factor(hs$layer_ID)))
mean(hs$oc, na.rm=TRUE) ## Organic Carbon in percent!
## estimate OC by correcting for CaCO3:
hs$ORCDRC <- 10*ifelse(!is.na(hs$c_tot), ifelse((hs$ph_h2o > 7)&!is.na(hs$caco3), hs$c_tot - .12 * hs$caco3, hs$c_tot), hs$oc)  
hs$ORCDRC <- ifelse(hs$ORCDRC < 0, 0, hs$ORCDRC) 
## stats:
# summary(hs$ORCDRC)


# ------------------------------------------------------------
# SoilProfileCollection
# ------------------------------------------------------------

NCSS_all <- list(sites=site.s, horizons=hs[sel,])
str(NCSS_all)
# specify classes:
NCSS_all$sites$user_site_id <- as.factor(NCSS_all$sites$user_site_id)
NCSS_all$sites$horizontal_datum_name <- as.factor(NCSS_all$sites$horizontal_datum_name)
NCSS_all$sites$correlated_taxon_name <- as.factor(NCSS_all$sites$correlated_taxon_name)
NCSS_all$sites$sampled_taxon_name <- as.factor(NCSS_all$sites$sampled_taxon_name)
NCSS_all$sites$correlated_taxon_kind <- as.factor(NCSS_all$sites$correlated_taxon_kind)
NCSS_all$sites$correlated_class_name <- as.factor(NCSS_all$sites$correlated_class_name)
NCSS_all$sites$SSL_class_name <- as.factor(NCSS_all$sites$SSL_class_name)
NCSS_all$horizons$hzn_desgn <- as.factor(NCSS_all$horizons$hzn_desgn)
NCSS_all$horizons$natural_key <- as.factor(NCSS_all$horizons$natural_key)
NCSS_all$horizons$pyr_col <- NULL
# round up the numbers:
NCSS_all$horizons$ORCDRC <- round(NCSS_all$horizons$ORCDRC, 1)
NCSS_all$horizons$clay_tot_psa <- round(NCSS_all$horizons$clay_tot_psa, 0)
NCSS_all$horizons$sand_tot_psa <- round(NCSS_all$horizons$sand_tot_psa, 0)
NCSS_all$horizons$silt_tot_psa <- round(NCSS_all$horizons$silt_tot_psa, 0)
NCSS_all$horizons$oc <- round(NCSS_all$horizons$oc, 2)
NCSS_all$horizons$c_tot <- round(NCSS_all$horizons$c_tot, 2)
NCSS_all$horizons$base_sum <- round(NCSS_all$horizons$base_sum, 1)
NCSS_all$horizons$caco3 <- round(NCSS_all$horizons$caco3, 2)
NCSS_all$horizons$cec_sum <- round(NCSS_all$horizons$cec_sum, 1)
NCSS_all$horizons$cec_nh4 <- round(NCSS_all$horizons$cec_nh4, 1)
NCSS_all$horizons$ph_h2o <- round(NCSS_all$horizons$ph_h2o, 1)
NCSS_all$horizons$ph_kcl <- round(NCSS_all$horizons$ph_kcl, 1)
NCSS_all$horizons$db_od <- round(NCSS_all$horizons$db_od, 2)
save(NCSS_all, file="NCSS_all.rda", compress="xz")

## Promote to SoilProfileCollection
# NCSS.spc <- join(NCSS_all$horizons, NCSS_all$sites, type='inner')
# depths(NCSS.spc) <- site_key ~ hzn_top + hzn_bot
## TAKES 4-5 mins!
## extract site data
# site(NCSS.spc) <- ~ LON + LAT + horizontal_datum_name + observation_date + sampled_taxon_name + correlated_taxon_name + correlated_taxon_kind + correlated_class_name + SSL_class_name
## generate SpatialPoints
# coordinates(NCSS.spc) <- ~ LON + LAT
## assign CRS data
# proj4string(NCSS.spc) <- "+proj=latlong +datum=NAD83"
# str(NCSS.spc)
# library(plotKML)
# kml(NCSS.spc, var.name="ORCDRC")

## subset to Indiana state:
sel2 <- NCSS_all$sites$LAT>36.9&NCSS_all$sites$LAT<42.5&NCSS_all$sites$LON< -84.5&NCSS_all$sites$LON> -91.7
site.s_in <- NCSS_all$sites[sel2,]
# remove some duplicated "user_site_id"s:
site.s_in <- site.s_in[!duplicated(site.s_in$user_site_id),]
hs_in <- NCSS_all$horizons[NCSS_all$horizons$site_key %in% site.s_in$site_key,]
# add the user_site_id column:
hs_in <- merge(hs_in, site.s_in[,c("site_key","user_site_id")])
ilin <- list(sites=site.s_in[,c("user_site_id","LAT","LON","observation_date","correlated_taxon_name","correlated_taxon_kind","correlated_class_name")], horizons=hs_in[,c("user_site_id", "layer_sequence", "hzn_desgn", "hzn_bot", "hzn_top", "clay_tot_psa", "sand_tot_psa", "silt_tot_psa", "wpg2", "ph_h2o", "ph_kcl", "ORCDRC", "db_od")])
# remove exta levels:
ilin$sites$user_site_id <- as.factor(paste(ilin$sites$user_site_id))
ilin$horizons$user_site_id <- as.factor(paste(ilin$horizons$user_site_id))
ilin$sites$correlated_taxon_name <- as.factor(paste(ilin$sites$correlated_taxon_name))
ilin$sites$correlated_class_name <- as.factor(paste(ilin$sites$correlated_class_name))
ilin$horizons$hzn_desgn <- as.factor(paste(ilin$horizons$hzn_desgn))
names(ilin$sites) <- c("SOURCEID", "LATNAD83", "LONNAD83", "TIMESTRR", "TAXSUSDA", "TAXKUSDA", "TAXNUSDA")
names(ilin$horizons) <- c("SOURCEID", "LSQINT", "HZDUSD", "UHDICM", "LHDICM", "CLYPPT", "SNDPPT", "SLTPPT", "CRFMAS", "PHIHO5", "PHIKCL", "ORCDRC", "BLDODR")
ilin$sites$TIMESTRR <- as.Date(ilin$sites$TIMESTRR)
str(ilin)

## Save object:
save(ilin, file="ilin.rda", compress="xz")
# save(NCSS, file="NCSS.rda", compress="gzip")

save.image("Repository2010.RData")

# end of script;