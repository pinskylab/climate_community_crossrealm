# Making publication-ready figures

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(nlme) # for ME models
library(glmmTMB) # for beta regression
library(maps) # for map
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together
library(RColorBrewer)
library(scales) # for defining signed sqrt axis transformation
library(here)
library(rcompanion) # for CIs on median
source(here('code', 'error.bar.R'))

# produce ggplot-style colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

signedsqrt = function(x) sign(x)*sqrt(abs(x))
signedsq = function(x) sign(x) * x^2
signedsqrttrans <- trans_new(name = 'signedsqrt', transform = signedsqrt, inverse = signedsq)


# from https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
addSmallLegend <- function(myPlot, pointSize = 0.5, textSize = 3, spaceLegend = 0.1) {
    newplot <- myPlot +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize), 
              legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
    return(newplot)
}

# binomial ci (95% by default)
binomci <- function(x, n){
    out <- as.numeric(binom.test(x, n)$conf.int)
    return(list(out[1], out[2]))
}

# given a duration, make a Gaussian white noise timeseries and return the slope
calcslopeGauss <- function(dur){
    x <- 1:dur
    y <- rnorm(dur)
    return(coef(lm(y~x))[2])
}

# convert back from scaled covariates. Assumes scaling is loaded as scalingall object
unscaleme <- function(x.sc, nm){
    if(!(nm %in% scalingall[,var])) stop('nm not found in scalingall')
    if(scalingall[var==nm, log]){
        x <- exp(x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center]) - scalingall[var==nm, plus]
    } else {
        x <- x.sc * scalingall[var == nm, scale] + scalingall[var == nm, center] - scalingall[var==nm, plus]
    }
    return(x)
}


### Dataset sizes ---------
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC
scalingall <- fread('output/turnover_w_covariates_scaling.csv') # covariate scaling data. From assemble_turnover_covariates.Rmd

bt[, .(nyears = length(unique(c(year1, year2)))), by = rarefyID][, sum(nyears)] # number of assemblage composition records
bt[, length(unique(STUDY_ID))] # number of studies
bt[, length(unique(rarefyID))] # number of timeseries
bt[, length(unique(paste(rarefyID_x, rarefyID_y)))] # number of unique locations
bt[, .(N = length(unique(rarefyID))), by = REALM] # numbers of time-series by realm
bt[, length(unique(rarefyID)), by = taxa_mod2] # number of time-series by taxon group
bt[, .(Nts = length(unique(rarefyID))), by = STUDY_ID][Nts >1, .N] # number of studies with >1 rarefyID
bt[, range(duration+1)] # range of years sampled (2 to 119)
bt[unscaleme(tempave.sc, 'tempave.sc') >10 & REALM=='Marine', length(unique(rarefyID))] # number of timeseries >10degC
bt[unscaleme(tempave.sc, 'tempave.sc') >10  & REALM=='Marine', length(unique(rarefyID))]/bt[REALM=='Marine', length(unique(rarefyID))] # proportion of timeseries >10degC

# number of time series
bt[, length(unique(rarefyID))]

# number of studies with negative slopes
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & rarefyID %in% bt$rarefyID & year2 - year1 >2,]
trendsw <- dcast(trends, rarefyID ~ measure, value.var = 'disstrend') # wide format for plotting
trendsw[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trendsw[, .(Horn = mean(Horn, na.rm=TRUE), Jbeta = mean(Jbeta), Jtu = mean(Jtu)), by = STUDY_ID]
trends_by_study[Jtu < 0, .N] # number of Jtu trends < 0
trends_by_study[Jtu < 0, .N/trends_by_study[!is.na(Jtu), .N]] # 15% < 0

# fraction of time series with thermal affinities for most species
cti <- fread(here('output', 'cti_byrarefyID.csv.gz'))
btcti <- merge(bt[!duplicated(rarefyID),], cti, by = 'rarefyID') # number of timeseries
btcti[(nspp_wdata/nspp) >= 0.95, .N] # number of ts for which at least 95% of species have thermal affinity
btcti[(nspp_wdata/nspp) >= 0.95, .N]/nrow(btcti) # fraction of ts for which at least 95% of species have thermal affinity
btcti[(nspp_wdata/nspp) >= 0.75, .N] # number of ts for which at least 75% of species have thermal affinity
btcti[(nspp_wdata/nspp) >= 0.75, .N]/nrow(btcti) # fraction of ts for which at least 75% of species have thermal affinity


### Miscellaneous statistics -----------

# temporal turnover for Swedish birds
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends[duration_group == 'All' & rarefyID =='339_1085477' & measure=='Jtu',]


# median temporal turnover across studies
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & rarefyID %in% bt$rarefyID & year2 - year1 >2,]
trendsw <- dcast(trends, rarefyID ~ measure, value.var = 'disstrend') # wide format for plotting
trendsw[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trendsw[, .(Horn = mean(Horn, na.rm=TRUE), Jbeta = mean(Jbeta), Jtu = mean(Jtu)), by = STUDY_ID]
trends_by_study[, median(Jtu)]
groupwiseMedian(Jtu ~ 1, data = trends_by_study, conf = 0.95, R = 5000, percentile = TRUE, 
                                  bca = FALSE, basic = FALSE, normal = FALSE, wilcox = FALSE, digits = 3)


#### Table 1: AICs --------------
if(!exists('modAllJtu')) modAllJtu <- readRDS(here('temp', 'modAllJtu.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
if(!exists('modRealmAllJtu')) modRealmAllJtu <- readRDS('temp/modRealmAllJtu.rds') # Realm. Fit by code/turnover_GLMM_fit.R
if(!exists('modTaxamod2AllJtu')) modTaxamod2AllJtu <- readRDS('temp/modTaxamod2AllJtu.rds') # Taxon. Fit by code/turnover_GLMM_fit.R
if(!exists('modsdTtsignAllJtu')) modsdTtsignAllJtu <- readRDS(here('temp', 'modsdTtsignAllJtu.rds')) # tsign, tempchange_abs. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R
if(!exists('modsdTRealmtsignAllJtu')) modsdTRealmtsignAllJtu <- readRDS(here('temp', 'modsdTRealmtsignAllJtu.rds')) # tsign:tempchange_abs by realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R
if(!exists('modabsLatsdTabsLatRealmtsignAllJtu')) modabsLatsdTabsLatRealmtsignAllJtu <- readRDS(here('temp', 'modabsLatsdTabsLatRealmtsignAllJtu.rds')) # tsign:tempchange_abs:absLat Fit by code/turnover_vs_temperature_GLMM_fit_modabsLatsdTabsLatRealmtsignAllJtu.R
if(!exists('modrawTsdTTRealmtsignAllJtu')) modrawTsdTTRealmtsignAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsignAllJtu.rds')) # tsign:tempchange_abs:tempave:realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R

# compare sdT amd TsdTT models against null
aics <- AIC(modAllJtu, modRealmAllJtu, modTaxamod2AllJtu, # simple models w/out tempchange
            modsdTtsignAllJtu, modsdTRealmtsignAllJtu, # tsign:tempchange_abs w/out or w/ realm
            modabsLatsdTabsLatRealmtsignAllJtu, modrawTsdTTRealmtsignAllJtu) # add tempave:tempchange_abs
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modAllJtu']
aics

write.csv(aics, here('figures', 'table1.csv'))






#### Figure 1: map and data --------
# load BioTime data
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
bt[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part A
trends <- fread('output/slope.csv.gz') # from calc_turnover.R


# load sampled temperature trends
temptrends <- fread('output/temperature_trends_sampled.csv.gz') # from sample_global_temp.Rmd

# make temperature trends data.table
temptrends$type <- 'Global'
temptrends <- temptrends[, .(tempchange, type, REALM)]
bt$type <- 'BioTime'
temptrends <- rbind(temptrends, bt[!is.na(tempchange), .(tempchange, type, REALM)])
temptrends[REALM %in% c('Terrestrial', 'Freshwater'), REALM := 'Terrestrial & Freshwater']

# make table of temporal trends by STUDY_ID
# use the slopes that use all data points and all pairs
trends <- trends[duration_group == 'All' & rarefyID %in% bt$rarefyID & year2 - year1 >2,]
trendsw <- dcast(trends, rarefyID ~ measure, value.var = 'disstrend') # wide format for plotting
trendsw[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trendsw[, .(Horn = mean(Horn, na.rm=TRUE), Jbeta = mean(Jbeta), Jtu = mean(Jtu)), by = STUDY_ID]

# average temperate change by realm, standard error of the mean, and standard deviation (degC per year)
temptrends[, .(mean = mean(tempchange), se = sd(tempchange)/sqrt(.N), sd = sd(tempchange)), by = REALM]

# range of trends in Fig. 1E
trends_by_study[, range(Jtu)]
trends_by_study[, median(Jtu)]

# make plot pieces
# a) map
world <- map_data('world')
p1 <- ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white', size = 0.1) +
    geom_point(data = bt, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.3, alpha = 0.4, shape = 16)  +
    scale_color_brewer(palette="Dark2", name = 'Realm') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8)) +
    labs(x = 'Longitude (°)', y = 'Latitude (°)', tag = 'A)') +
    guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))

# b) temperature trends on land and freshwater
p2 <- ggplot(temptrends[REALM == 'Terrestrial & Freshwater'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans, 
                       breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2)) +
    scale_fill_brewer(palette="Set2", name = 'Realm') +
    labs(tag = 'B)', x = 'Temperature trend [°C/yr]', title = 'Terrestrial & Freshwater') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = 'none',
          axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# c) temperature trends at sea
p3 <- ggplot(temptrends[REALM == 'Marine'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans, 
                       breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2)) +
    scale_fill_brewer(palette="Set2", name = 'Realm') +
    labs(tag = 'C)', x = 'Temperature trend [°C/yr]', title = 'Marine') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = c(0.8, 0.95),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))
p3 <- addSmallLegend(p3, pointSize = 0.7, spaceLegend = 0.15, textSize = 7)

# d) example of Jtu trend calculation
p4 <- ggplot(bt[rarefyID=='339_1085477', .(dY = year2 - year1, Jtu.sc)], aes(dY, Jtu.sc)) +
    geom_point(alpha = 0.2, size = 0.5, shape = 16) +
    geom_smooth(method = 'glm', method.args = list(family = beta_family(link='logit')), color = '#a6cee3') + # a beta regression
    labs(tag = 'D)', x = 'Temporal distance [years]', y = 'Turnover\n[proportion species]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# e) distribution of Jtu trends
p5 <- ggplot(trends_by_study, aes(x = Jtu)) +
    geom_density(color = NA, alpha = 0.5, fill = 'grey') +
    scale_y_sqrt(breaks = c(0.1,1,3)) +
    geom_vline(xintercept = 0, linetype = 'solid', size = 0.5) +
    geom_vline(xintercept = 0.008, linetype = 'dashed', size = 0.5) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.05, 0, 0.05, 0.2, 0.4)) +
    labs(tag = 'E)', x = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

fig1 <- arrangeGrob(p1, p2, p3, p4, p5, ncol = 4, 
                   layout_matrix = rbind(c(1,1,1,1), c(2,2,3,3), c(4,4,5,5)),
                   heights=c(unit(0.5, "npc"), unit(0.25, "npc"), unit(0.25, "npc")))

ggsave('figures/fig1.png', fig1, width = 6, height = 6, units = 'in')


### Figure 2: main effects ---------
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & measure == 'Jtu',]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, REALM, tempchange)])
trends[, duration := year2 - year1]
trends[, REALM := factor(REALM, levels = c('Marine', 'Terrestrial', 'Freshwater'))] # re-order for nicer plotting in part B
trends_by_study <- trends[duration>2, .(disstrend = mean(disstrend, na.rm=TRUE), tempchange = mean(tempchange, na.rm=TRUE), duration = mean(duration)), by = .(STUDY_ID, REALM)] # average by studyID. Can't use 2-year trends since they assume dissimilarity at y0 is 0.
trends_by_study[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part A

# average slopes by realm
ave_by_realm <- trends_by_study[, .(disstrend = mean(disstrend), se = sd(disstrend)/sqrt(.N), duration = mean(duration), duration.se = sd(duration)/sqrt(.N)), by = REALM]
ave_by_realm[, offset := c(-1, 0, 1)] # amount to vertically dodge the lines in part a
write.csv(ave_by_realm, file='output/ave_by_realm.csv')

# max tempchange by realm, for plotting limits
tempchange_by_realm <- trends[, .(max = max(tempchange, na.rm=TRUE), min = min(tempchange, na.rm=TRUE)), by = REALM]

# predicted slopes from the tsign model (no tempave)
slopespredsdT <- readRDS(here('temp', 'slopes_modsdTRealmtsignAllJtu.rds')) # from pred_modrawXAllJtu.sh
slopespredsdT <- merge(slopespredsdT, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespredsdT <- slopespredsdT[tempchange > min & tempchange < max & !duplicated(cbind(tempchange, REALM)), ] # trim to min & max by realm

# predicted slopes from the tempave interaction model
slopespred <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsign.rds')) # from pred_GLMMmodrawTsdTTRealmtsignAllJtu.R
vals <- slopespred[c(which.min(abs(tempave - 0)), which.min(abs(tempave-25))), unique(tempave)] # show 0 and 25degC (+/-2SD from the mean)
slopespred <- slopespred[tempave %in% vals,]
slopespred <- merge(slopespred, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespred <- slopespred[tempchange > min & tempchange < max, ] # trim to min & max by realm
slopespred[, tempave := factor(as.character(round(tempave)), levels = c('0', '25'))] # re-order factor for nicer plotting


# fastest turnover at highest observed rate of temperature change
slopespredsdT[, max(slope)] # just looking at tempchange
slopespred[, max(slope_realmtsign)] # also considering tempave

# a) across realms
ht <- 6.3
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = REALM, fill = REALM)) +
    geom_density(color = NA, alpha = 0.4) +
    scale_y_sqrt(breaks = c(0.1, 1, 2, 6)) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4)) +
    geom_segment(data = ave_by_realm, aes(x=disstrend - 1.96*se, xend = disstrend + 1.96*se, y= ht+offset, yend = ht+offset, color = REALM), alpha = 1) +
    geom_segment(data = ave_by_realm, aes(x = disstrend, y = 0, xend = disstrend, yend = ht+offset, color = REALM), size=0.5, linetype = 'dashed') +
    labs(tag = 'A)', x = expression('Turnover rate ['~Delta~'Turnover/year]'), y = 'Density', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_fill_brewer(palette = 'Set2') +
    scale_color_brewer(palette = 'Set2')
p1 <- addSmallLegend(p1, pointSize = 0.5, spaceLegend = 0.1, textSize = 6)

# b) plot of change vs. dT
p2 <- ggplot() +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(tempchange, disstrend, size = duration), 
             color='#AAAAAA', alpha = 0.1, stroke = 0) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_line(data = slopespredsdT, mapping=aes(tempchange, slope), size=0.5) +
    geom_ribbon(data = slopespredsdT, alpha = 0.2, color = NA, 
                aes(tempchange, slope,
                    ymin=slope - slope.se, 
                    ymax=slope + slope.se)) +
    geom_line(data = slopespred, mapping=aes(tempchange, slope_realmtsign, color = tempave, group = tempave), linetype = 'dashed') +
    geom_ribbon(data = slopespred, alpha = 0.2, color = NA,
                aes(tempchange, slope_realmtsign, fill = tempave,
                    ymin=slope_realmtsign - slope_realmtsign.se,
                    ymax=slope_realmtsign + slope_realmtsign.se)) +
    scale_color_brewer(palette = 'Dark2') +
    scale_fill_brewer(palette = 'Dark2') +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Average temperature [°C]', 
         color = 'Average temperature [°C]',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(seq(-1,-0.2, by = 0.2), -0.1, 0, 0.1, seq(0.2, 1, by=0.2))) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend
    
p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)
fig2 <- arrangeGrob(p1, p2, nrow = 2, heights = c(1,2))

ggsave('figures/fig2.png', fig2, width = 6, height = 4, units = 'in')


# w/out T predictions
p2noT <- ggplot() +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(tempchange, disstrend, size = duration), 
               color='#AAAAAA', alpha = 0.1, stroke = 0) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_line(data = slopespred, mapping=aes(tempchange, slope_realmtsign, color = tempave, group = tempave), linetype = 'dashed') +
    geom_ribbon(data = slopespred, alpha = 0.2, color = NA,
                aes(tempchange, slope_realmtsign, fill = tempave,
                    ymin=slope_realmtsign - slope_realmtsign.se,
                    ymax=slope_realmtsign + slope_realmtsign.se)) +
    scale_color_brewer(palette = 'Set2') +
    scale_fill_brewer(palette = 'Set2') +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Average temperature [°C]', 
         color = 'Average temperature [°C]',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(seq(-1,-0.2, by = 0.2), -0.1, 0, 0.1, seq(0.2, 1, by=0.2))) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend
p2noT <- addSmallLegend(p2noT, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)
fig2noT <- arrangeGrob(p1, p2noT, nrow = 2, heights = c(1,2))
ggsave('figures/fig2_nopredsT.png', fig2noT, width = 6, height = 4, units = 'in')


# w/out Tempav x Temptrend predictions
p2noTT <- ggplot() +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(tempchange, disstrend, size = duration), 
               color='#AAAAAA', alpha = 0.1, stroke = 0) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_line(data = slopespredsdT, mapping=aes(tempchange, slope), size=0.5) +
    geom_ribbon(data = slopespredsdT, alpha = 0.2, color = NA, 
                aes(tempchange, slope,
                    ymin=slope - slope.se, 
                    ymax=slope + slope.se)) +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Average temperature [°C]', 
         color = 'Average temperature [°C]',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(seq(-1,-0.2, by = 0.2), -0.1, 0, 0.1, seq(0.2, 1, by=0.2))) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend
p2noTT <- addSmallLegend(p2noTT, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)
fig2noTT <- arrangeGrob(p1, p2noTT, nrow = 2, heights = c(1,2))

ggsave('figures/fig2_nopredsTxT.png', fig2noTT, width = 6, height = 4, units = 'in')


### Figure 3: interactions ---------
slopes2 <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignCovariate.rds')) # from code/pred_GLMMmodrawTsdTTRealmCovariateAllJtu.R with argument tsign
slopes2 <- slopes2[tempave == 10, ]
slopes2[, ':='(microclim = as.factor(signif(microclim,2)),
               npp = as.factor(signif(npp,2)),
               seas = as.factor(signif(seas,2)),
               human_bowler = as.factor(signif(human_bowler,2)))] # set as factors for plotting

# max rates by realm and covariate
slopes2[tempave==10 & tempchange %in% c(2,-1.5), .(slope_microclim, slope_microclim.se, slope_human, slope_human.se), 
        by = .(REALM, tempchange, microclim, human_bowler)]

# plots
p1 <- ggplot(slopes2, aes(tempchange, slope_microclim, color = microclim, fill = microclim, group = microclim,
                                                 ymin=slope_microclim-slope_microclim.se,  ymax=slope_microclim+slope_microclim.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'A)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Microclimate') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + # to keep missing factor levels in the plot if we drop them above
    scale_fill_brewer(drop = FALSE, palette = 'Dark2')

p2 <- ggplot(slopes2, aes(tempchange, slope_human, color = human_bowler, fill = human_bowler, group = human_bowler,
                                           ymin=slope_human-slope_human.se,  ymax=slope_human+slope_human.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Human       ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + 
    scale_fill_brewer(drop = FALSE, palette = 'Dark2')

fig3 <- arrangeGrob(p1, p2, ncol = 1)

ggsave('figures/fig3.png', fig3, width = 4, height = 4, units = 'in')


# only low values
slopes2low <- slopes2[human_bowler == 0.055,] # to manually make a plot with only the low factor levels
# slopes2 <- slopes2[human_bowler == 10,] # to manually make a plot with only the high factor levels. Also have to add , fill = '#00BFC4' to geom_ribbon and color= '#00BFC4' to geom_line

p1low <- ggplot(slopes2low, aes(tempchange, slope_microclim, color = microclim, fill = microclim, group = microclim,
                                           ymin=slope_microclim-slope_microclim.se,  ymax=slope_microclim+slope_microclim.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'A)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Microclimate') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + # to keep missing factor levels in the plot if we drop them above
    scale_fill_brewer(drop = FALSE, palette = 'Dark2') +
    ylim(-0.044, 0.102)

p2low <- ggplot(slopes2low, aes(tempchange, slope_human, color = human_bowler, fill = human_bowler, group = human_bowler,
                                           ymin=slope_human-slope_human.se,  ymax=slope_human+slope_human.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Human       ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + 
    scale_fill_brewer(drop = FALSE, palette = 'Dark2') +
    ylim(-0.052, 0.094)

fig3low <- arrangeGrob(p1low, p2low, ncol = 1)
ggsave('figures/fig3_onlylow.png', fig3low, width = 4, height = 4, units = 'in')


# only high values
slopes2high <- slopes2[human_bowler == 10,] # to manually make a plot with only the high factor levels. 

p1high <- ggplot(slopes2high, aes(tempchange, slope_microclim, color = microclim, fill = microclim, group = microclim,
                                                 ymin=slope_microclim-slope_microclim.se,  ymax=slope_microclim+slope_microclim.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'A)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Microclimate') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + # to keep missing factor levels in the plot if we drop them above
    scale_fill_brewer(drop = FALSE, palette = 'Dark2') +
    ylim(-0.044, 0.102)

p2high <- ggplot(slopes2high, aes(tempchange, slope_human, color = human_bowler, fill = human_bowler, group = human_bowler,
                                                 ymin=slope_human-slope_human.se,  ymax=slope_human+slope_human.se)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Human       ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_color_brewer(drop = FALSE, palette = 'Dark2') + 
    scale_fill_brewer(drop = FALSE, palette = 'Dark2') +
    ylim(-0.052, 0.094)

fig3high <- arrangeGrob(p1high, p2high, ncol = 1)
ggsave('figures/fig3_onlyhigh.png', fig3high, width = 4, height = 4, units = 'in')


#### Table S1: random effects for main model --------------
if(!exists('modrawTsdTTRealmtsignAllJtu')) modrawTsdTTRealmtsignAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsignAllJtu.rds')) # tsign:tempchange_abs:tempave:realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R
if(!exists('sum_modrawTsdTTRealmtsignAllJtu')) sum_modrawTsdTTRealmtsignAllJtu <- summary(modrawTsdTTRealmtsignAllJtu)
capture.output(print(sum_modrawTsdTTRealmtsignAllJtu$varcor), file = 'figures/tableS1.txt')

#### Table S2: fixed effects for main model --------------
if(!exists('modrawTsdTTRealmtsignAllJtu')) modrawTsdTTRealmtsignAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsignAllJtu.rds')) # tsign:tempchange_abs:tempave:realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R
if(!exists('sum_modrawTsdTTRealmtsignAllJtu')) sum_modrawTsdTTRealmtsignAllJtu <- summary(modrawTsdTTRealmtsignAllJtu)
out <- as.data.frame(sum_modrawTsdTTRealmtsignAllJtu$coefficients$cond)

# get term names
out$term <- gsub('Terrestrial|Marine|Freshwater|1|-1', '', rownames(out))
out$term <- gsub('duration', 'Years', out$term)
out$term <- gsub('REALM', 'Realm', out$term)
out$term <- gsub('tsign', 'sign', out$term)
out$term <- gsub('tempchange_abs.sc', 'Ttrend', out$term)
out$term <- gsub('tempave.sc', 'Tave', out$term)
out$term[out$term == 'Years:Realm'] <- 'Realm ✕ Years'
out$term[out$term == 'Years:Realm:sign:Ttrend'] <- 'sign ✕ |Ttrend| ✕ Realm ✕ Years'
out$term[out$term == 'Years:Realm:sign:Tave'] <- 'sign ✕ Tave ✕ Realm ✕ Years'
out$term[out$term == 'Years:Realm:sign:Ttrend:Tave'] <- 'sign ✕ |Ttrend| ✕ Tave ✕ Realm ✕ Years'

# get realm
out$realm <- rownames(out)
out$realm[grepl('Freshwater', out$realm)] <- 'Freshwater'
out$realm[out$term == 'Years'] <- 'Freshwater' # the baseline Year term is for Freshwater
out$term[out$term == 'Years'] <- 'Years:Realm'
out$realm[grepl('Marine', out$realm)] <- 'Marine'
out$realm[grepl('Terrestrial', out$realm)] <- 'Terrestrial'
out$realm[!grepl('Terrestrial|Marine|Freshwater', out$realm)] <- ''

# get tsign
out$tsign <- rownames(out)
out$tsign[grepl('tsign-1', out$tsign)] <- 'cooling'
out$tsign[grepl('tsign1', out$tsign)] <- 'warming'
out$tsign[!grepl('cooling|warming', out$tsign)] <- ''

# reorder columns
out <- out[,c('term', 'realm', 'tsign', 'Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')]

# write out
write.csv(format(out, digits=2), file = 'figures/tableS2.csv', row.names=FALSE)


#### Table S3: Horn AICs --------------
modAllHorn <- readRDS(here('temp', 'modAllHorn.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
modRealmAllHorn <- readRDS('temp/modRealmAllHorn.rds') # Realm. Fit by code/turnover_GLMM_fit.R
modTaxamod2AllHorn <- readRDS('temp/modTaxamod2AllHorn.rds') # Taxon. Fit by code/turnover_GLMM_fit.R
modsdTtsignAllHorn <- readRDS(here('temp', 'modsdTtsignAllHorn.rds')) # tsign, tempchange_abs. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllHorn.R
modsdTRealmtsignAllHorn <- readRDS(here('temp', 'modsdTRealmtsignAllHorn.rds')) # tsign:tempchange by realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllHorn.R
modrawTsdTTRealmtsignAllHorn <- readRDS(here('temp','modrawTsdTTRealmtsignAllHorn.rds')) # adds tsign to tempave:tempchange:realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllHorn.R

# compare sdT amd TsdTT models against null
aics <- AIC(modAllHorn, modRealmAllHorn, modTaxamod2AllHorn, # simple models w/out tempchange
            modsdTtsignAllHorn, modsdTRealmtsignAllHorn, # tsign:tempchange_abs, also by realm
            modrawTsdTTRealmtsignAllHorn) # add tempave:tempchange_abs
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modAllHorn']
aics

write.csv(aics, here('figures', 'tableS3.csv'))


#### Table S4: covariate AICs --------------
# load models for AICs
if(!exists('modAllJtu')) modAllJtu <- readRDS(here('temp', 'modAllJtu.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
if(!exists('modrawTsdTTRealmtsignAllJtu')) modrawTsdTTRealmtsignAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsignAllJtu.rds')) # adds tsign to tempave:tempchange:realm. Fit by code/turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmtsignAllJtu.R
if(!exists('modrawTsdTTRealmtsignmicroclimAllJtu')) modrawTsdTTRealmtsignmicroclimAllJtu <- readRDS('temp/modrawTsdTTRealmtsignmicroclimAllJtu.rds') # has microclimates. Fit by turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmmicroclimAllJtu.R.
if(!exists('modrawTsdTTRealmtsignhumanAllJtu')) modrawTsdTTRealmtsignhumanAllJtu <- readRDS('temp/modrawTsdTTRealmtsignhumanAllJtu.rds') # has human impact. Fit by turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmhumanAllJtu.R

# compare covariate models against null
aics <- AIC(modAllJtu, modrawTsdTTRealmtsignAllJtu, modrawTsdTTRealmtsignmicroclimAllJtu,
            modrawTsdTTRealmtsignhumanAllJtu) 
aics$dAICTsdTT <- aics$AIC - aics$AIC[rownames(aics)=='modrawTsdTTRealmtsignAllJtu']
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modAllJtu']
aics

write.csv(aics, here('figures', 'tableS4.csv'))




#### Table S5: thermal_bias AICs --------------
# load models
modrawTsdTTRealmtsignAllJtu_thermal_biasdata <- readRDS(here('temp','modrawTsdTTRealmtsignAllJtu_thermal_biasdata.rds')) # has thermal bias:tempchange_abs:tsign. Fit by turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmthermal_biasAllJtu.R
modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata <- readRDS(here('temp','modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata.rds')) # has thermal bias:tempchange_abs:tsign. Fit by turnover_vs_temperature_GLMM_fit_modrawTsdTTRealmthermal_biasAllJtu.R

# compare thermal_bias against tsign model
aics <- AIC(modrawTsdTTRealmtsignAllJtu_thermal_biasdata, modrawTsdTTRealmthermal_biassdTAllJtu_thermal_biasdata) 
aics$dAIC <- aics$AIC - min(aics$AIC)
aics

write.csv(aics, here('figures', 'tableS5.csv'))







#### Table S6: Dispersion estimates ----------------
modnms <- c('modAllJtu.rds', 'modRealmAllJtu.rds', 
          'modTaxamod2AllJtu.rds', 'modsdTtsignAllJtu.rds', 
          'modsdTRealmtsignAllJtu.rds', 'modrawTsdTTRealmtsignAllJtu.rds', 
          'modAllHorn.rds', 'modRealmAllHorn.rds', 
          'modTaxamod2AllHorn.rds', 'modsdTtsignAllHorn.rds', 
          'modsdTRealmtsignAllHorn.rds', 'modrawTsdTTRealmtsignAllHorn.rds')
out <- data.frame(modnms = gsub('.rds', '', modnms), freshwater = numeric(12), marine = numeric(12), terrestrial = numeric(12))
for(i in 1:length(modnms)){ # a bit slow to load each model
    cat(i)
    mod <- readRDS(here('temp', modnms[i]))
    out[i, 2:4] <- fixef(mod)$disp
}

out[,2:4] <- signif(out[,2:4], 3)
write.csv(out, here('figures', 'tableS6.csv'))



### Figure S1: time-series info----------
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
btts <- bt[, .(year1 = min(year1), year2 = max(year2), nsamp = length(unique(c(year1, year2)))), by = .(rarefyID, STUDY_ID)] # summarize by time-series

labpos <- -0.3 # horizontal position for the subfigure label

png(file = 'figures/figS1.png', width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(2.5, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a: start dates
btts[, hist(year1, main = '', xlab = '', col = 'grey')]
mtext('Start year', side = 1, line = 1.5, cex=0.8)
mtext('A)', side = 3, line = -0.5, adj = labpos)

# part b: end dates
btts[, hist(year2, main = '', xlab = '', col = 'grey')]
mtext('End year', side = 1, line = 1.5, cex=0.8)
mtext('B)', side = 3, line = -0.5, adj = labpos)

# part c: durations
btts[, hist(year2-year1+1, main = '', xlab = '', col = 'grey')]
mtext('Number of years', side = 1, line = 1.5, cex=0.8)
mtext('C)', side = 3, line = -0.5, adj = labpos)

# part d: number of samples
btts[, hist(nsamp, main = '', xlab = '', col = 'grey')]
mtext('Number of samples', side = 1, line = 1.5, cex=0.8)
mtext('D)', side = 3, line = -0.5, adj = labpos)

dev.off()



### Figure S2: duration problem ----------
# load raw BioTime
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians, which has all pairwise dissimilarities
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, dY := year2 - year1]
bt[, Horn := 1-Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]

# load biotime trends
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends[, duration := year2 - year1]
trends <- trends[measure == 'Jtu' & duration>2, ] # trim to Jaccard turnover

# load temperature slopes
tempchange <- fread('output/turnover_w_covariates.csv.gz') # covariate data. From assemble_turnover_covariates.Rmd
tempchange <- tempchange[, .(tempchange.sc = mean(tempchange.sc, na.rm=TRUE), duration = max(duration)), by = rarefyID] # summarize by rarefyID
tempchange[, tempchange := tempchange.sc * scalingall[var == 'tempchange.sc', scale] + scalingall[var == 'tempchange.sc', center]]

# load simulations
cors <- fread(here('output', 'simulated_ts.csv.gz')) # created by duration_sim.R
corsl <- melt(cors, id.vars = c('n', 'minduration', 'maxduration', 'name', 'range'), measure.vars = c('cor.p', 'cor.cor', 'lm.m', 'glmmwgt.p', 'glmmwgt.beta', 'glmmonegauss.p', 'glmmonegauss.beta', 'glmmonebeta.p', 'glmmonebeta.beta'))
prop <- corsl[variable %in% c('cor.p', 'glmmwgt.p', 'glmmonegauss.p', 'glmmonebeta.p'), 
              .(nsims = sum(!is.na(value)), prop = sum(value < 0.05, na.rm=TRUE)/sum(!is.na(value), na.rm=TRUE)), by = c("range", "n", "variable")]
prop[, c("lower", "upper") := binomci(nsims*prop, nsims), by = .(range, n, variable)]

# make slopes of Gaussian white noise timeseries
set.seed(10)
trends[, gauss.slope := calcslopeGauss(duration), by = rarefyID]

# make mean predictions of turnover rate and temperature change rate
modloess <- trends[, loess(disstrend~duration)] # loess fit
predsloess <- data.table(duration = 2:118)
predsloess[, c('disstrend', 'se') := predict(modloess, newdata = predsloess, se.fit = TRUE)]

modloessgauss <- trends[, loess(gauss.slope~duration)] # loess fit
predsloess[, c('gauss.slope', 'gauss.se') := predict(modloessgauss, newdata = predsloess, se.fit = TRUE)]


modloesstemp <- tempchange[, loess(tempchange~duration)] # loess fit
predsloesstemp <- data.table(duration = 2:118)
predsloesstemp[, c('tempchange', 'se') := predict(modloesstemp, newdata = predsloesstemp, se.fit = TRUE)]


# make plots of dissimilarity vs. duration with different durations plotted
png(file = 'figures/figS2.png', width = 6.5, height = 5, units = 'in', res = 300)
par(mfrow=c(2,3), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a
bt[rarefyID == '339_1085477', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Turnover [proportion of species]', col = '#00000044', bty = 'l', ylim = c(0,1))]
bt[rarefyID == '339_1085477', abline(lm(Jtu~dY), col = '#a6cee3', lwd = 3)]
mod5 <- bt[rarefyID == '339_1085477' & dY <=5, lm(Jtu~dY)] # calc trendline
preds <- data.table(dY = 1:20)
preds$Jtu5 <- predict(mod5, preds)
preds[dY <=10, lines(dY, Jtu5, col = '#b2df8a', lwd = 3)]
mtext('A)', side = 3, line = -0.5, adj = -0.28)

# part a inset
oldpar <- par(no.readonly=TRUE)
par(fig = c(0.05,0.3, 0.8, 1), new = T, mgp = c(0.7, 0.12, 0), cex.lab = 0.7, cex.axis = 0.5, tcl = -0.1)
plot(-1, -1, xlim=c(0,20), ylim=c(0,1), xlab = 'Temporal difference', ylab = 'Tturnover', bty = 'l')
abline(h = 1, lty= 2)
segments(0,0,5,1, col = '#b2df8a', lwd = 3)
segments(5,0,5,1, lty = 2)
segments(0,0,20,1, col = '#a6cee3', lwd = 3)
segments(20,0,20,1, lty = 2)

par(oldpar) # go back to original figure settings
par(mfg = c(1,2)) # start with top-right

# part b: turnover by duration
trends[, plot(duration, disstrend, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), bty = 'l')]
abline(h = 0, lty = 2)
predsloess[, lines(duration, disstrend, col = 'red')]
mtext('B)', side = 3, line = -0.5, adj = -0.28)


# part c: tempchange by duration
tempchange[, plot(duration, tempchange, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Temperature trend [°C/yr]', bty = 'l')]
abline(h = 0, lty = 2)
predsloesstemp[, lines(duration, tempchange, col = 'red')]
mtext('C)', side = 3, line = -0.5, adj = -0.28)


# part d: Gaussian white noise slope by duration
trends[, plot(duration, gauss.slope, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Slope of Gaussian white noise', bty = 'l')]
abline(h = 0, lty = 2)
predsloess[, lines(duration, gauss.slope, col = 'red')]
mtext('D)', side = 3, line = -0.5, adj = -0.28)


# part e: type I error
cols <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c')
dg = c(-3, -1, 1, 3)
prop[variable == 'cor.p', plot(range+dg[1], prop, xlab = 'Range of durations', ylab = 'Proportion false positive', xlim=c(0,102), ylim = c(0,1), col = cols[1], type = 'o', bty = 'l')]
prop[variable == 'cor.p', error.bar(range+dg[1], prop, lower = prop-lower, upper = upper-prop, length = 0.02, col = cols[1])]
prop[variable == 'glmmwgt.p', points(range+dg[2], prop, xlab = 'Range of durations', ylab = 'Proportion false positive', ylim = c(0,1), col = cols[2], type = 'o')]
prop[variable == 'glmmwgt.p', error.bar(range+dg[2], prop, lower = prop-lower, upper = upper-prop, length = 0.02, col = cols[2])]
prop[variable == 'glmmonegauss.p', points(range+dg[3], prop, xlab = 'Range of durations', ylab = 'Proportion false positive', ylim = c(0,1), col = cols[3], type = 'o')]
prop[variable == 'glmmonegauss.p', error.bar(range+dg[3], prop, lower = prop-lower, upper = upper-prop, length = 0.02, col = cols[3])]
prop[variable == 'glmmonebeta.p', points(range+dg[4], prop, xlab = 'Range of durations', ylab = 'Proportion false positive', ylim = c(0,1), col = cols[4], type = 'o')]
prop[variable == 'glmmonebeta.p', error.bar(range+dg[4], prop, lower = prop-lower, upper = upper-prop, length = 0.02, col = cols[4])]
abline(h = 0.05, lty = 2, col = 'red')
legend('topleft', legend = c('Pearson correlation', 'Meta-analysis', 'One-stage Gaussian ME', 'One-stage beta ME'), col = cols, pch = 1, cex=0.5)
mtext('E)', side = 3, line = -0.5, adj = -0.28)


# part f: an example negative slope
bt[rarefyID == '213_435199', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Turnover [proportion of species]', col = '#00000044', bty = 'l', 
                                  ylim = c(0,1))]
mod <- bt[rarefyID == '213_435199', lm(Jtu~dY)] # calc trendline
preds <- data.table(dY = 1:35)
preds[, c('Jtu', 'se', 'df', 'residual.scale') := predict(mod, preds, se.fit=TRUE)]
preds[, polygon(x = c(dY, rev(dY)), y= c(Jtu+se, rev(Jtu-se)), col = '#88888855', border = NA)]
preds[, lines(dY, Jtu, col = '#a6cee3', lwd = 3)]
mtext('F)', side = 3, line = -0.5, adj = -0.28)

dev.off()


### Figure S3: turnover by taxon ----------
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # covariate data
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & measure == 'Jtu',] # trim to those we use
trends[, duration := year2 - year1]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, taxa_mod2)])
trends[taxa_mod2 == 'All', taxa_mod2 := 'Multiple taxa'] # rename a level so more intuitive
trends_by_study <- trends[duration>2, .(disstrend = mean(disstrend, na.rm=TRUE), duration = mean(duration)), by = .(STUDY_ID, taxa_mod2)] # average by studyID

# average slopes by taxon
ave_by_taxon <- trends_by_study[, .(disstrend = mean(disstrend), se = sd(disstrend)/sqrt(.N), duration = mean(duration), duration.se = sd(duration)/sqrt(.N)), by = taxa_mod2][order(taxa_mod2),]
ave_by_taxon[, offset := seq(1, -1, length.out = 9)] # amount to vertically dodge the lines in part a
write.csv(ave_by_taxon, file='output/ave_by_taxon.csv')


# plot
ht <- 6
pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d") # 15-level colorblind-friendly palette from https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = taxa_mod2, fill = taxa_mod2)) +
    geom_density(color = NA, alpha = 0.25) +
    scale_y_sqrt(limits = c(0,7)) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4)) +
    geom_segment(data = ave_by_taxon, aes(x=disstrend - 1.96*se, xend = disstrend + 1.96*se, y= ht+offset, yend = ht+offset, color = taxa_mod2), alpha = 1) +
    geom_segment(data = ave_by_taxon, aes(x = disstrend, y = 0, xend = disstrend, yend = ht+offset, color = taxa_mod2), size=0.5, linetype = 'dashed') +
    labs(x = expression('Turnover rate ['~Delta~'Turnover/year]'), y = 'Density', title = '', fill = 'Taxon', color = 'Taxon') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal)
p1 <- addSmallLegend(p1, pointSize = 0.8, spaceLegend = 0.2, textSize = 8)
ggsave('figures/figS3.png', p1, width = 6, height = 4, units = 'in')



### Figure S4: T trend x Ave T interaction ---------

# read in slopes
slopesTsdTTRealmtsignJtu <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsign.rds')) # made by pred_GLMMmodrawTsdTTRealmtsignAllJtu.R

# plot
p1 <- ggplot(slopesTsdTTRealmtsignJtu, aes(tempchange, tempave, z = slope_realmtsign)) +
    geom_raster(aes(fill = slope_realmtsign)) +
    labs(x = 'Temperature trend (°C/yr)', y = 'Average Temperature (°C)') +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0, name = 'Turnover rate') +
    facet_grid(cols = vars(REALM)) +
    theme(axis.text = element_text(size = 10), 
          axis.title = element_text(size = 12),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -5, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 8),
          legend.title=element_text(size= 12),
          legend.title.align = 1)

ggsave('figures/figS4.png', p1, width = 6, height = 3, units = 'in')
    


### Figure S5: Horn main effects ---------
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & measure == 'Horn',]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, REALM, tempchange)])
trends[, duration := year2 - year1]
trends[, REALM := factor(REALM, levels = c('Marine', 'Terrestrial', 'Freshwater'))] # re-order for nicer plotting in part B
trends_by_study <- trends[duration>2, .(disstrend = mean(disstrend, na.rm=TRUE), tempchange = mean(tempchange, na.rm=TRUE), duration = mean(duration)), by = .(STUDY_ID, REALM)] # average by studyID. Can't use 2-year trends since they assume dissimilarity at y0 is 0.
trends_by_study[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part A

# average slopes by realm
ave_by_realm <- trends_by_study[, .(disstrend = mean(disstrend), se = sd(disstrend)/sqrt(.N), duration = mean(duration), duration.se = sd(duration)/sqrt(.N)), by = REALM]
ave_by_realm[, offset := c(-1, 0, 1)] # amount to vertically dodge the lines in part a

# max tempchange by realm, for plotting limits
tempchange_by_realm <- trends[, .(max = max(tempchange, na.rm=TRUE), min = min(tempchange, na.rm=TRUE)), by = REALM]

# predicted slopes from the tsign model (no tempave)
slopespredsdT <- readRDS(here('temp', 'slopes_modsdTRealmtsignAllHorn.rds')) # from pred_modrawXAllHorn.sh/.R
slopespredsdT <- merge(slopespredsdT, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespredsdT <- slopespredsdT[tempchange > min & tempchange < max & !duplicated(cbind(tempchange, REALM)), ] # trim to min & max by realm

# predicted slopes from the tempave interaction model
slopespred <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsignHorn.rds')) # from pred_GLMMmodrawTsdTTRealmtsignAllHorn.R
vals <- slopespred[c(which.min(abs(tempave - 0)), which.min(abs(tempave-25))), unique(tempave)] # show 0 and 25degC (+/-2SD from the mean)
slopespred <- slopespred[tempave %in% vals,]
slopespred <- merge(slopespred, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespred <- slopespred[tempchange > min & tempchange < max, ] # trim to min & max by realm
slopespred[, tempave := factor(as.character(round(tempave)), levels = c('25', '0'))] # re-order factor for nicer plotting



# a) across realms
ht <- 6.3
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = REALM, fill = REALM)) +
    geom_density(color = NA, alpha = 0.4) +
    scale_y_sqrt(breaks = c(0.1, 1, 2, 6)) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4)) +
    geom_segment(data = ave_by_realm, aes(x=disstrend - 1.96*se, xend = disstrend + 1.96*se, y= ht+offset, yend = ht+offset, color = REALM), alpha = 1) +
    geom_segment(data = ave_by_realm, aes(x = disstrend, y = 0, xend = disstrend, yend = ht+offset, color = REALM), size=0.5, linetype = 'dashed') +
    labs(tag = 'A)', x = expression('Turnover rate ['~Delta~'Turnover/year]'), y = 'Density', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_fill_brewer(palette = 'Set2') +
    scale_color_brewer(palette = 'Set2')
p1 <- addSmallLegend(p1, pointSize = 0.5, spaceLegend = 0.1, textSize = 6)

# b) plot of change vs. dT
p2 <- ggplot() +
    #geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(tempchange, disstrend, size = duration), 
               color='#AAAAAA', alpha = 0.1, stroke = 0) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_line(data = slopespredsdT, mapping=aes(tempchange, slope), size=1) +
    geom_ribbon(data = slopespredsdT, alpha = 0.5, color = NA, 
                aes(tempchange, slope,
                    ymin=slope - slope.se, 
                    ymax=slope + slope.se)) +
    geom_line(data = slopespred, mapping=aes(tempchange, slope_realmtsign, color = tempave, group = tempave), linetype = 'dashed') +
    geom_ribbon(data = slopespred, alpha = 0.2, color = NA,
                aes(tempchange, slope_realmtsign, fill = tempave,
                    ymin=slope_realmtsign - slope_realmtsign.se,
                    ymax=slope_realmtsign + slope_realmtsign.se)) +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'B)', x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Average temperature [°C]', 
         color = 'Average temperature [°C]',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(seq(-1,-0.2, by = 0.2), -0.1, 0, 0.1, seq(0.2, 1, by=0.2))) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend

p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)
figS5 <- arrangeGrob(p1, p2, nrow = 2, heights = c(1,2))

ggsave('figures/figS5.png', figS5, width = 6, height = 4, units = 'in')



### Figure S6: thermal bias ---------
slopesTB <- readRDS(here('temp', 'slopes_rawTsdTTRealmthermal_bias.rds')) # from pred_GLMMmodrawTsdTTRealmthermal_biasAllJtu.R
slopesTB[, ':='(thermal_bias = as.factor(signif(thermal_bias,2)))] # set as factors for plotting

p1 <- ggplot(slopesTB[tempave == 30, ], aes(tempchange, slope_thermal_biassdT, color = thermal_bias, fill = thermal_bias, group = thermal_bias,
                                            ymin=slope_thermal_biassdT-slope_thermal_biassdT.se,  ymax=slope_thermal_biassdT+slope_thermal_biassdT.se)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(x = 'Temperature trend [°C/year]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), color = 'Thermal bias') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  


ggsave('figures/figS6.png', p1, width = 6, height = 3, units = 'in')
