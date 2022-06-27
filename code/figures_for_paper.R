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
library(mgcv) # for gam smoother
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

# binomial ci
binomci <- function(x, n){
    out <- as.numeric(binom.test(x, n)$conf.int)
    return(list(out[1], out[2]))
}

### Dataset sizes ---------
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC

nrow(bt)
bt[, .(N = length(unique(rarefyID))), by = REALM]
bt[, length(unique(STUDY_ID))]
bt[, length(unique(rarefyID)), by = taxa_mod2]
bt[, range(duration+1)] # range of years sampled (2 to 119)

#### Figure 1: map --------
# load BioTime data
bt <- fread('output/turnover_w_covariates.csv.gz')
trends <- readRDS('temp/trendstemp.rds') # the slope for all trends

# load sampled temperature trends
temptrends <- fread('output/temperature_trends_sampled.csv.gz')

# make temperature trends data.table
temptrends$type <- 'Global'
temptrends <- temptrends[, .(tempchange, type, REALM)]
bt$type <- 'BioTime'
temptrends <- rbind(temptrends, bt[!is.na(tempchange), .(tempchange, type, REALM)])
temptrends[REALM %in% c('Terrestrial', 'Freshwater'), REALM := 'Terrestrial & Freshwater']

# make table of temporal trends by STUDY_ID
# use the slopes that use all data points and all pairs
trends <- trends[duration_group == 'All' & rarefyID %in% bt$rarefyID,]
trendsw <- dcast(trends, rarefyID ~ measure, value.var = 'disstrend') # wide format for plotting
trendsw[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trendsw[, .(Horn = mean(Horn, na.rm=TRUE), Jbeta = mean(Jbeta), Jtu = mean(Jtu)), by = STUDY_ID]

# correlation across dissimilarity metrics
trendsw[, cor.test(Jtu, Jbeta)]
trendsw[, cor.test(Jtu, Horn)]

# make plot pieces
# a) map
world <- map_data('world')
p1 <- ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white') +
    geom_point(data = bt, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.3, alpha = 0.4, shape = 16)  +
    scale_color_brewer(palette="Set1", name = 'Realm') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8)) +
    labs(x = 'Longitude (°)', y = 'Latitude (°)', tag = 'A)') +
    guides(color = guide_legend(override.aes = list(size=2)))

# b) temperature trends on land and freshwater
p2 <- ggplot(temptrends[REALM == 'Terrestrial & Freshwater'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans, 
                       breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2)) +
    labs(tag = 'B)', x = 'Temperature trend (°C/yr)', title = 'Terrestrial & Freshwater') +
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
    labs(tag = 'C)', x = 'Temperature trend (°C/yr)', title = 'Marine') +
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
    geom_smooth(method = 'glm', method.args = list(family = beta_family(link='logit'))) + # a beta regression
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
    scale_y_sqrt() +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.05, 0, 0.05, 0.2, 0.4)) +
    labs(tag = 'E)', x = 'Turnover\n[proportion species/year]', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# f) Jtu vs. Jbeta trends
p6 <- ggplot(trendsw, aes(x = Jtu, y = Jbeta)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    labs(tag = 'F)', y = 'Total Jaccard', 
         x = 'Turnover\n[prop spp/yr]', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# g) Horn vs. Jtu trends
p7 <- ggplot(trendsw[!is.na(Horn),], aes(x = Jtu, y = Horn)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    labs(tag = 'G)', x = 'Turnover\n[prop spp/yr]', y = 'Horn-Morisita', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

fig1 <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, ncol = 6, 
                   layout_matrix = rbind(c(1,1,1,1,1,1), c(2,2,3,3,4,4), c(5,5,6,6,7,7)),
                   heights=c(unit(0.5, "npc"), unit(0.25, "npc"), unit(0.25, "npc")))

ggsave('figures/fig1.png', fig1, width = 6, height = 6, units = 'in')


### Figure 2: main effects ---------
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC
trends <- fread('output/slope_w_covariates.csv.gz') # the lm fit of dissimilarity vs. time, from assemble_slope_covariates.Rmd
trends <- trends[duration_group == 'All' & measure == 'Jtu' & rarefyID %in% bt$rarefyID,] # use the slopes that use 5 years of data points (to standardize length)
trends[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part B
trends_by_study <- trends[, .(disstrend = mean(disstrend, na.rm=TRUE), temptrend = mean(temptrend, na.rm=TRUE)), by = .(STUDY_ID, REALM)] # average by studyID
trends_by_study[, REALM := factor(REALM, levels = c('Freshwater', 'Terrestrial', 'Marine'))] # re-order for nicer plotting in part A

# average slopes by realm
ave_by_realm <- trends_by_study[, .(disstrend = mean(disstrend), se = sd(disstrend)/sqrt(.N)), by = REALM]
ave_by_realm[, offset := c(-1, 0, 1)] # amount to vertically dodge the lines in part a

# max temptrend by realm, for plotting limits
temptrend_by_realm <- trends[, .(max = max(temptrend, na.rm=TRUE), min = min(temptrend, na.rm=TRUE)), by = REALM]

# predicted slopes from the dT-only model
slopespredTdTT <- readRDS(here('temp', 'slopes_TdTTRealm.rds'))
slopespredTdTT <- slopespredTdTT[round(tempave_metab,1) %in% c(10.1, 25.0),]
slopespredTdTT <- merge(slopespredTdTT, temptrend_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespredTdTT <- slopespredTdTT[tempchange > min & tempchange < max, ] # trim to min & max by realm
slopespredTdTT[, tempave_metab := factor(as.character(round(tempave_metab)), levels = c('25', '10'))] # re-order factor for nicer plotting

# predicted slopes from the main model
slopespred <- readRDS(here('temp', 'slopes_TsdTTRealm.rds'))
slopespred <- slopespred[round(tempave_metab,1) %in% c(10.1, 25.0),]
slopespred2 <- slopespred # make the negative temperature change points
slopespred2[, tempchange := -tempchange_abs]
slopespred <- rbind(slopespred[, .(tempave_metab, REALM, tempchange = tempchange_abs, slope, slope.se)], 
                    slopespred2[, .(tempave_metab, REALM, tempchange, slope, slope.se)]) # merge neg and pos temp change points
slopespred <- merge(slopespred, temptrend_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespred <- slopespred[tempchange > min & tempchange < max, ] # trim to min & max by realm
slopespred[, tempave_metab := factor(as.character(round(tempave_metab)), levels = c('25', '10'))] # re-order factor for nicer plotting


# a) across realms
ht <- 6.3
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = REALM, fill = REALM)) +
    geom_density(color = NA, alpha = 0.25) +
    scale_y_sqrt(breaks = c(0.1, 1, 2, 6)) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.4)) +
    geom_segment(data = ave_by_realm, aes(x=disstrend - 1.96*se, xend = disstrend + 1.96*se, y= ht+offset, yend = ht+offset, color = REALM), alpha = 1) +
    geom_segment(data = ave_by_realm, aes(x = disstrend, y = 0, xend = disstrend, yend = ht+offset, color = REALM), size=0.5, linetype = 'dashed') +
    labs(tag = 'A)', x = 'Turnover [proportion species/year]', y = 'Density', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  
p1 <- addSmallLegend(p1, pointSize = 0.5, spaceLegend = 0.1, textSize = 6)

# b) plot of change vs. dT
p2 <- ggplot() +
    geom_point(data = trends[!is.na(temptrend)], mapping = aes(temptrend, disstrend, size = duration), 
             color='#000000', alpha = 0.1, stroke = 0) +
    geom_line(data = slopespredTdTT, mapping=aes(tempchange, slope, color = tempave_metab, group = tempave_metab), linetype = 'dashed', alpha = 0.5) +
    geom_line(data = slopespred, mapping=aes(tempchange, slope, color = tempave_metab, group = tempave_metab)) +
    geom_ribbon(data = slopespred, alpha = 0.25, color = NA, 
                aes(tempchange, slope, fill = tempave_metab, ymin=slope-slope.se, ymax=slope+slope.se)) +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'B)', x = 'Temperage change [°C/year]', y = 'Turnover [proportion species/year]', 
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
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50))
    
p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)
fig2 <- arrangeGrob(p1, p2, nrow = 2, heights = c(1,2))

ggsave('figures/fig2.png', fig2, width = 6, height = 4, units = 'in')

### Figure 3: interactions ---------
slopes2 <- readRDS(here('temp', 'slopes_interactions2.rds'))
slopes2[, ':='(microclim = as.factor(signif(microclim,2)),
               npp = as.factor(signif(npp,2)),
               seas = as.factor(signif(seas,2)),
               human_bowler = as.factor(signif(human_bowler,2)))] # set as factors for plotting

p1 <- ggplot(slopes2[tempave_metab == 30, ], aes(tempchange_abs, slope_microclim, color = microclim, fill = microclim, group = microclim,
                                                 ymin=slope_microclim-slope_microclim.se,  ymax=slope_microclim+slope_microclim.se)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'A)', x = '|Temperage change| (°C/year)', y = 'Slope', color = 'Microclimate') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  

p2 <- ggplot(slopes2[tempave_metab == 30, ], aes(tempchange_abs, slope_npp, color = npp, fill = npp, group = npp,
                                                 ymin=slope_npp-slope_npp.se,  ymax=slope_npp+slope_npp.se)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'B)', x = '|Temperage change| (°C/year)', y = 'Slope', color = 'NPP         ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  

p3 <- ggplot(slopes2[tempave_metab == 30, ], aes(tempchange_abs, slope_seas, color = seas, fill = seas, group = seas,
                                                 ymin=slope_seas-slope_seas.se,  ymax=slope_seas+slope_seas.se)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'C)', x = '|Temperage change| (°C/year)', y = 'Slope', color = 'Seasonality ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  

p4 <- ggplot(slopes2[tempave_metab == 30, ], aes(tempchange_abs, slope_human, color = human_bowler, fill = human_bowler, group = human_bowler,
                                                 ymin=slope_human-slope_human.se,  ymax=slope_human+slope_human.se)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    facet_grid(cols = vars(REALM)) +
    labs(tag = 'D)', x = '|Temperage change| (°C/year)', y = 'Slope', color = 'Human       ') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  


fig3 <- arrangeGrob(p1, p2, p3, p4, ncol = 1)

ggsave('figures/fig3.png', fig3, width = 4, height = 6, units = 'in')



### Figure S1: duration problem ----------
# load raw BioTime
load(here::here('data', 'biotime_blowes', 'all_pairs_beta.Rdata')) # load rarefied_beta_medians, which has all pairwise dissimilarities
bt <- data.table(rarefied_beta_medians); rm(rarefied_beta_medians)
bt[, year1 := as.numeric(year1)] # not sure why it gets read in as character
bt[, year2 := as.numeric(year2)]
bt[, dY := year2 - year1]
bt[, Horn := 1-Hornsim] # convert similarity to dissimilarity
bt[, Hornsim := NULL]

# load biotime trends
trends <- fread('output/slope_w_covariates.csv.gz')

# load simulations
cors <- fread(here('output', 'simulated_ts.csv.gz'))
corsl <- melt(cors, id.vars = c('n', 'minduration', 'maxduration', 'name', 'range'), measure.vars = c('cor.p', 'cor.cor', 'lm.m', 'glmmwgt.p', 'glmmwgt.beta', 'glmmonegauss.p', 'glmmonegauss.beta', 'glmmonebeta.p', 'glmmonebeta.beta'))
prop <- corsl[variable %in% c('cor.p', 'glmmwgt.p', 'glmmonegauss.p', 'glmmonebeta.p'), 
              .(nsims = sum(!is.na(value)), prop = sum(value < 0.05, na.rm=TRUE)/sum(!is.na(value), na.rm=TRUE)), by = c("range", "n", "variable")]
prop[, c("lower", "upper") := binomci(nsims*prop, nsims), by = .(range, n, variable)]


# make plots of dissimilarity vs. duration with different durations plotted
png(file = 'figures/figS1.png', width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a
bt[rarefyID == '339_1085477', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Jaccard turnover dissimilarity', col = '#00000044', bty = 'l', ylim = c(0,1))]
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

# part b
modgam <- trends[measure == 'Jtu' & duration_group == 'All', gam(disstrend~s(duration))] # gam fit
predsgam <- data.table(duration = 1:120)
predsgam[, c('disstrend', 'se') := predict(modgam, newdata = predsgam, se.fit = TRUE)]

trends[measure == 'Jtu' & duration_group == 'All', plot(duration, disstrend, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Jaccard turnover slope', bty = 'l')]
predsgam[, lines(duration, disstrend, col = 'red')]
abline(h = 0, lty = 2)
mtext('B)', side = 3, line = -0.5, adj = -0.28)

# part c
plot(-1, -1, xlim = c(0,120), ylim = c(-0.04, 0.04), xlab = 'Duration', ylab = 'Jaccard turnover slope', bty = 'l')
predsgam[, polygon(c(duration, rev(duration)), c(disstrend+se, rev(disstrend - se)), col = '#00000044', border = NA)]
predsgam[, lines(duration, disstrend, col = 'red')]
abline(h = 0, lty = 2)
mtext('C)', side = 3, line = -0.5, adj = -0.28)

# part d
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
legend('topleft', legend = c('Pearson correlation', 'Two-stage ME', 'One-stage Gaussian ME', 'One-stage Beta ME'), col = cols, pch = 1, cex=0.5)
mtext('D)', side = 3, line = -0.5, adj = -0.28)

dev.off()

