## Making publication-ready figures

#################
# Functions
#################
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

####################
## Figure 1: map
####################
# load BioTime data
bt <- fread('output/turnover_w_covariates.csv.gz')
trends <- readRDS('temp/trendstemp.rds') # the glmmTMB model fit for all trends

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
    labs(x = 'Longitude (°)', y = 'Latitude (°)', tag = 'A)')

# b) temperature trends on land and freshwater
p2 <- ggplot(temptrends[REALM == 'Terrestrial & Freshwater'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans) +
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
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans) +
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
    labs(tag = 'D)', x = 'Year', y = 'Jaccard turnover') +
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
    scale_x_continuous(trans = signedsqrttrans) +
    labs(tag = 'E)', x = 'Slope', title = 'Jaccard') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# Jtu vs. Jbeta trends
p6 <- ggplot(trendsw, aes(x = Jbeta, y = Jtu)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    labs(tag = 'F)', x = 'Slope', y = 'Slope', title = 'Jaccard turnover vs. total') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# Horn vs. Jtu trends
p7 <- ggplot(trendsw[!is.na(Horn),], aes(x = Jtu, y = Horn)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    labs(tag = 'G)', x = 'Slope', y = 'Slope', title = 'H-M vs. Jaccard total') +
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


#########################
## Figure 2: main effects
#########################
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC
trends <- readRDS('temp/trendstemp.rds') # the lm fit of dissimilarity vs. time, from calc_turnover.Rmd
trends <- trends[duration_group == 'All' & measure == 'Jtu' & rarefyID %in% bt$rarefyID,] # use the slopes that use all data points and all pairs
trends[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trends[, .(disstrend = mean(disstrend, na.rm=TRUE)), by = STUDY_ID] # average by studyID
trends_by_study <- merge(trends_by_study, bt[!duplicated(STUDY_ID), .(STUDY_ID = as.character(STUDY_ID), REALM)]) # add REALM

# predicted slopes from the model
slopespred <- readRDS(here('temp', 'slopes_TsdTTRealm.rds'))
slopespred <- slopespred[round(tempave_metab,1) %in% c(10.1, 30.2),]
slopespred2 <- slopespred
slopespred2[, tempchange := -tempchange_abs]
slopespred <- rbind(slopespred[, .(tempave_metab, REALM, tempchange = tempchange_abs, slope)], slopespred2[, .(tempave_metab, REALM, tempchange, slope)])
slopespredtsign <- readRDS(here('temp', 'slopes_tsign.rds'))
slopespredtsign <- slopespredtsign[round(tempave_metab,1) %in% c(10.1, 30.2),]


# a) across realms
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = REALM, fill = REALM)) +
    geom_density(color = NA, alpha = 0.25) +
    scale_y_sqrt() +
    scale_x_continuous(trans = signedsqrttrans) +
    labs(tag = 'A)', x = 'Slope', title = 'Jaccard') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  

# b) plot of change vs. dT
p2 <- ggplot(slopespred, aes(tempchange, slope, color = factor(round(tempave_metab)), group = tempave_metab)) +
    geom_line() +
    facet_grid(cols = vars(REALM))  +
    labs(tag = 'B)', x = '|dT|') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  


p2 <- ggplot(slopespredtsign, aes(tempchange, slope, color = factor(round(tempave_metab)), group = interaction(tempave_metab))) +
    geom_line() +
    facet_grid(cols = vars(REALM))  +
    labs(tag = 'B)', x = '|dT|') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))  


#########################
## Figure 3: interactions
#########################
newdat <- fread('temp/interactions.csv')

intstoplot <- c('tempave_metab', 'microclim', 'npp', 'human_bowler')
titles <- c('Metabolic temperature', 'Microclimates', 'NPP', 'Human impacts')
realms <- c(NA, NA, NA, 'TerrFresh')
logs <- c(FALSE, FALSE, TRUE, FALSE)
tags <- c('A', 'B', 'C', 'D')
# prep the plots
intplots <- vector('list', length(intstoplot))
for(j in 1:length(intstoplot)){
    subs <- newdat$var == intstoplot[j] & newdat$temptrend < 0 # select one side side
    xvar <- 'temptrend_abs'
    title <- titles[j]
    if(intstoplot[j] %in% c('tsign')){
        subs <- newdat$var == intstoplot[j]
    } 
    if(intstoplot[j] %in% c('thermal_bias')){
        subs <- newdat$var == intstoplot[j]
        xvar <- 'temptrend'
    } 
    if(intstoplot[j] %in% c('human_bowler')){
        subs <- newdat$var == intstoplot[j] & newdat$temptrend < 0 & newdat$REALM2 == realms[j]
    } 
    
    thisplot <- ggplot(newdat[subs, ], 
                       aes_string(x = xvar, y = 'preds', 
                                  group = intstoplot[j], 
                                  color = intstoplot[j])) +
        geom_line() +
        coord_cartesian(ylim = c(-0.05, 0.4)) +
        theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
        labs(title = title, y = 'Slope', x = 'Temperature trend (°C/yr)', tag = tags[j]) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key=element_blank(),
              legend.title=element_blank(),
              axis.text=element_text(size=10),
              axis.title=element_text(size=12))
    
    if(logs[j]){
        intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'log')
    }
    if(!logs[j]){
        intplots[[j]] <- thisplot + scale_color_distiller(palette = "YlGnBu", trans = 'identity')
    }

}

fig2 <- arrangeGrob(grobs = intplots, ncol = 2)

ggsave('figures/fig3.png', fig2, width = 7, height = 6, units = 'in')



##############################
## Figure S1: duration problem
##############################
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
prop <- cors[, .(prop = sum(p < 0.05)/length(p), n = length(p)), by = .(range,n)]
prop[, se := sqrt((prop * (1-prop))/n)]

# make plots of dissimilarity vs. duration with different durations plotted
png(file = 'figures/fig1.png', width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a
bt[rarefyID == '339_1085477', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Jaccard turnover', col = '#00000044', bty = 'l', ylim = c(0,1))]
bt[rarefyID == '339_1085477', abline(lm(Jtu~dY), col = '#a6cee3', lwd = 3)]
mod5 <- bt[rarefyID == '339_1085477' & dY <=5, lm(Jtu~dY)] # calc trendline
preds <- data.table(dY = 1:20)
preds$Jtu5 <- predict(mod5, preds)
preds[dY <=10, lines(dY, Jtu5, col = '#b2df8a', lwd = 3)]

# part a inset
oldpar <- par()
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
prop[n==1000, plot(range, prop, xlab = 'Range of durations', ylab = 'Proportion false positive', ylim = c(0,1))]
prop[n==1000, error.bar(range, prop, lower = se, upper = se, length = 0.02)]
prop[n==10000, points(range, prop, col = 'grey')]
prop[n==10000, error.bar(range, prop, lower = se, upper = se, length = 0.02, col = 'grey')]
abline(h = 0.05, lty = 2, col = 'red')

# part c
modgam <- trends[measure == 'Jtu' & duration_group == 'All', gam(disstrend~s(duration))]
predsgam <- data.table(duration = 1:120)
predsgam[, c('disstrend', 'se') := predict(modgam, newdata = predsgam, se.fit = TRUE)]

trends[measure == 'Jtu' & duration_group == 'All', plot(duration, disstrend, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Jaccard turnover slope', bty = 'l')]
predsgam[, lines(duration, disstrend, col = 'red')]
abline(h = 0, lty = 2)

# part d
plot(-1, -1, xlim = c(0,120), ylim = c(-0.04, 0.04), xlab = 'Duration', ylab = 'Jaccard turnover slope', bty = 'l')
predsgam[, polygon(c(duration, rev(duration)), c(disstrend+se, rev(disstrend - se)), col = '#00000044', border = NA)]
predsgam[, lines(duration, disstrend, col = 'red')]
abline(h = 0, lty = 2)

dev.off()
