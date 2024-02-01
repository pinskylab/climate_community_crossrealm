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
source(here('code', 'util.R'))


### Dataset sizes and descriptive stats ---------
bt <- fread('output/turnover_w_covariates.csv.gz') # the timeseries that pass QA/QC. from assemble_turnover_covariates.Rmd
scalingall <- fread('output/turnover_w_covariates_scaling.csv') # covariate scaling data. From assemble_turnover_covariates.Rmd

bt[, .(nyears = length(unique(c(year1, year2)))), by = rarefyID][, sum(nyears)] # number of assemblage composition records
bt[, length(unique(STUDY_ID))] # number of studies
bt[, length(unique(rarefyID))] # number of timeseries
bt[, length(unique(paste(rarefyID_x, rarefyID_y)))] # number of unique locations
bt[, .(N = length(unique(rarefyID))), by = REALM] # numbers of time-series by realm
bt[, length(unique(rarefyID)), by = taxa_mod2] # number of time-series by taxon group
bt[, .(Nts = length(unique(rarefyID))), by = STUDY_ID][Nts >1, .N] # number of studies with >1 rarefyID
bt[, range(duration+1)] # range of years sampled (2 to 119)
bt[, median(duration+1)] # median time series length
bt[unscaleme(tempave.sc, 'tempave.sc') >10 & REALM=='Marine', length(unique(rarefyID))] # number of timeseries >10degC
bt[unscaleme(tempave.sc, 'tempave.sc') >10  & REALM=='Marine', length(unique(rarefyID))]/bt[REALM=='Marine', length(unique(rarefyID))] # proportion of timeseries >10degC

# temporal turnover for Swedish birds
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends[duration_group == 'All' & rarefyID =='339_1085477' & measure=='Jtu',]


# median and range of temporal turnover across studies
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & rarefyID %in% bt$rarefyID & year2 - year1 >2,]
trendsw <- dcast(trends, rarefyID ~ measure, value.var = 'disstrend') # wide format for plotting
trendsw[, STUDY_ID := vapply(strsplit(rarefyID,"_"), `[`, 1, FUN.VALUE=character(1))] # extract STUDY_ID from rarefyID
trends_by_study <- trendsw[, .(Horn = mean(Horn, na.rm=TRUE), Jbeta = mean(Jbeta), Jtu = mean(Jtu)), by = STUDY_ID]
trends_by_study[, median(Jtu)]
groupwiseMedian(Jtu ~ 1, data = trends_by_study, conf = 0.95, R = 5000, percentile = TRUE, 
                                  bca = FALSE, basic = FALSE, normal = FALSE, wilcox = FALSE, digits = 3)

trends_by_study[, range(Jtu)]



### Likelihood ratio tests among models ----------------
if(!exists('modRealmInitAllJtu')) modRealmInitAllJtu <- readRDS(here('temp', 'modRealmInitAllJtu.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
if(!exists('modsdTRealmtsigninitAllJtu')) modsdTRealmtsigninitAllJtu <- readRDS(here('temp', 'modsdTRealmtsigninitAllJtu.rds')) # Tchange model Fit by code/turnover_GLMM_fit.R
if(!exists('modrawTsdTTRealmtsigninitAllJtu')) modrawTsdTTRealmtsigninitAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R
if(!exists('modrawTsdTTRealmtsignmicroclimInitAllJtu')) modrawTsdTTRealmtsignmicroclimInitAllJtu <- readRDS('temp/modrawTsdTTRealmtsignmicroclimInitAllJtu.rds') # Microclimate model. Fit by code/turnover_GLMM_fit.R
if(!exists('modrawTsdTTRealmtsignhumanInitAllJtu')) modrawTsdTTRealmtsignhumanInitAllJtu <- readRDS('temp/modrawTsdTTRealmtsignhumanInitAllJtu.rds') # Human impact model. Fit by code/turnover_GLMM_fit.R

anova(modRealmInitAllJtu, modsdTRealmtsigninitAllJtu) # Realm vs. Tchange model
anova(modsdTRealmtsigninitAllJtu, modrawTsdTTRealmtsigninitAllJtu) # Tchange vs. Tchange x Tave model
anova(modrawTsdTTRealmtsigninitAllJtu, modrawTsdTTRealmtsignmicroclimInitAllJtu) # Tchange x Tave model vs. microclimate
anova(modrawTsdTTRealmtsigninitAllJtu, modrawTsdTTRealmtsignhumanInitAllJtu) # Tchange x Tave model vs. human



### Table 1: AICs --------------

# use models
modInit <- readRDS(here('temp', 'modOBInitAllJtu.rds')) # Null with only duration. Fit by code/turnover_?
modRealm <- readRDS('temp/modOBRRealmInitAllJtu.rds') # Realm. Fit by code/turnover_GLMM_?
modTaxamod2 <- readRDS('temp/modOBTTaxamod2InitAllJtu.rds') # Taxon. Fit by code/turnover_GLMM_?
modTchange <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllJtu.rds')) # Tchange model. Fit by code/turnover_GLMM_fit.R
modLat <- readRDS(here('temp', 'modOBabsLatsdTabsLatMERtsRealmtsignInitAllJtu.rds'))
modTchangeTave <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R


# compare Tchange amd Tchange x Tave models against null
aics <- AIC(modInit, modRealm, modTaxamod2, # simple models w/out Tchange
            modTchange, # Tchange
            modLat, # latitude
            modTchangeTave) # add Tchange x Tave
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInit']
aics

write.csv(aics, here('figures', 'table1.csv'))






### Figure 1: map and data --------
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

# average temperature change by realm, standard error of the mean, and standard deviation (degC per year)
temptrends[, .(mean = mean(tempchange), se = sd(tempchange)/sqrt(.N), sd = sd(tempchange)), by = REALM]

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
          plot.tag=element_text(face='bold'),
          axis.title=element_text(size=8)) +
    labs(x = 'Longitude (°)', y = 'Latitude (°)', tag = 'a)') +
    guides(color = guide_legend(override.aes = list(size=2, alpha = 1)))

# b) temperature trends on land and freshwater
p2 <- ggplot(temptrends[REALM == 'Terrestrial & Freshwater'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans, 
                       breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2)) +
    scale_fill_manual(values=c('#E69F00', '#56B4E9'), name = 'Dataset') +
    labs(tag = 'b)', x = 'Temperature change [°C/yr]', title = 'Terrestrial & Freshwater') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = 'none',
          axis.text=element_text(size=7),
          plot.tag=element_text(face='bold'),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# c) temperature trends at sea
p3 <- ggplot(temptrends[REALM == 'Marine'], aes(x = tempchange, fill = type)) +
    geom_density(alpha = 0.4, color = NA) +
    scale_y_sqrt() +
    scale_x_continuous(limits = c(-2, 2.5), trans = signedsqrttrans, 
                       breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2)) +
    scale_fill_manual(values=c('#E69F00', '#56B4E9'), name = 'Dataset') +
    labs(tag = 'c)', x = 'Temperature change [°C/yr]', title = 'Marine') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = c(0.8, 0.95),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))
p3 <- addSmallLegend(p3, pointSize = 0.7, spaceLegend = 0.15, textSize = 7)

# d) example of Jtu trend calculation
p4 <- ggplot(bt[rarefyID=='339_1085477', .(dY = year2 - year1, Jtu.sc)], aes(dY, Jtu.sc)) +
    geom_point(alpha = 0.2, size = 0.5, shape = 16) +
    geom_smooth(method = 'glm', method.args = list(family = beta_family(link='logit')), color = '#a6cee3') + # a beta regression
    labs(tag = 'd)', x = 'Temporal distance [years]', y = 'Turnover\n[proportion species]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

# e) conceptual figure
conceptual <- data.table(expand.grid(duration = 1:40, tempchange = c('Fast', 'Slow')))
conceptual[tempchange == 'Fast', tchange := 2]
conceptual[tempchange == 'Slow', tchange := 0]
conceptual[, Jtu := 0.1 + duration*0.005 + duration*tchange/320]

p5 <- ggplot(conceptual, aes(duration, Jtu, color = tempchange)) +
    geom_line(size=1) +
    labs(tag = 'e)', x = 'Temporal distance [years]', y = 'Turnover\n[proportion species]') +
    scale_color_manual(values=c('#0072B2', '#D55E00'), name = expression(T[change])) +
    scale_y_continuous(limits = c(0, 0.5)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=7),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8))

fig1 <- arrangeGrob(p1, p2, p3, p4, p5, ncol = 4, 
                   layout_matrix = rbind(c(1,1,1,1), c(2,2,3,3), c(4,4,5,5)),
                   heights=c(unit(0.5, "npc"), unit(0.25, "npc"), unit(0.25, "npc")))

ggsave('figures/fig1.png', fig1, width = 6, height = 6, units = 'in')


### Figure 2: main effects ---------
# files to read in
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
slopespredsdT <- readRDS(here('temp', 'slopes_modOBsdTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMMmodrawXAllJtu.R
slopespred <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMMmodrawXAllJtu.R
senspred <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMMmodrawXAllJtu.R

# slopes for all timeseries
trends <- trends[duration_group == 'All' & measure == 'Jtu',]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, REALM, tempchange)])
trends[, duration := year2 - year1]
trends[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
trends[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part B
trends_by_study <- trends[duration>2, .(disstrend = mean(disstrend, na.rm=TRUE), tempchange = mean(tempchange, na.rm=TRUE), duration = mean(duration)), by = .(STUDY_ID, REALM)] # average by studyID. Can't use 2-year trends since they assume dissimilarity at y0 is 0.
trends_by_study[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part A

# average slopes by realm
ave_by_realm <- trends_by_study[, .(disstrend = mean(disstrend), se = sd(disstrend)/sqrt(.N), duration = mean(duration), duration.se = sd(duration)/sqrt(.N)), by = REALM]
ave_by_realm[, offset := c(-1, 0, 1)] # amount to vertically dodge the lines in part a
write.csv(ave_by_realm, file='output/ave_by_realm.csv')

# min and max tempchange by realm, for plotting limits
tempchange_by_realm <- trends[, .(max = quantile(tempchange, na.rm=TRUE, probs = 0.999), min = quantile(tempchange, na.rm=TRUE, probs = 0.001)), by = REALM]

# predicted slopes from the Tchange model
slopespredsdT <- merge(slopespredsdT, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespredsdT <- slopespredsdT[tempchange > min & tempchange < max & !duplicated(cbind(tempchange, REALM)), ] # trim to min & max by realm
slopespredsdT[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]

# predicted turnover and sensitivity of turnover rate to temperature change from the Tchange x Tave model
senspred <- senspred[tempave %in% c(0,25), ]
senspred[, tsign := factor(tsign, levels = c('-1', '1'), labels = c('cooling', 'warming'))]

# fastest turnover at highest observed rate of temperature change
slopespredsdT[, .SD[which.max(slope), .(tempchange, slope, slope.se)], by = REALM] # just looking at tempchange
slopespred[, .SD[which.max(slope), .(tempave, tempchange, slope, slope.se)], by = REALM] # also considering tempave

# predicted turnover at moderate rates of temperature change
slopespredsdT[, .SD[which.min(abs(tempchange - 0.5)), .(tempchange, slope, slope.se)], by = REALM] # just looking at tempchange
slopespred[, .SD[which.max(slope), .(tempave, tempchange, slope, slope.se)], by = REALM] # also considering tempave


# a) across realms
ht <- 6.3
p1 <- ggplot(trends_by_study, aes(x=disstrend, group = REALM, fill = REALM)) +
    geom_density(color = NA, alpha = 0.4) +
    scale_y_sqrt(breaks = c(0.1, 1, 2, 6)) +
    scale_x_continuous(trans = signedsqrttrans, breaks = c(-0.2, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.2, 0.4)) +
    geom_segment(data = ave_by_realm, aes(x=disstrend - 1.96*se, xend = disstrend + 1.96*se, y= ht+offset, yend = ht+offset, color = REALM), alpha = 1) +
    geom_segment(data = ave_by_realm, aes(x = disstrend, y = 0, xend = disstrend, yend = ht+offset, color = REALM), linewidth=0.5, linetype = 'dashed') +
    labs(tag = 'a)', x = expression('Turnover rate ['~Delta~'Turnover/year]'), y = 'Density', title = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_fill_brewer(palette = 'Dark2') +
    scale_color_brewer(palette = 'Dark2')
p1 <- addSmallLegend(p1, pointSize = 0.5, spaceLegend = 0.1, textSize = 6)

# b) plot of turnover rate vs. Tchange
p2 <- ggplot() +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(abs(tempchange), disstrend, size = duration, color = as.factor(tsign)), 
             color='#AAAAAA', alpha = 0.1, stroke = 0) +
    #geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_line(data = slopespredsdT, mapping=aes(abs(tempchange), slope, color = tsign, group = tsign), linewidth=0.5) +
    geom_ribbon(data = slopespredsdT, alpha = 0.2, color = NA, 
                aes(abs(tempchange), slope,
                    ymin=slope - slope.se, 
                    ymax=slope + slope.se,
                    fill = tsign,
                    group = tsign)) +
    scale_color_manual(values=c('#0072B2', '#D55E00')) +
    scale_fill_manual(values=c('#0072B2', '#D55E00')) +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'b)', x = 'Temperature change rate [|°C/year|]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Direction', 
         color = 'Direction',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(-1,-0.5, -0.1, -0.05, -0.01, 0, 0.01, 0.05, 0.1, 0.5, 1)) +
    scale_x_continuous(trans = signedsqrttrans,
                       breaks = c(-1, -0.5, -0.1, 0, 0.1, 0.5, 1)) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend
    
p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

# c) in three parts
p3 <- ggplot(senspred[REALM=='Terrestrial'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
    geom_point()+
    geom_line(linetype='dashed')+
    geom_errorbar()+
    facet_grid(col = vars(REALM))+
    labs(tag='c)',
         x = '', 
         y = '') +
    scale_color_manual(values=c('#0072B2', '#D55E00')) +
    scale_fill_manual(values=c('#0072B2', '#D55E00')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          legend.position='none', # no legend
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    coord_cartesian(clip = 'off') + # solution for multi-line y-axis from https://stackoverflow.com/questions/13223846/ggplot2-two-line-label-with-expression
    annotation_custom(textGrob(expression("Sensitivity of turnover rate"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -40, xmax = -40, ymin = 0.01, ymax = 0.01) +
    annotation_custom(textGrob(expression("to temperature change"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -35, xmax = -35, ymin = 0.01, ymax = 0.01) +
    annotation_custom(textGrob(expression('[('~Delta~'Turnover rate)/'~Delta~'°C/year)]'), rot = 90, gp = gpar(fontsize=6.5)), xmin = -30, xmax = -30, ymin = 0.01, ymax = 0.01) +
    scale_x_continuous(breaks=c(0,25), labels=c(0,25), limits=c(-10,35))


p4 <- ggplot(senspred[REALM=='Freshwater'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
    geom_point()+
    geom_line(linetype='dashed')+
    geom_errorbar()+
    facet_grid(col = vars(REALM))+
    labs(tag = '',
         x = 'Ave. temperature [°C]', 
         y = '') +
    scale_color_manual(values=c('#0072B2', '#D55E00')) +
    scale_fill_manual(values=c('#0072B2', '#D55E00')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position='none', # no legend
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    coord_cartesian(clip = 'off') + 
    scale_x_continuous(breaks=c(0,25), labels=c(0,25), limits=c(-10,35))

p5 <- ggplot(senspred[REALM=='Marine'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
    geom_point()+
    geom_line(linetype='dashed')+
    geom_errorbar()+
    facet_grid(col = vars(REALM))+
    labs(tag = '', x = '', 
         y = '',
         fill = 'Direction',
         color = 'Direction') +
    scale_color_manual(values=c('#0072B2', '#D55E00')) +
    scale_fill_manual(values=c('#0072B2', '#D55E00')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_x_continuous(breaks=c(0,25), labels=c(0,25), limits=c(-10,35))

p5 <- addSmallLegend(p5, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

fig2 <- arrangeGrob(p1, p2, p3, p4, p5, nrow = 3, ncol = 3, layout_matrix = matrix(c(1,1,1,2,2,2,3,4,5), nrow=3, byrow=TRUE), 
                    heights = c(1,2,1), widths = c(3,3,4))

ggsave('figures/fig2.png', fig2, width = 6, height = 6, units = 'in')


# w/out model predictions
p2noT <- ggplot() +
    geom_point(data = trends[!is.na(tempchange)], mapping = aes(abs(tempchange), disstrend, size = duration, color = as.factor(tsign)), 
               color='#AAAAAA', alpha = 0.1, stroke = 0) +
    #geom_hline(yintercept = 0, linetype = 'dotted') +
    #geom_line(data = slopespredsdT, mapping=aes(abs(tempchange), slope, color = tsign, group = tsign), size=0.5) +
    #geom_ribbon(data = slopespredsdT, alpha = 0.2, color = NA, 
    #            aes(abs(tempchange), slope,
    #                ymin=slope - slope.se, 
    #                ymax=slope + slope.se,
    #                fill = tsign,
    #                group = tsign)) +
    scale_color_brewer(palette = 'Dark2') +
    scale_fill_brewer(palette = 'Dark2') +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'b)', x = 'Temperature change rate [|°C/year|]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         fill = 'Direction', 
         color = 'Direction',
         size = 'Duration [years]') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=8),
          plot.title=element_text(size=8)) +
    scale_y_continuous(trans = signedsqrttrans, 
                       breaks = c(-1,-0.3, -0.1, -0.03, -0.01, 0, 0.01, 0.03, 0.1, 0.3, 1)) +
    scale_x_continuous(trans = signedsqrttrans,
                       breaks = c(-1, -0.5, -0.1, 0, 0.1, 0.5, 1)) +
    scale_size(trans='log', range = c(0.8,3), breaks = c(2, 5, 20, 50)) +
    guides(size = guide_legend(override.aes = list(alpha=1))) # set alpha to 1 for points in the legend
p2noT <- addSmallLegend(p2noT, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

fig2noT <- arrangeGrob(p1, p2noT, p3, p4, p5, nrow = 3, ncol = 3, layout_matrix = matrix(c(1,1,1,2,2,2,3,4,5), nrow=3, byrow=TRUE), 
                    heights = c(1,2,1), widths = c(3,3,4))
ggsave('figures/fig2_nopredsT.png', fig2noT, width = 6, height = 6, units = 'in')



### Figure 3: environmental interactions ---------
slopes2 <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu.rds')) # turnover rates from code/pred_GLMMmodrawCovariate.R
sensitivity2 <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu.rds')) # sensitivities from code/pred_GLMMmodrawCovariate.R
sensitivity2[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting

# max turnover rate by realm and covariate
slopes2[tempave==10 & tempchange == 2 & human_bowler %in% c(0,10), .(slope_microclim, slope_microclim.se, slope_human, slope_human.se), 
        by = .(REALM, microclim, human_bowler, tempchange)][order(REALM, microclim, tempchange)]

# plots
ylims.microclimate <- c(-0.05, 0.18)
ylims.human <- c(0, 0.18)

p1 <- ggplot(sensitivity2[REALM %in% c('Marine', 'Terrestrial')], 
             aes(microclim, sensitivity_microclim, 
                 ymin=sensitivity_microclim-1.96*sensitivity_microclim.se,  
                 ymax=sensitivity_microclim+1.96*sensitivity_microclim.se,
                 color = REALM,
                 group = REALM,
                 fill = REALM)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    labs(tag = 'a)', 
         x = 'Microclimate availability',
         y = '') +
    coord_cartesian(clip = 'off') + # solution for multi-line y-axis from https://stackoverflow.com/questions/13223846/ggplot2-two-line-label-with-expression
    annotation_custom(textGrob(expression("Sensitivity of turnover rate"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -2.65, xmax = -2.65, ymin = 0.005, ymax = 0.005) + # note x-axis is in log10 units
    annotation_custom(textGrob(expression("to temperature change"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -2.45, xmax = -2.45, ymin = 0.005, ymax = 0.005) +
    annotation_custom(textGrob(expression('[('~Delta~'Turnover rate)/'~Delta~'°C/year)]'), rot = 90, gp = gpar(fontsize=6.5)), xmin = -2.25, xmax = -2.25, ymin = 0.005, ymax = 0.005) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          legend.position='none', # no legend
          axis.text=element_text(size=8),
          axis.title=element_text(size=7),
          axis.title.y=element_text(vjust = 6),
          plot.title=element_text(size=8),
          plot.margin = unit(c(0.05, 0.07, 0.05, 0.05), 'in')) +
    #scale_x_log10(limits = c(0.03, 1)) +
    lims(y = ylims.microclimate) +
    scale_fill_manual(values = c('#1b9e77', '#7570b3')) + # ColorBrewer2 Dark2 green and blue
    scale_color_manual(values = c('#1b9e77', '#7570b3'))

p2 <- ggplot(sensitivity2[REALM %in% c('Marine', 'Terrestrial')], 
             aes(human_bowler, sensitivity_human, 
                 ymin=sensitivity_human-1.96*sensitivity_human.se,  
                 ymax=sensitivity_human+1.96*sensitivity_human.se,
                 color = REALM,
                 group = REALM,
                 fill = REALM)) +
    geom_ribbon(alpha = 0.25, color = NA, show.legend = FALSE) +
    geom_line() +
    labs(tag = 'b)', 
         x = 'Human impact', 
         y = '') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          plot.tag=element_text(face='bold'),
          axis.text=element_text(size=8),
          axis.title=element_text(size=7),
          axis.title.y=element_blank(),
          plot.title=element_text(size=8),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.01), 'in')) +
    scale_x_log10()  +
    lims(y = ylims.human) +
    scale_fill_manual(values = c('#66c2a5', '#8da0cb')) +
    scale_color_manual(values = c('#66c2a5', '#8da0cb'))
p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

fig3 <- arrangeGrob(p1, p2, ncol = 2, widths = c(3,4))

ggsave('figures/fig3.png', fig3, width = 4, height = 2, units = 'in')




### Table S1: random effects for main model --------------
modTchangeTave <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_?
if(!exists('sum_modTchangeTave')) sum_modTchangeTave <- summary(modTchangeTave)
sum_modTchangeTave$varcor
capture.output(print(sum_modTchangeTave$varcor), file = 'figures/tableS1.txt')



### Table S2: fixed effects for Tchange x Tave model --------------
modTchangeTave <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_?
if(!exists('sum_modTchangeTave')) sum_modTchangeTave <- summary(modTchangeTave)
out <- as.data.frame(sum_modTchangeTave$coefficients$cond)

# get term names
out$term <- gsub('Terrestrial|Marine|Freshwater|1|-1', '', rownames(out))
out$term <- gsub('duration', 'Years', out$term)
out$term <- gsub('REALM', 'Realm', out$term)
out$term <- gsub('tsign', 'sign', out$term)
out$term <- gsub('tempchange_abs.sc', 'Tchange', out$term)
out$term <- gsub('tempave.sc', 'Tave', out$term)
out$term[out$term == 'Years:Realm'] <- 'Realm ✕ Years'
out$term[out$term == 'Years:Realm:sign:Tchange'] <- '|Tchange| ✕ Years'
out$term[out$term == 'Years:Realm:sign:Tave'] <- 'Tave ✕ Years'
out$term[out$term == 'Years:Realm:sign:Tchange:Tave'] <- '|Tchange| ✕ Tave ✕ Years'

# get realm
out$realm <- rownames(out)
out$realm[grepl('Freshwater', out$realm)] <- 'Freshwater'
out$term[out$term == 'Years:Jtu.init'] <- 'Years:Dinit'
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

# reorder rows
out <- out[order(out$realm, out$tsign), ]

out

# write out
write.csv(format(out, digits=2), file = 'figures/tableS2.csv', row.names=FALSE)


### Table S3: Horn AICs --------------
modInitHorn <- readRDS(here('temp', 'modOBInitAllHorn.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
modRealmHorn <- readRDS('temp/modOBRealmInitAllHorn.rds') # Realm. Fit by code/turnover_GLMM_fit.R
modTaxamod2Horn <- readRDS('temp/modOBTaxamod2InitAllHorn.rds') # Taxon. Fit by code/turnover_GLMM_fit.R
modTchangeHorn <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllHorn.rds')) # Tchange model. Fit by code/turnover_GLMM_fit.R
modTchangeTaveHorn <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllHorn.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R

# compare Tchange amd Tchange x Tave models against null
aics <- AIC(modInitHorn, modRealmHorn, modTaxamod2Horn, # simple models w/out Tchange
            modTchangeHorn, # Tchange
            modTchangeTaveHorn) # Tchange x Tave
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInitHorn']
aics

write.csv(aics, here('figures', 'tableS3.csv'))


### Table S4: covariate AICs --------------
# load models for AICs
modInit <- readRDS(here('temp', 'modOBInitAllJtu.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
modTchange <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R
modmicroclim <- readRDS('temp/modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu.rds') # Microclimates. Fit by turnover_GLMM_fit.R.
modhuman <- readRDS('temp/modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu.rds') # Human impact. Fit by turnover_GLMM_fit.R
 
# compare covariate models against null
aics <- AIC(modInit, 
            modTchange, 
            modmicroclim,
            modhuman) 
aics$dAICtemp <- aics$AIC - aics$AIC[rownames(aics)=='modTchange']
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInit']
aics

write.csv(aics, here('figures', 'tableS4.csv'))





### Table S5: AICs for initgainloss models ------------------
# with Jtu.init:gainlossprop for Table 1. All fit by code/turnover_GLMMgainloss_fit.R
if(!exists('modInitGainLossAllJtu')) modInitGainLossAllJtu <- readRDS(here('temp', 'modInitGainLossAllJtu.rds')) # Null
if(!exists('modRealmInitGainLossAllJtu')) modRealmInitGainLossAllJtu <- readRDS('temp/modRealmInitGainLossAllJtu.rds') # Realm. 
if(!exists('modTaxamod2InitGainLossAllJtu')) modTaxamod2InitGainLossAllJtu <- readRDS('temp/modTaxamod2InitGainLossAllJtu.rds') # Taxon. 
if(!exists('modsdTRealmtsignInitGainLossAllJtu')) modsdTRealmtsignInitGainLossAllJtu <- readRDS(here('temp', 'modsdTRealmtsignInitGainLossAllJtu.rds')) # Tchange model
if(!exists('modabsLatsdTabsLatRealmtsignInitGainLossAllJtu')) modabsLatsdTabsLatRealmtsignInitGainLossAllJtu <- readRDS(here('temp', 'modabsLatsdTabsLatRealmtsignInitGainLossAllJtu.rds')) # Tchange x Latitude model
if(!exists('modrawTsdTTRealmtsignInitGainLossAllJtu')) modrawTsdTTRealmtsignInitGainLossAllJtu <- readRDS(here('temp','modrawTsdTTRealmtsignInitGainLossAllJtu.rds')) # Tchange x Tave model

aicsIGL <- AIC(modInitGainLossAllJtu, 
            modRealmInitGainLossAllJtu, 
            modTaxamod2InitGainLossAllJtu, # simple models w/out Tchange
            modsdTRealmtsignInitGainLossAllJtu, # Tchange model
            modabsLatsdTabsLatRealmtsignInitGainLossAllJtu,
            modrawTsdTTRealmtsignInitGainLossAllJtu) # Tchange x Tave model
aicsIGL$dAIC <- aicsIGL$AIC - min(aicsIGL$AIC)
aicsIGL$dAICnull <- aicsIGL$AIC - aicsIGL$AIC[rownames(aicsIGL)=='modInitGainLossAllJtu']
aicsIGL

write.csv(aicsIGL, here('figures', 'tableS5.csv'))



### Figure S1: time-series info----------
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
btts <- bt[, .(year1 = min(year1), year2 = max(year2), nsamp = length(unique(c(year1, year2)))), by = .(rarefyID, STUDY_ID)] # summarize by time-series

labpos <- -0.2 # horizontal position for the subfigure label

png(file = 'figures/figS1.png', width = 183, height = 183, units = 'mm', res = 300, pointsize=7)
par(mfrow=c(2,2), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(2.5, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a: start dates
btts[, hist(year1, main = '', xlab = '', col = 'grey')]
mtext('Start year', side = 1, line = 1.5, cex=0.8)
mtext('a)', side = 3, line = -0.5, adj = labpos, font = 2)

# part b: end dates
btts[, hist(year2, main = '', xlab = '', col = 'grey')]
mtext('End year', side = 1, line = 1.5, cex=0.8)
mtext('b)', side = 3, line = -0.5, adj = labpos, font = 2)

# part c: durations
btts[, hist(year2-year1+1, main = '', xlab = '', col = 'grey')]
mtext('Number of years', side = 1, line = 1.5, cex=0.8)
mtext('c)', side = 3, line = -0.5, adj = labpos, font = 2)

# part d: number of samples
btts[, hist(nsamp, main = '', xlab = '', col = 'grey')]
mtext('Number of samples', side = 1, line = 1.5, cex=0.8)
mtext('d)', side = 3, line = -0.5, adj = labpos, font = 2)

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
scalingall <- fread('output/turnover_w_covariates_scaling.csv') # covariate scaling data. From assemble_turnover_covariates.Rmd
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
png(file = 'figures/figS2.png', width = 183, height = 141, units = 'mm', res = 300, pointsize=7)
par(mfrow=c(2,3), mai = c(0.7, 0.7, 0.1, 0.1), las = 1, mgp = c(1.9, 0.5, 0), tcl = -0.2, cex.axis = 0.8)

# part a
bt[rarefyID == '339_1085477', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Turnover [proportion of species]', col = '#00000044', bty = 'l', ylim = c(0,1))]
bt[rarefyID == '339_1085477', abline(lm(Jtu~dY), col = '#a6cee3', lwd = 3)]
mod5 <- bt[rarefyID == '339_1085477' & dY <=5, lm(Jtu~dY)] # calc trendline
preds <- data.table(dY = 1:20)
preds$Jtu5 <- predict(mod5, preds)
preds[dY <=10, lines(dY, Jtu5, col = '#b2df8a', lwd = 3)]
mtext('a)', side = 3, line = -0.5, adj = -0.28, font = 2)

# part a inset
oldpar <- par(no.readonly=TRUE)
par(fig = c(0.05,0.3, 0.8, 1), new = T, mgp = c(1, 0.2, 0), cex.lab = 0.8, cex.axis = 0.8, tcl = -0.1)
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
mtext('b)', side = 3, line = -0.5, adj = -0.28, font = 2)


# part c: tempchange by duration
tempchange[, plot(duration, tempchange, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Temperature change [°C/yr]', bty = 'l')]
abline(h = 0, lty = 2)
predsloesstemp[, lines(duration, tempchange, col = 'red')]
mtext('c)', side = 3, line = -0.5, adj = -0.28, font = 2)


# part d: Gaussian white noise slope by duration
trends[, plot(duration, gauss.slope, cex=0.1, col = '#0000000F', xlab = 'Duration', ylab = 'Slope of Gaussian white noise', bty = 'l')]
abline(h = 0, lty = 2)
predsloess[, lines(duration, gauss.slope, col = 'red')]
mtext('d)', side = 3, line = -0.5, adj = -0.28, font = 2)


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
legend('topleft', legend = c('Pearson correlation', 'Meta-analysis', 'One-stage Gaussian ME', 'One-stage beta ME'), col = cols, pch = 1, cex=0.8)
mtext('e)', side = 3, line = -0.5, adj = -0.28, font = 2)


# part f: an example negative slope
bt[rarefyID == '213_435199', plot(dY, Jtu, xlab = 'Temporal difference (years)', ylab = 'Turnover [proportion of species]', col = '#00000044', bty = 'l', 
                                  ylim = c(0,1))]
mod <- glmmTMB(I(transform01(Jtu))~dY, data = bt[rarefyID == '213_435199',], family = beta_family(link = 'logit')) # calc trendline
preds <- data.table(dY = 1:35)
preds[, c('Jtu', 'se') := predict(mod, preds, se.fit=TRUE, type = 'response')]
preds[, polygon(x = c(dY, rev(dY)), y= c(Jtu+se, rev(Jtu-se)), col = '#88888855', border = NA)]
preds[, lines(dY, Jtu, col = '#a6cee3', lwd = 3)]
mtext('f)', side = 3, line = -0.5, adj = -0.28, font = 2)

dev.off()


### Figure S4: turnover by taxon ----------
# slopes for all timeseries
bt <- fread('output/turnover_w_covariates.csv.gz') # covariate data
trends <- fread('output/slope.csv.gz') # from calc_turnover.R
trends <- trends[duration_group == 'All' & measure == 'Jtu',] # trim to those we use
trends[, duration := year2 - year1]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, taxa_mod2)])
trends[taxa_mod2 == 'All', taxa_mod2 := 'Multiple taxa'] # rename a level so more intuitive
trends[taxa_mod2 == 'Benthos', taxa_mod2 := 'Benthic species'] # rename a level so more intuitive
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
          axis.text=element_text(size=7),
          axis.title=element_text(size=7),
          plot.title=element_text(size=7)) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal)
p1 <- addSmallLegend(p1, pointSize = 0.8, spaceLegend = 0.2, textSize = 8)
ggsave('figures/figS4.png', p1, width = 183, height = 122, units = 'mm')




### Figure S5: T_change x T_ave interaction ---------

# read in slopes
slopesTsdTTRealmtsigninit <- readRDS(here('temp', 'slopes_rawTsdTTRealmtsigninit.rds')) # made by pred_modrawXAllJtu.R

# plot
p1 <- ggplot(slopesTsdTTRealmtsigninit, aes(tempchange, tempave, z = slope)) +
    geom_raster(aes(fill = slope)) +
    labs(x = 'Temperature change (°C/yr)', y = 'Average Temperature (°C)') +
    scale_fill_gradient2(high= "#B2182B", mid = "white", low= "#2166AC", midpoint = 0, name = 'Turnover rate') +
    facet_grid(cols = vars(REALM)) +
    theme(axis.text = element_text(size = 5), 
          axis.title = element_text(size = 7),
          strip.text = element_text(size=7),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "top",
          legend.margin = margin(c(0, 0.1, -5, 0)),
          legend.justification = c(0.92, 0.9),
          legend.text=element_text(size= 5),
          legend.title=element_text(size= 7),
          legend.title.align = 1)

ggsave('figures/figS5.png', p1, width = 183, height = 92, units = 'mm')
