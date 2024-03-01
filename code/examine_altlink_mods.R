# Compare ordbeta and gaussian error models
# Uses output from turnover_GLMMaltlink_fit.R and pred_GLMMaltlink.R

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(here)
source(here('code', 'util.R'))


### AIC like Table 1 -------------------------
# load models
modInitlin <- readRDS(here('temp', 'modOBRInitAllJtu_lin.rds')) # Null with only duration and realm.
modRealmlin <- readRDS('temp/modOBRRealmInitAllJtu_lin.rds') # Realm:duration. 
modTaxamod2lin <- readRDS('temp/modOBTTaxamod2InitAllJtu_lin.rds') # Taxon:duration. 
modTchangelin <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllJtu_lin.rds')) # Tchange model. 
modTchangeTavelin <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds')) # Tchange x Tave model. 

aics <- AIC(modInitlin, modRealmlin, modTaxamod2lin, # simple models w/out Tchange
            modTchangelin, # Tchange
            modTchangeTavelin) # add Tchange x Tave
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInitlin']
aics

write.csv(aics, here('figures', 'table1_gaussian.csv'))



### Plot main effects like Fig. 2 --------------------------
slopespredsdT <- readRDS(here('temp', 'slopes_modOBsdTMERtsRealmtsigninitAllJtu_lin.rds'))
slopespred <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds'))
senspred <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds'))
bt <- fread('output/turnover_w_covariates.csv.gz') # from assemble_turnover_covariates.Rmd
trends <- fread('output/slope.csv.gz')

# process slopes from the Tchange model
trends <- trends[duration_group == 'All' & measure == 'Jtu',]
trends <- merge(trends, bt[!duplicated(rarefyID),. (rarefyID, STUDY_ID, REALM, tempchange)])
trends[, duration := year2 - year1]
trends[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
trends[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting in part B
trends_by_study <- trends[duration>2, .(disstrend = mean(disstrend, na.rm=TRUE), tempchange = mean(tempchange, na.rm=TRUE), duration = mean(duration)), by = .(STUDY_ID, REALM)] # average by studyID. Can't use 2-year trends since they assume dissimilarity at y0 is 0.
#tempchange_by_realm <- trends[, .(max = quantile(tempchange, na.rm=TRUE, probs = 0.999), min = quantile(tempchange, na.rm=TRUE, probs = 0.001)), by = REALM] # min and max tempchange by realm, for plotting limits
tempchange_by_realm <- trends[, .(max = max(tempchange, na.rm=TRUE), min = min(tempchange, na.rm=TRUE)), by = REALM] # min and max tempchange by realm, for plotting limits
slopespredsdT <- merge(slopespredsdT, tempchange_by_realm, all.x = TRUE, by = "REALM") # add min and max by realm
slopespredsdT <- slopespredsdT[tempchange > min & tempchange < max & !duplicated(cbind(tempchange, REALM)), ] # trim to min & max by realm
slopespredsdT[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]

# predicted turnover and sensitivity of turnover rate to temperature change from the Tchange x Tave model
senspred <- senspred[tempave %in% c(0,25), ]
senspred[, tsign := factor(tsign, levels = c('-1', '1'), labels = c('cooling', 'warming'))]


# plot of turnover rate vs. Tchange
p1 <- ggplot() +
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
    labs(x = 'Temperature change rate [|°C/year|]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
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

p1 <- addSmallLegend(p1, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

ggsave('figures/fig2_gaussian.png', p1, width = 6, height = 3, units = 'in')




### Environmental covariate plots ------------------
slopes2 <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_lin.rds')) 
sensitivity2 <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu_lin.rds')) 
sensitivity2[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting

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
#    lims(y = ylims.microclimate) +
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
#    lims(y = ylims.human) +
    scale_fill_manual(values = c('#66c2a5', '#8da0cb')) +
    scale_color_manual(values = c('#66c2a5', '#8da0cb'))
p2 <- addSmallLegend(p2, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

p1
p2

fig3 <- arrangeGrob(p1, p2, ncol = 2, widths = c(3,4))

ggsave('figures/fig3_gaussian.png', fig3, width = 4, height = 2, units = 'in')

