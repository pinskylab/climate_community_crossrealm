# Compare ordbeta and gaussian error models

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(here)
source(here('code', 'util.R'))


### AIC -------------------------
# load models
modInit <- readRDS(here('temp', 'modOBRInitAllJtu.rds')) # Null with only duration and realm. Fit by code/turnover_GLMM_fit.R
modRealm <- readRDS('temp/modOBRRealmInitAllJtu.rds') # Realm:duration. Fit by code/turnover_GLMM_fit.R
modTaxamod2 <- readRDS('temp/modOBTTaxamod2InitAllJtu.rds') # Taxon:duration. Fit by code/turnover_GLMM_fit.R
modTchange <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllJtu.rds')) # Tchange model. Fit by code/turnover_GLMM_fit.R
modLat <- readRDS(here('temp', 'modOBabsLatsdTabsLatMERtsRealmtsignInitAllJtu.rds')) # Lat model. Fit by code/turnover_GLMM_fit.R
modTchangeTave <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R

modInitlin <- readRDS(here('temp', 'modOBRInitAllJtu_lin.rds')) # Null with only duration and realm.
modRealmlin <- readRDS('temp/modOBRRealmInitAllJtu_lin.rds') # Realm:duration. 
modTaxamod2lin <- readRDS('temp/modOBTTaxamod2InitAllJtu_lin.rds') # Taxon:duration. 
modTchangelin <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllJtu_lin.rds')) # Tchange model. 
modLatlin <- readRDS(here('temp', 'modOBabsLatsdTabsLatMERtsRealmtsignInitAllJtu_lin.rds')) # Lat model. 
modTchangeTavelin <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds')) # Tchange x Tave model. 


# compare beta and gaussian
aics <- AIC(modInit, modInitlin,
            modRealm, modRealmlin,
            modTaxamod2, modTaxamod2lin, # simple models w/out Tchange
            modTchange, modTchangelin, # Tchange
            modLat, #modLatlin, # latitude
            modTchangeTave, modTchangeTavelin) # add Tchange x Tave
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInit']
aics

# compare within gaussian
aics <- AIC(modInitlin, modRealmlin, modTaxamod2lin, # simple models w/out Tchange
            modTchangelin, # Tchange
            #modLatlin, # latitude
            modTchangeTavelin) # add Tchange x Tave
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInitlin']
aics



### Plot main effects --------------------------
slopespredsdT <- readRDS(here('temp', 'slopes_modOBsdTMERtsRealmtsigninitAllJtu_lin.rds'))
slopespred <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds'))
senspred <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu_lin.rds'))

# predicted slopes from the Tchange model
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

p1 <- addSmallLegend(p1, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

# sensitivity by Tave
p2 <- ggplot(senspred[REALM=='Terrestrial'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
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


p3 <- ggplot(senspred[REALM=='Freshwater'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
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

p4 <- ggplot(senspred[REALM=='Marine'], aes(tempave, sensitivity, ymin = sensitivity-sensitivity.se, ymax = sensitivity+sensitivity.se, group = tsign, color = tsign, fill = tsign)) +
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

p4 <- addSmallLegend(p4, pointSize = 0.8, spaceLegend = 0.1, textSize = 6)

