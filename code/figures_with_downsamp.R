# Making tables and figures with the downsampled model fits

### Functions -----------
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(glmmTMB) # for beta regression
library(RColorBrewer)
library(scales) # for defining signed sqrt axis transformation
library(here)
source(here('code', 'util.R'))






### Tchange effects ---------
# files to read in
slopes <- list.files(path = 'temp', pattern = glob2rx('slopes_modOBsdTMERtsRealmtsigninitAllJtu_boot*.rds'), full.names=TRUE) # from pred_GLMM_downsamp.R or fit_pred_turnover_GLMM_downsamp.R
n <- length(slopes)
n # number of downsample slopes made
n <- min(1000, n) # take first 1000
for(i in 1:n){
    if(i %% 20 == 0) cat(paste0(i,','))
    temp <- readRDS(slopes[i])
    if(i==1){ 
        slopespredsdT <- cbind(data.table(boot =i, type = 'downsamp'), readRDS(slopes[i])) 
    } 
    else{
        slopespredsdT <- rbind(slopespredsdT, cbind(data.table(boot =i, type = 'downsamp'), 
                                                    readRDS(slopes[i])))
    }
}
slopespredsdT[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
slopespredsdT <- slopespredsdT[tempave ==15,] # no Tave effect, so trim out the various levels
slopespredsdTfull <- readRDS(here('temp', 'slopes_modOBsdTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMM.R
slopespredsdTfull[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
slopespredsdTfull[, ':='(boot = NA_integer_, type = 'full')]
slopespredsdT <- rbind(slopespredsdT, slopespredsdTfull)

slopeave <- slopespredsdT[type=='downsamp', .(boot=1, slope = mean(slope), slope.u95 = quantile(slope, 0.975), slope.l95 = quantile(slope, 0.025)), by = .(tempchange, REALM)]

# b) plot of turnover rate vs. Tchange
ggplot(slopespredsdT[type=='downsamp',], aes(x=tempchange, y=slope,
                          group = boot)) +
    geom_ribbon(data = slopespredsdT[type=='full',], alpha = 0.2, fill = '#D55E00', aes(ymin=slope-slope.se, ymax=slope+slope.se)) +
    geom_line(data = slopespredsdT[type=='full',], color = '#D55E00', linewidth = 1) +
    geom_ribbon(data = slopeave, alpha = 0.25, fill = '#0072B2', aes(ymin=slope.l95, ymax=slope.u95)) +
    geom_line(alpha = 0.15, linewidth = 0.1, color = '#0072B2') +
    geom_line(data = slopeave, color = '#0072B2', linewidth = 1) +
    facet_grid(cols = vars(REALM), scales = 'free')  +
    labs(tag = 'b)', x = 'Temperature change rate [|°C/year|]', y = expression(atop('Turnover rate','['~Delta~'Turnover/year]')), 
         color = 'Type') +
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
    


### Tchange x Tave effects ---------

# downsampled files to read in
sensfiles <- list.files(path = 'temp', pattern = glob2rx('sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu_boot*.rds'), full.names=TRUE) # from fit_pred_turnover_GLMM_downsamp.R
n <- length(sensfiles)
n # number of downsample sensitivities made
n <- min(1000, n) # take first 1000
for(i in 1:n){
    if(i %% 20 == 0) cat(paste0(i,','))
    temp <- readRDS(sensfiles[i])
    if(i==1){ 
        temp <- readRDS(sensfiles[i])
        temp <- temp[tempave %in% c(0,25), ]
        senspred <- cbind(data.table(boot =i, type = 'downsamp'), temp) 
    } 
    else{
        temp <- readRDS(sensfiles[i])
        temp <- temp[tempave %in% c(0,25), ]
        senspred <- rbind(senspred, cbind(data.table(boot =i, type = 'downsamp'), temp))
    }
}
senspred[, tsign := factor(tsign, levels = c('-1', '1'), labels = c('cooling', 'warming'))]
senspred[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting

sensave <- senspred[type=='downsamp', .(boot=1, sensitivity = mean(sensitivity), sensitivity.u95 = quantile(sensitivity, 0.975), sensitivity.l95 = quantile(sensitivity, 0.025)), 
                    by = .(tempave, tsign, REALM)] # average across the downsamples

# read in fit to the full dataset
sensfull <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMM.R, Tchange x Tave x Realm x Year model
sensfull <- sensfull[tempave %in% c(0,25), ]
sensfull[, tsign := factor(tsign, levels = c('-1', '1'), labels = c('cooling', 'warming'))]



# plot of sensitivity vs. Tave
ggplot(senspred, aes(tempave, sensitivity, group = interaction(tsign, boot), color = tsign, fill = tsign)) +
    geom_line(alpha = 0.15, linewidth = 0.1, position = position_dodge(width=2)) + # the downsampled fits
    geom_errorbar(alpha = 0.15, linewidth = 0.1, 
                  aes(ymin = sensitivity-1.96*sensitivity.se, ymax = sensitivity+1.96*sensitivity.se),
                  position = position_dodge(width=2)) +
    geom_line(data = sensfull, alpha = 1, linewidth = 1, linetype = 'dashed', aes(group = tsign)) + # the fit to the full dataset
    geom_errorbar(data = sensfull, alpha = 1, linewidth = 1, width = 0,
                  aes(group = tsign, ymin = sensitivity-1.96*sensitivity.se, ymax = sensitivity+1.96*sensitivity.se)) +
    geom_line(data = sensave, alpha = 1, linewidth = 1, position = position_dodge(width=3)) + # the average across downsampled fits
    geom_errorbar(data = sensave, alpha = 1, linewidth = 1, width = 0,
                  aes(ymin = sensitivity.l95, ymax = sensitivity.u95),
                  position = position_dodge(width=3)) +
    facet_grid(col = vars(REALM)) +
    labs(tag = '',
         x = 'Ave. temp. [°C]', 
         y = 'Sensitivity') +
    scale_color_manual(values=c('#0072B2', '#D55E00')) +
    scale_fill_manual(values=c('#0072B2', '#D55E00')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position='none', # no legend
          axis.text=element_text(size=7),
          axis.title=element_text(size=7),
          plot.title=element_text(size=7)) +
    coord_cartesian(clip = 'off') + 
    scale_x_continuous(breaks=c(0,25), labels=c(0,25), limits=c(-10,35))

