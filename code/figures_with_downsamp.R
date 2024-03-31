# Making tables and figures with the bootstrap model fits

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
source(here('code', 'util.R'))





### AIC, BIC, and Likelihood ratio tests among models ----------------
# read in each bootstrapped model fit
n1 <- 10L # number of bootstraps
modTchange <- vector(mode = 'list', length = n1)
modTchangeYear <- vector(mode = 'list', length = n1)
modTchangeTave <- vector(mode = 'list', length = n1)
n2 <- 3L
modTchangeTaveYear <- vector(mode = 'list', length = n2)
modmicroclim <- vector(mode = 'list', length = n2)
modhuman <- vector(mode = 'list', length = n2)
for(i in 1:n1){
    cat(i)
    modTchange[[i]] <- readRDS(here('temp', paste0('modOBMERtsRealmtsignTchangeinitAllJtu_boot', i, '.rds'))) # Tchange x Realm model. Fit by code/turnover_GLMM_fit.R
    modTchangeYear[[i]] <- readRDS(here('temp', paste0('modOBsdTMERtsRealmtsigninitAllJtu_boot', i, '.rds'))) # Tchange x Realm x Year model. Fit by code/turnover_GLMM_fit.R
    modTchangeTave[[i]] <- readRDS(here('temp', paste0('modOBMERtsRealmtsignTchangeTaveinitAllJtu_boot', i, '.rds'))) # Tchange x Tave x Realm model. Fit by code/turnover_GLMM_fit.R
}

for(i in 1:n2){
    cat(i)
    modTchangeTaveYear[[i]] <- readRDS(here('temp', paste0('modOBrawTsdTTMERtsRealmtsigninitAllJtu_boot', i, '.rds'))) # Tchange x Tave x Year x Realm model. Fit by code/turnover_GLMM_fit.R
    modmicroclim[[i]] <- readRDS(here('temp', paste0('modOBrawTsdTTMERtsRealmtsignmicroclimInitAllJtu_boot', i, '.rds'))) # Microclimates. Fit by turnover_GLMM_fit.R.
    modhuman[[i]] <- readRDS(here('temp', paste0('modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_boot', i, '.rds'))) # Human impact. Fit by turnover_GLMM_fit.R
}

# modTchangeTaveYearTerr <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllJtu_terr.rds')) # Tchange x Tave x Year terrestrial model. Fit by code/turnover_GLMM_fit.R
# modhumanTerr <- readRDS('temp/modOBrawTsdTTMERtsRealmtsignhumanInitAllJtu_terr.rds') # Human impact terrestrial model. Fit by turnover_GLMM_fit.R

anovamodTchangemodTchangeYear <- data.frame(boot = 1:n1, dAIC = numeric(n1), dBIC = numeric(n1), Chisq = numeric(n1), ChiDf = numeric(n1), p = numeric(n1)) # for Tchange vs. TchangeYear
for(i in 1:n1){
    temp <- anova(modTchange[[i]], modTchangeYear[[i]]) # Tchange x Realm vs. Tchange x Year x Realm model
    anovamodTchangemodTchangeYear[i, ] <- c(i, diff(temp$AIC), diff(temp$BIC), temp$Chisq[2], temp$`Chi Df`[2], temp$`Pr(>Chisq)`[2])
}
anovamodTchangemodTchangeYear

anovamodTchangeTavemodTchangeTaveYear <- data.frame(boot = 1:n2, dAIC = numeric(n2), dBIC = numeric(n2), Chisq = numeric(n2), ChiDf = numeric(n2), p = numeric(n2)) # TchangeTave vs. TchangeTaveYear
for(i in 1:n2){
    temp <- anova(modTchangeTave[[i]], modTchangeTaveYear[[i]]) # Tchange x Tave x Realm vs. Tchange x Tave x Year x Realm model
    anovamodTchangeTavemodTchangeTaveYear[i, ] <- c(i, diff(temp$AIC), diff(temp$BIC), temp$Chisq[2], temp$`Chi Df`[2], temp$`Pr(>Chisq)`[2])
}
anovamodTchangeTavemodTchangeTaveYear

anovamodTchangeTaveYearmodmicroclim <- data.frame(boot = 1:n2, dAIC = numeric(n2), dBIC = numeric(n2), Chisq = numeric(n2), ChiDf = numeric(n2), p = numeric(n2)) # TchangeTaveYear vs. microclim
for(i in 1:n2){
    temp <- anova(modTchangeTaveYear[[i]], modmicroclim[[i]]) # Tchange x Tave x Realm x Year model vs. microclimate
    anovamodTchangeTaveYearmodmicroclim[i, ] <- c(i, diff(temp$AIC), diff(temp$BIC), temp$Chisq[2], temp$`Chi Df`[2], temp$`Pr(>Chisq)`[2])
}
anovamodTchangeTaveYearmodmicroclim

anovamodTchangeTaveYearmodhuman <- data.frame(boot = 1:n2, dAIC = numeric(n2), dBIC = numeric(n2), Chisq = numeric(n2), ChiDf = numeric(n2), p = numeric(n2)) # TchangeTaveYear vs. human
for(i in 1:n2){
    temp <- anova(modTchangeTaveYear[[i]], modhuman[[i]]) # Tchange x Tave x Realm x Year model vs. human
    anovamodTchangeTaveYearmodhuman[i, ] <- c(i, diff(temp$AIC), diff(temp$BIC), temp$Chisq[2], temp$`Chi Df`[2], temp$`Pr(>Chisq)`[2])
}
anovamodTchangeTaveYearmodhuman


# anova(modTchangeTaveYearTerr, modhumanTerr) # Tchange x Tave x Year model vs. human, terrestrial only



### Tchange effects ---------
# files to read in
n1 <- 10L # number of bootstraps
for(i in 1:n1){
    cat(i)
    if(i==1){
        slopespredsdT <- cbind(data.table(boot =i, type = 'downsamp'), 
                               readRDS(here('temp', paste0('slopes_modOBsdTMERtsRealmtsigninitAllJtu_boot', i, '.rds')))) # from pred_GLMM_boot.R  
    } 
    else{
        slopespredsdT <- rbind(slopespredsdT,
                               cbind(data.table(boot =i, type = 'downsamp'), 
                                     readRDS(here('temp', paste0('slopes_modOBsdTMERtsRealmtsigninitAllJtu_boot', i, '.rds'))))) # from pred_GLMM_boot.R  
    } 
}
slopespredsdT[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
slopespredsdT <- slopespredsdT[tempave ==15,] # no Tave effect, so trim out the various levels
slopespredsdTfull <- readRDS(here('temp', 'slopes_modOBsdTMERtsRealmtsigninitAllJtu.rds')) # from pred_GLMM.R
slopespredsdTfull[, tsign := factor(sign(tempchange), levels = c('-1', '1'), labels = c('cooling', 'warming'))]
slopespredsdTfull[, ':='(boot = NA_integer_, type = 'full')]
slopespredsdT <- rbind(slopespredsdT, slopespredsdTfull)

# b) plot of turnover rate vs. Tchange
ggplot(slopespredsdT, aes(x=tempchange, y=slope,
                          ymin=slope - slope.se,
                          ymax=slope + slope.se,
                          group = boot)) +
    geom_line(linewidth=0.5, aes(linetype = type, color = type)) +
    geom_ribbon(alpha = 0.2, aes(fill = type)) +
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
    




### Figure 3: environmental interactions ---------
slopes2 <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu.rds')) # turnover rates from code/pred_GLMMmodrawCovariate.R
sensitivity2 <- readRDS(here('temp', 'sensitivity_modOBrawTsdTTMERtsRealmtsignCovariateInitAllJtu.rds')) # sensitivities from code/pred_GLMMmodrawCovariate.R
sensitivity2[, REALM := factor(REALM, levels = c('Terrestrial', 'Freshwater', 'Marine'))] # re-order for nicer plotting

# max turnover rate by realm and covariate
slopes2[tempave==10 & abs(tempchange - 0.3) < 0.02 & (abs(human_bowler) < 0.1 | abs(human_bowler - 10) < 0.1), .(slope_microclim, slope_microclim.se, slope_human, slope_human.se), 
        by = .(REALM, microclim, human_bowler, tempchange)][order(REALM, microclim, tempchange)]

# ratio of homogenous vs. heterogeneous
slopes2[REALM == 'Terrestrial' & tempave==10 & abs(tempchange - 0.3) < 0.02 & abs(microclim  - 0.02) < 0.01, .(slope_microclim)] / 
    slopes2[REALM == 'Terrestrial' & tempave==10 & abs(tempchange - 0.3) < 0.02 & abs(microclim  - 1.14) < 0.01, .(slope_microclim)]

# ratio of human impacted vs. not
slopes2[REALM == 'Terrestrial' & tempave==10 & abs(tempchange - 0.3) < 0.02 & abs(human_bowler  - 10) < 0.01, .(slope_human)] / 
    slopes2[REALM == 'Terrestrial' & tempave==10 & abs(tempchange - 0.3) < 0.02 & abs(human_bowler  - 0.055) < 0.01, .(slope_human)]

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
    annotation_custom(textGrob(expression("Sensitivity of turnover rate"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -0.75, xmax = -0.75, ymin = 0.05, ymax = 0.05) + # note x-axis is in log10 units
    annotation_custom(textGrob(expression("to temperature change"), rot = 90, gp = gpar(fontsize=6.5)), xmin = -0.625, xmax = -0.625, ymin = 0.05, ymax = 0.05) +
    annotation_custom(textGrob(expression('[('~Delta~'Turnover rate)/'~Delta~'°C/year)]'), rot = 90, gp = gpar(fontsize=6.5)), xmin = -0.5, xmax = -0.5, ymin = 0.05, ymax = 0.05) +
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
modInitHorn <- readRDS(here('temp', 'modOBRInitAllHorn.rds')) # Null with only duration. Fit by code/turnover_GLMM_fit.R
modRealmHorn <- readRDS('temp/modOBRRealmInitAllHorn.rds') # Realm. Fit by code/turnover_GLMM_fit.R
modTaxamod2Horn <- readRDS('temp/modOBTTaxamod2InitAllHorn.rds') # Taxon. Fit by code/turnover_GLMM_fit.R
modTchangeHorn <- readRDS(here('temp', 'modOBMERtsRealmtsignTchangeinitAllHorn.rds')) # Tchange model. Fit by code/turnover_GLMM_fit.R
modTchangeYearHorn <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitAllHorn.rds')) # Tchange model. Fit by code/turnover_GLMM_fit.R
modTchangeTaveHorn <- readRDS(here('temp','modOBMERtsRealmtsignTchangeTaveinitAllHorn.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R
modTchangeTaveYearHorn <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitAllHorn.rds')) # Tchange x Tave model. Fit by code/turnover_GLMM_fit.R

# compare Tchange amd Tchange x Tave models against null
aics <- AIC(modInitHorn, modRealmHorn, modTaxamod2Horn, # simple models w/out Tchange
            modTchangeHorn, # Tchange
            modTchangeYearHorn, # Tchange x Year
            modTchangeTaveHorn, # Tchange x Tave
            modTchangeTaveYearHorn) # Tchange x Tave x Year
aics$dAIC <- aics$AIC - min(aics$AIC)
aics$dAICnull <- aics$AIC - aics$AIC[rownames(aics)=='modInitHorn']
aics

write.csv(aics, here('figures', 'tableS3.csv'))


### Table S4: covariate AICs --------------
# load models for AICs
modInit <- readRDS(here('temp', 'modOBRInitAllJtu.rds')) # Null with duration and realm. Fit by code/turnover_GLMM_fit.R
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
# with Jtu.init:gainlossprop for Table 1. All fit by code/turnover_GLMMordbeta_fit.R
modInitGL <- readRDS(here('temp', 'modOBRInitGLAllJtu.rds')) # Null. Switch to modOBRInitGLAllJtu
modRealmGL <- readRDS('temp/modOBRRealmInitGLAllJtu.rds') # Realm.
modTaxaGL <- readRDS('temp/modOBTTaxamod2InitGLAllJtu.rds') # Taxon. 
modTchangeGL <- readRDS(here('temp', 'modOBsdTMERtsRealmtsigninitGLAllJtu.rds')) # Tchange model
modTchangeTaveGL <- readRDS(here('temp','modOBrawTsdTTMERtsRealmtsigninitGLAllJtu.rds')) # Tchange x Tave model

aicsIGL <- AIC(modInitGL, 
            modRealmGL, 
            modTaxaGL, # simple models w/out Tchange
            modTchangeGL, # Tchange model
            modTchangeTaveGL) # Tchange x Tave model
aicsIGL$dAIC <- aicsIGL$AIC - min(aicsIGL$AIC)
aicsIGL$dAICnull <- aicsIGL$AIC - aicsIGL$AIC[rownames(aicsIGL)=='modInitGL']
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
modloess <- trends[, loess(disstrend~duration)] # loess fit (slow)
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
mod <- glmmTMB(Jtu~dY, data = bt[rarefyID == '213_435199',], family = ordbeta(link = 'logit')) # calc trendline
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
slopesTchangeTave <- readRDS(here('temp', 'slopes_modOBrawTsdTTMERtsRealmtsigninitAllJtu.rds')) # made by pred_modrawXAllJtu.R

# plot
p1 <- ggplot(slopesTchangeTave, aes(tempchange, tempave, z = slope)) +
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



