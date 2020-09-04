## Making publication-ready figures

#################
# Functions
#################
library(data.table) # for handling large datasets
library(ggplot2) # for some plotting
library(nlme) # for ME models
library(maps) # for map
library(gridExtra) # to combine ggplots together
library(grid) # to combine ggplots together
library(RColorBrewer)
library(scales) # for defining signed sqrt axis transformation

# produce ggplot-style colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}



####################
## Figure 1: map
####################
# load BioTime data
trends <- fread('output/turnover_w_covariates.csv.gz')

# load sampled temperature trends
temptrends <- fread('output/temperature_trends_sampled.csv.gz')

# trim BT to >= 3 yrs
trends <- trends[nyrBT >= 3, ]

# make trends lists
temptrends$type <- 'Global'
temptrends <- temptrends[, .(trend, type, REALM)]
trends$type <- 'BioTime'
temptrends <- rbind(temptrends, trends[!is.na(temptrend), .(trend = temptrend, type, REALM)])
temptrends[REALM %in% c('Terrestrial', 'Freshwater'), REALM := 'Terrestrial & Freshwater']

# make plot
world <- map_data('world')
p1 <- ggplot(world, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill = 'lightgray', color = 'white') +
    geom_point(data = trends, aes(rarefyID_x, rarefyID_y, group = REALM, color = REALM), size = 0.5, alpha = 0.4)  +
    scale_color_brewer(palette="Set1", name = 'Realm') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16)) +
    labs(x = 'Longitude (°)', y = 'Latitude (°)', tag = 'A)')

p2 <- ggplot(temptrends[REALM == 'Terrestrial & Freshwater'], aes(x = trend, fill = type)) +
    geom_density(alpha = 0.25) +
    scale_y_sqrt(limits = c(0,16)) +
    scale_x_continuous(limits = c(-2, 2.5)) +
    labs(tag = 'B)', x = 'Temperature trend (°C/yr)', title = 'Terrestrial & Freshwater') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = 'none',
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

p3 <- ggplot(temptrends[REALM == 'Marine'], aes(x = trend, fill = type)) +
    geom_density(alpha = 0.25) +
    scale_y_sqrt(limits = c(0,16)) +
    scale_x_continuous(limits = c(-2, 2.5)) +
    labs(tag = 'C)', x = 'Temperature trend (°C/yr)', title = 'Marine') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = c(0.8, 0.95),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

p4 <- ggplot(trends, aes(x = Jbetatrendrem0)) +
    geom_density() +
    scale_y_sqrt() +
    labs(tag = 'D)', x = 'Slope', title = 'Jaccard') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

p5 <- ggplot(trends, aes(x = Jtutrendrem0)) +
    geom_density() +
    scale_y_sqrt() +
    labs(tag = 'E)', x = 'Slope', title = 'Jaccard turnover') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

p6 <- ggplot(trends, aes(x = Horntrendrem0)) +
    geom_density() +
    scale_y_sqrt() +
    labs(tag = 'F)', x = 'Slope', title = 'Horn-Morisita') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

fig1 <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 6, 
                   layout_matrix = rbind(c(1,1,1,1,1,1), c(2,2,2,3,3,3), c(4,4,5,5,6,6)),
                   heights=c(unit(0.5, "npc"), unit(0.25, "npc"), unit(0.25, "npc")))

ggsave('figures/fig1.png', fig1, width = 6, height = 6, units = 'in')


#########################
## Figure 2: main effects
#########################
trendsum <- fread('output/trendsummary.csv')
newdat <- read.csv('output/maineffects.csv')
newdat <- newdat[newdat$var == 'temptrend_abs',]
aicsfromfull <- read.csv('output/aics_from_full.csv')
aicsfromfull$row <- nrow(aicsfromfull):1


p1 <- ggplot(trendsum[text %in% c('Changing', 'Stable')], aes(text, ave, group = type, color = type)) +
    geom_point(position = position_dodge(width = 0.25), size = 0.5) +
    geom_errorbar(aes(ymin=ave-se, ymax=ave+se), width=.1, position = position_dodge(width = 0.25)) +
    labs(x = '', y = 'Slope', title = 'Raw data', tag = 'A') +
    geom_abline(intercept = 0, slope = 0, linetype = 'dashed') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.position = 'none',
          axis.text=element_text(size=8),
          axis.title=element_text(size=10)) +
    coord_cartesian(ylim = c(0,0.04))



cols = gg_color_hue(3)
colors <- c("Jaccard turnover" = cols[3], "Jaccard total" = cols[2], "Horn-Morisita" = cols[1])

p2 <- ggplot(newdat, aes(x = temptrend_abs)) +
    geom_line(aes(y = predsJtu, color = 'Jaccard turnover')) +
    geom_ribbon(aes(ymin = predsJtu - 1.96*SE_Jtu, ymax = predsJtu + 1.96*SE_Jtu), alpha = 0.2, fill = "black") +
    geom_line(aes(y = predsJbeta, color = 'Jaccard total')) +
    geom_ribbon(aes(ymin = predsJbeta - 1.96*SE_Jbeta, ymax = predsJbeta + 1.96*SE_Jbeta), alpha = 0.2, fill = "red") +
    geom_line(aes(y = predsHorn, color = 'Horn-Morisita')) +
    geom_ribbon(aes(ymin = predsHorn - 1.96*SE_Horn, ymax = predsHorn + 1.96*SE_Horn), alpha = 0.2, fill = "blue") +
    theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
    labs(tag = 'B', x = 'Temperature trend (°C/yr)', y = 'Slope', title = 'Model effect') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.title=element_blank(),
          legend.position='none',
          axis.text=element_text(size=10),
          axis.title=element_text(size=12)) +
    scale_color_manual(values = colors)


aicsfromfulllong <- reshape(aicsfromfull, direction = 'long',
                            varying = c('dAIC_Jtu', 'dAIC_Jbeta', 'dAIC_Horn'),
                            v.names = 'dAIC',
                            idvar = 'mod',
                            timevar = 'type',
                            times = c('Jtu', 'Jbeta', 'Horn'))

signedsqrt = function(x) sign(x)*sqrt(abs(x))
signedsq = function(x) sign(x) * x^2
newtrans <- trans_new(name = 'signedsqrt', transform = signedsqrt, inverse = signedsq)
aicsfromfulllong$row[aicsfromfulllong$type == 'Jtu'] <- aicsfromfulllong$row[aicsfromfulllong$type == 'Jtu'] + 0.1
aicsfromfulllong$row[aicsfromfulllong$type == 'Horn'] <- aicsfromfulllong$row[aicsfromfulllong$type == 'Horn'] - 0.1

# plot
p3 <- ggplot(aicsfromfulllong, aes(x = dAIC, y = row, group = type, color = type)) +
    scale_x_continuous(trans = newtrans) +
    scale_y_continuous(breaks = nrow(aicsfromfull):1, labels = aicsfromfull$mod) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') + 
    geom_hline(yintercept = nrow(aicsfromfull):1, linetype = 'dashed', color = 'grey', size = 0.2) + 
    geom_point() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(),
          legend.title=element_blank(),
          axis.text=element_text(size=6),
          axis.title=element_text(size=12)) +
    labs(tag = 'C', x = 'Change in AIC', y = '', title = 'Removing terms')


fig2 <- arrangeGrob(p1, p2, p3, 
                    ncol = 2, 
                    layout_matrix = rbind(c(1,2), c(3,3)),
                    heights=c(unit(0.5, "npc"), unit(0.5, "npc")))
ggsave('figures/fig2.png', fig2, width = 6, height = 5, units = 'in')


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
