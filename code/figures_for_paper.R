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
## Figure 2: interactions
#########################
newdat <- fread('temp/interactions.csv')

intstoplot <- c('tempave_metab', 'microclim', 'npp', 'human_bowler', 'human_bowler')
realms <- c(NA, NA, NA, 'Marine', 'TerrFresh')
logs <- c(FALSE, FALSE, TRUE, FALSE, FALSE)
tags <- c('A', 'B', 'C', 'D', 'E')
# prep the plots
intplots <- vector('list', length(intstoplot))
for(j in 1:length(intstoplot)){
    subs <- newdat$var == intstoplot[j] & newdat$temptrend > 0 # select warming side
    xvar <- 'temptrend_abs'
    title <- intstoplot[j]
    if(intstoplot[j] %in% c('tsign')){
        subs <- newdat$var == intstoplot[j]
    } 
    if(intstoplot[j] %in% c('thermal_bias')){
        subs <- newdat$var == intstoplot[j]
        xvar <- 'temptrend'
    } 
    if(intstoplot[j] %in% c('human_bowler')){
        subs <- newdat$var == intstoplot[j] & newdat$temptrend > 0 & newdat$REALM2 == realms[j]
        title <- paste0('human:', realms[j])
    } 
    
    thisplot <- ggplot(newdat[subs, ], 
                       aes_string(x = xvar, y = 'preds', 
                                  group = intstoplot[j], 
                                  color = intstoplot[j])) +
        geom_line() +
        coord_cartesian(ylim = c(-0.6, 0.6)) +
        theme(plot.margin = unit(c(0.5,0,0.5,0), 'cm')) +
        labs(title = title, y = 'Slope', tag = tags[j]) +
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

fig2 <- arrangeGrob(grobs = intplots, ncol = 3)

ggsave('figures/fig2.png', fig2, width = 10, height = 6, units = 'in')
