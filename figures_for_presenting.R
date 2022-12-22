require(RColorBrewer)

# Narrower thermal niches at higher temperatures

T <- seq(0,33,length.out=5000) # x-axis
w <- seq(4, 0.1, length.out = 12) # widths of the parabolas
Topt <- rep(2, length(w))
for(i in 2:length(w)) Topt[i] <- Topt[i-1] + sqrt(w[i-1]) + sqrt(w[i]) # find the midpoints so that they touch at the x-axis

# make the parabolas
y <- matrix(data = NA, nrow=length(w), ncol = length(T)) 
for(i in 1:nrow(y)){
    y[i,] <- -(T-Topt[i])^2/w[i] + 1
}
y[y<0] <- NA

# colors and plot parameters
cols <- brewer.pal(11, "RdYlBu")
cols <- colorRampPalette(rev(cols))(length(w))
lwd <- 2 # line width

# plot
png(filename = 'temp/niches.png', res = 300, width = 5, height = 2, units = 'in')
par(mai = c(1,0.1, 0.1, 0.1))
plot(T, y[1,], type='l', col =cols[1], bty = 'l', lwd = lwd, yaxt = 'n', ylab = '')
for(i in 2:nrow(y)){
    lines(T, y[i,], col = cols[i], lwd = lwd)
}

dev.off()