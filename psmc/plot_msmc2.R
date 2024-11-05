#this is a function to show the ticks in the plot
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  labels <- sapply(major.ticks,function(i)
    as.expression(bquote(10^ .(i)))
  )
  axis(ax,at=major.ticks,labels=labels,...)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

# here you group individuals per population
bel=c('GC129388','GC129395')
fra=c('GC129400','GC129406')
por=c('GC129413','GC129417')
spa=c('GC136116','GC136110')
hei=c('GC136123','GC136117')

#and here you group them per ecotype
tidal=c('GC129388','GC129400','GC129413','GC136116','GC136123')
seasonal=c('GC129395','GC129406','GC129417','GC136110','GC136117')

pops=list(bel,fra,por,spa,hei)
popNames=c('Belgium',
           'France',
           'Portugal',
           'Spain',
           'Heist')
color = c('hotpink2','mediumseagreen','steelblue2','orange','purple')

linetype=c(1,2,3,4,5)

par(mfrow=c(4,4),mai=c(0.05,0,0.05,0.05), oma=c(5,5,1,3)+0)

layout(matrix(c(1:6), ncol=2)) 

for(i in 1:length(pops)){
  plot(0, pch = "", ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE,  xlim=c(5000,2000000),ylim=c(5,8))
  
  for(name in pops[[i]]){
    fileN <- read.table(paste(name,'.final.Rin.txt', sep=""),h=T)
    par(new=TRUE)
    plot(fileN$leftX, log10(fileN$popSize), type = "S", col=color[[i]],yaxt='n',xaxt='n',ylab="",xlab="",log="x",  xlim=c(5000,2000000),ylim=c(5,8))
    text(100000,7.6,substitute(paste(italic(nn)), list(nn=popNames[i])))
    
  }
  if(i%%4 ==0){
    at.x <- outer(1:9, 10^(0:7))
    lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
    axis(1, at=at.x, labels=lab.x, las=1)
    mtext(side = 1, text = expression(paste("Years (g=1, µ=2.1x10"^"-9",")")), cex=0.7, line=3) 
  }
  if(i<=4){
    minor.ticks.axis(2,9,mn=0,mx=8)
    mtext(side = 2, text = expression(paste("Ne")), cex=0.7, line =3) 
  }
}


# Set up margins, outer margins, and axis labels
par(mai=c(0.7, 0.7, 0.2, 0.2), oma=c(5,5,1,3)+0)

# Initialize the plot with the desired axes, no data points yet
plot(0, pch = "", ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE,  
     xlim=c(5000,2000000), ylim=c(5,8), log="x")

# Customize x-axis and y-axis labels
at.x <- outer(1:9, 10^(0:7))
lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
axis(1, at=at.x, labels=lab.x, las=1)
mtext(side = 1, text = expression(paste("Years (g=1, µ=2.1x10"^"-9",")")), cex=0.7, line=3) 

minor.ticks.axis(2,9,mn=0,mx=8)
mtext(side = 2, text = expression(paste("Ne")), cex=0.7, line=3)

# Loop through each sample to plot them all in the same graph
for(i in 1:length(pops)){
  for(name in pops[[i]]){
    fileN <- read.table(paste(name, '.final.Rin.txt', sep=""), h=T)
    par(new=TRUE)  # Allow adding new plots to the existing graph
    # Use linetype from the linetype vector
    plot(fileN$leftX, log10(fileN$popSize), type = "S", col=color[[i]], lty=linetype[i], yaxt='n', xaxt='n', 
         ylab="", xlab="", log="x", xlim=c(5000,2000000), ylim=c(5,8))
  }
}

# add a legend
legend("topright", legend=popNames, col=color, lty=linetype, cex=0.8)


# Set up margins, outer margins, and axis labels
par(mfrow=c(1,1),mai=c(0.05,0,0.05,0.05), oma=c(5,5,1,3)+0)

layout(matrix(c(1:6), ncol=2))

# Loop through each sample to plot them all in the same graph
for(i in 1:length(pops)){
  plot(0, pch = "", ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE,  xlim=c(5000,2000000),ylim=c(5,8))
  
  for(name in pops[[i]]){
    fileN <- read.table(paste(name,'.final.Rin.txt', sep=""),h=T)
    
    # Determine the line type based on ecotype
    if(name %in% tidal) {
      line_type <- 1  # Solid line type for tidal (short-winged)
    } else if(name %in% seasonal) {
      line_type <- 2  # Dashed line type for seasonal (long-winged)
    }
    
    par(new=TRUE)
    plot(fileN$leftX, log10(fileN$popSize), type = "S", col=color[[i]], lty=line_type, yaxt='n', xaxt='n', ylab="", xlab="", log="x", xlim=c(5000,2000000), ylim=c(5,8))
    text(100000,7.6,substitute(paste(italic(nn)), list(nn=popNames[i])))
  }
  
  if(i%%4 == 0){
    at.x <- outer(1:9, 10^(0:7))
    lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
    axis(1, at=at.x, labels=lab.x, las=1)
    mtext(side = 1, text = expression(paste("Years (g=1, µ=2.1x10"^"-9",")")), cex=0.7, line=3) 
  }
  
  if(i <= 4){
    minor.ticks.axis(2,9,mn=0,mx=8)
    mtext(side = 2, text = expression(paste("Ne")), cex=0.7, line = 3) 
  }
}

# Set up one plot area for all samples
par(mai=c(0.9, 0.9, 0.5, 0.1))  # Adjust margins for a single plot

# Initialize the plot without points to set the axes
plot(0, pch = "", ylab = "Ne", xlab = expression(paste("Years (g=1, µ=2.1x10"^"-9",")")), 
     xlim = c(5000, 2000000), ylim = c(5, 8), log = "x", xaxt = "n", yaxt = "n", bty = "n")

# Add log-scaled axis labels
at.x <- outer(1:9, 10^(0:7))
lab.x <- ifelse(log10(at.x) %% 1 == 0, at.x, NA)
axis(1, at = at.x, labels = lab.x, las = 1)

# Add y-axis with minor ticks
minor.ticks.axis(2, 9, mn = 0, mx = 8)

# Loop through each population and plot them all on the same figure
for(i in 1:length(pops)){
  for(name in pops[[i]]){
    fileN <- read.table(paste(name,'.final.Rin.txt', sep=""),h=T)
    
    # Determine the line type based on ecotype
    if(name %in% tidal) {
      line_type <- 1  # Solid line type for tidal (short-winged)
    } else if(name %in% seasonal) {
      line_type <- 2  # Dashed line type for seasonal (long-winged)
    }
    
    # Plot the population data
    lines(fileN$leftX, log10(fileN$popSize), type = "S", col = color[[i]], lty = line_type, xlim = c(5000, 2000000), ylim = c(5, 8), log = "x")
  }
}

# add a legend
legend("topright", legend=popNames, col=color, lty=line_type, cex=0.8)
