######## Example plots function-on-function relationship #########

# libraries
library(ggplot2)
library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

# TODO
# * flipppen
# * nicht time
# * conc, hist, full, banded hist, lag, banded

nr_timepoints <- 8
s <- t <- 1:nr_timepoints
this_data <- expand.grid(s = s, t = t)

# plot function
relplot <- function(fill_fun, 
                    main="", 
                    breaks = 1:nr_timepoints,
                    expand = c(0,0),
                    this_widths = c(1,3),
                    this_heights = c(3,1),
                    magic_nr = 0.2
                    ){
  
  # the plot combiner
  paperPlot <- function(plotUpperLeft,
                        plotUpperRight,
                        plotLowerRight, 
                        widths = this_widths, 
                        heights = this_heights)
  {
    

    
    isULcomb <- class(plotUpperLeft)[1]!="gg"
    
    gr1 <- ggplotGrob(plotLowerRight) 
    gr0 <- ggplotGrob(plotUpperRight)
    gr2 <- if(!isULcomb) ggplotGrob(plotUpperLeft) else plotUpperLeft
    
    gr1$widths <- gr0$widths
    if(!isULcomb) gr2$heights <- gr0$heights else{ 
      
      for(i in 1:length(gr2$grobs)) gr2$grobs[[i]]$heights <- gr0$heights
      
    }
    
    #grid.newpage()
    #grid.draw(
    gtable_matrix("demo",matrix(list(gr2,nullGrob(),gr0,gr1),nrow=2), 
                  widths = unit(widths,"null"), heights = unit(heights,"null"))
    #)
    
  }
  
  # the actual plot
  
  this_data$z <- fill_fun(this_data$s, this_data$t) * 3
  # times 3 for different colour
  this_data <- this_data[this_data$z != 0,]
  
  gg_main <- ggplot(this_data, aes(x = s, y = t, fill = z)) + 
    geom_tile(colour = "black") + 
    theme_bw() + 
    ggtitle(main) + 
#    xlab("Covariate time points s") + 
    # ylab("Response time points t") + 
    scale_y_continuous(breaks = breaks, 
                       expand = expand
                       ) +
    scale_x_continuous(expand = expand
                       ) + 
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
  
  
  gglower <- ggplot(data.frame(s=1:nr_timepoints, x=sin(1:nr_timepoints)), aes(x=s, y=x)) + 
    geom_line() + geom_point() + theme_bw() + xlab("Covariate function (time s)") + ylab("") + 
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    ) + scale_x_continuous(breaks = breaks,
                           limits = c(1-magic_nr,nr_timepoints+magic_nr))
  
  ggleft <- ggplot(data.frame(t=1:nr_timepoints, y=tanh(1:nr_timepoints - 5)), aes(x=y, y=t)) + 
    geom_line() + geom_point() + theme_bw() + ylab("Response function (time t)")  + xlab("") + 
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    ) + scale_y_continuous(breaks = breaks,
                           limits = c(1-magic_nr,nr_timepoints+magic_nr))
  
  pp <- paperPlot(ggleft, gg_main, gglower)
  
  return(pp)
  
}

grid.newpage()
grid.draw(ggc <- relplot(fill_fun = function(s,t) s==t, main = "Concurrent Model"))
grid.newpage()
grid.draw(ggh <- relplot(fill_fun = function(s,t) t<=s, main = "Historical Model"))
grid.newpage()
grid.draw(ggl <- relplot(fill_fun = function(s,t) t<=s & t+3>=s, main = "Lag-Lead Model"))
grid.newpage()
grid.draw(ggf <- relplot(fill_fun = function(s,t) TRUE, main = "Full"))
grid.newpage()
grid.draw(ggb <- relplot(fill_fun = function(s,t) t-2<=s & t+2>=s, main = "Band Model"))
