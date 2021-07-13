#' R squared
#'
#' @param x vector of x values
#' @param y vector of y values
#'
#' @return the R squared between the x and y variables
#' @export
rsq <- function(x, y) {
  ind <- which(is.infinite(x) | is.infinite(y))
  if (length(ind)>0){
    x <- x[-ind]
    y <- y[-ind]
  }
  if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3) && (length(unique(x[is.na(x) ==F]))>1)){
    rs <- summary(lm(y~x))$r.squared
  }else(rs <- NA)
  rs
}

#' #' Intercept, slope of the linear model fit and Shapiro-Wilk normality test of the residuals
#' #'
#' #' @param x vector of x values
#' #' @param y vector of y values
#' #'
#' #' @return the intercept and slope of the linear fit and p value of Shapiro-Wilk normality test
#' #' @export
#' linFit <- function(x, y) {
#'   ind <- which(is.infinite(x) | is.infinite(y))
#'   if (length(ind)>0){
#'     x <- x[-ind]
#'     y <- y[-ind]
#'   }
#'   if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3) && (length(unique(x[is.na(x) ==F]))>1)){
#'     mod <- lm(y~x)
#'     coef <- summary(mod)$coefficients
#'     tst <- shapiro.test(mod$residuals)#p-value > 0.05 implies that the distribution of the data is not significantly different from normal distribution
#'     out <- c(coef[1,1], coef[2,1], tst$p.value)# intercept and slope
#'   }else(out <- c(NA,NA,NA))
#'   out
#' }


#' RMSE
#'
#' @param val vector of x values
#' @param meas vector of y values
#'
#' @return the RMSE of the two vectors
#' @export
rmse <- function(val, meas){
  sqrt(mean((val - meas)^2,na.rm=TRUE))
}

#' MAPE
#'
#' @param val vector of x values
#' @param meas vector of y values
#'
#' @return the MAPE of the two vectors
#' @export
mape <- function(val, meas){
  mean(abs(val - meas)/abs(val),na.rm=TRUE)
}

#' Derive performance (R2, MAPE or RMSE) from recovery indicators
#'
#' @param val  ground truth
#' @param meas measured
#' @param sttngs simulation settings file
#' @param recSttngs recovery settings file
#' @param metr recovery metric being evaluated
#' @param perf performance indicator
#'
#' @return list of performance indicator
#' @export
#' @import reshape2
#' @import plyr
#' # calcPerf <- function(val, meas, sttngs, recSttngs, metr, perf){#   lst <- list()#   simcases <- names(meas)##   for(sci in 1:length(meas)){# evaluated parameters#     vls <- list()#     for(rpi in 1:length(meas[[1]])){# recovery settings#       if((metr == 'SL') & ((recSttngs$input[rpi] == 'raw') | (recSttngs$input[rpi] == 'smooth'))){#       }else{#         val[[sci]][[rpi]][is.infinite(val[[sci]][[rpi]])] <- NA#         meas[[sci]][[rpi]][is.infinite(meas[[sci]][[rpi]])]<-NA#         if(perf == 'R2'){#           vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) rsq(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))#         }#         if(perf == 'RMSE'){#           vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) rmse(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))#         }#         if(perf == 'MAPE'){#           vl <- sapply(1:dim(meas[[sci]][[rpi]])[2], function(ii) mape(val[[sci]][[rpi]][,ii], meas[[sci]][[rpi]][,ii]))#         }#         tmp <- data.frame(t(vl))#         names(tmp) <- sttngs[[simcases[sci]]]#         tmp2 <- melt(tmp)#         tmp2$Metric <- factor(metr)#         tmp2$Dense <- factor(recSttngs$freq[rpi])#         tmp2$Smooth <- factor(recSttngs$input[rpi])#revalue(factor(recSttngs$input[rpi]), c("BFAST"="segmented", 'smooth'='smoothed'))#         tmp2$Period <- factor(recSttngs$nDist[rpi])#revalue(factor(recSttngs$nDist[rpi]), c("1"="Short", "12"="Long"))#         vls <- rbind(vls,tmp2)#       }#     }#     vls$Period <- revalue(vls$Period, c("1"="Short", "12"="Long"))##     levels(vls$Period) <- c('Short', 'Long')#     vls$Period[vls$Dense == 'annual'] <- "Long"#     lst[[names(meas)[sci]]] <- vls#   }#   lst#

#' Plot results sensitivity analysis
#'
#' @param data data to be plotted, should be a dataframe with the following headers: Metric, Dense, Smooth, Period, Breaks, Seas, variable, value
#' @param lbls legend labels
#' @param xlbl title x axis
#' @param ylbl title y axis
#' @param scales should the x and y axis of the subplots have fixed ranges? Can be set to 'fixed', 'free_x' and 'free_y'
#'
#' @return ggplot object
#' @import ggplot2
#' @import colorspace
#' @export
#'
plotSens  <- function(data,  xlbl, ylbl, scales = 'fixed'){
  ggplot(data, aes(variable,value,color=interaction(Dense), group = interaction(Smooth,Dense))) +
    geom_line(aes(linetype=Smooth),size=3, alpha = 1)+#linetype=interaction(Dense,Smooth)+
    # geom_point(aes(shape=Smooth),size=1.2)+
    scale_color_discrete_qualitative(palette = 'Dark 3',  name = 'Noise removal')+#
    # ,labels=lbls
    facet_grid(vars(Metric),vars(Period), scales = scales)+
    # scale_y_continuous(trans='log2')+
    # labs(color = "Preprocessing")+
    scale_y_continuous(trans='log2',labels=scaleFUN)+
    xlab(xlbl) +
    ylab(ylbl)+
    theme_light(base_size = 50)+
    guides(linetype = guide_legend(title = 'Noise removal'),
           colour=guide_legend(title="Temporal aggregation"))+
             theme(legend.key.width = unit(3, 'cm'))
}

#' Plot results sensitivity analysis for the environmental parameters
#'
#' @param data data to be plotted, should be a dataframe with the following headers: Metric, Dense, Smooth, Period, Breaks, Seas, variable, value, param
#' @param xlbl title x axis
#' @param ylbl title y axis
#' @param scales should the x and y axis of the subplots have fixed ranges? Can be set to 'fixed', 'free_x' and 'free_y'
#'
#' @return ggplot object
#' @import ggplot2
#' @import colorspace
#' @export
#'
plotEnv  <- function(data, xlbl, ylbl, scales = 'fixed'){
  ggplot(data, aes(variable,value,color=param, group = param)) +
    geom_line(aes(),size=1.2, alpha = 1)+#linetype=interaction(Dense,Smooth)+
    geom_ribbon(aes(ymin=data$min_conf, ymax=data$max_conf), linetype=2, alpha=0.1)+
    scale_color_discrete_qualitative(palette = 'Dark 3')+
    facet_grid(vars(Metric), vars(paramType), scales = scales)+
    labs(color = "Parameter")+
    xlab(xlbl) +
    ylab(ylbl)+
    theme(axis.text.x = element_text(color = "grey50", size = 20),
          axis.text.y = element_text(color = "grey50", size = 20),
          axis.title.x = element_text(color = "grey20", size = 25),
          axis.title.y = element_text(color = "grey20", size = 25),
          plot.title = element_text(size=25),
          legend.title = element_text(size=25),
          legend.text = element_text(color = "grey50",size=25),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20,color = "grey20"))
}

#' Compare the reliability of the recovery metrics
#'
#' @param data data to be plotted, should be a dataframe with the following headers: Metric, value
#' @param xlbl title x axis
#' @param ylbl title y axis
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import colorspace
#'
plotMet <- function(data, xlbl, ylbl){
  ggplot(data, aes(Metric, value,color=Metric))+
    geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=T, lwd=2.5)+
    scale_color_discrete_qualitative(palette = 'Dark 3', name = 'Metric', guide=FALSE)+
    # + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    scale_y_continuous(trans='log2',labels=scaleFUN)+
    xlab(xlbl) +
    ylab(ylbl)+
    theme_light(base_size = 50)
}

#' Round number to 3 digits
#'
#' @return rounded number
#' @export
#'
scaleFUN <- function(x){
  formatC(x, format = "e", digits = 2)#format(x, scientific = TRUE)#sprintf(paste0("%.3f"), x)
}


#' plot the effect of each parameter on the performance using box plots
#'
#' @param data data to be plotted, should be a dataframe with the following headers: variable, value, param
#' @param xlbl x label
#' @param ylbl y label
#' @param scales should the x and y axis of the subplots have fixed ranges? Can be set to 'fixed', 'free_x' and 'free_y'
#'
#' @return ggplot plot object
#' @import colorspace
#' @export
#'
plotSensBar <- function(data, xlbl, ylbl, scales = 'free_y'){
  ggplot(data, aes(variable,value,color=param, group = interaction(variable,param)))+
    geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=T, alpha = 0.8)+
    stat_summary(aes(group = param),fun = median, geom = "point", size = 3, position = position_dodge(width = .75)) +
    stat_summary(aes(group = param),fun = median, geom = "line", size = 5, position = position_dodge(width = .75))+
    scale_color_discrete_qualitative(palette = 'Dark 3', name = 'Parameter')+#, guide=FALSE
    facet_grid(vars(Metric), vars(paramType), scales = scales)+
    scale_y_continuous(trans='log2',labels=scaleFUN)+
    xlab(xlbl) +
    ylab(ylbl)+
    # theme(axis.text.x = element_text(color = "grey50", size = 20),
    #       axis.text.y = element_text(color = "grey50", size = 20),
    #       axis.title.x = element_text(color = "grey20", size = 25),
    #       axis.title.y = element_text(color = "grey20", size = 25),
    #       plot.title = element_text(size=25),
    #       legend.title = element_text(size=25),
    #       legend.text = element_text(color = "grey50",size=25),
    #       strip.text.x = element_text(size = 20),
    #       strip.text.y = element_text(size = 20,color = "grey20"))
  theme_light(base_size = 50)
}

#' Plot the effect of the preprocessing techniques on the performance of the recovery indicators
#'
#' @param data data to be plotted, should be a dataframe with the following headers: Dense, Smooth, Metric, Period
#' @param xlbl x label
#' @param ylbl y label
#' @param scales should the x and y axis of the subplots have fixed ranges? Can be set to 'fixed', 'free_x' and 'free_y'
#'
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import colorspace
#'
pltPrepBox <- function(data, xlbl, ylbl, scales = 'free_y'){
  ggplot(data, aes(interaction(Dense),value,color=interaction(Smooth))) +
    geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=2, notch=F, lwd=2.5)+#linetype=interaction(Dense,Smooth)+
    facet_grid(vars(Metric), vars(Period), scales = scales)+
    scale_y_continuous(trans='log2',labels=scaleFUN)+
    xlab(xlbl) +
    ylab(ylbl)+
    # theme(axis.text.x = element_text(color = "grey50", size = 20),
    #       axis.text.y = element_text(color = "grey50", size = 20),
    #       axis.title.x = element_text(color = "grey20", size = 25),
    #       axis.title.y = element_text(color = "grey20", size = 25),
    #       plot.title = element_text(size=25),
    #       legend.title = element_text(size=25),
    #       legend.text = element_text(color = "grey50",size=25),
    #       strip.text.x = element_text(size = 20),
    #       strip.text.y = element_text(size = 20,color = "grey20"))+
    theme_light(base_size = 50)+
    guides(colour=guide_legend(title="Noise removal"))}
