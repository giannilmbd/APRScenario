#' plot_cond_forc function; Data should conatain the variable "variable", the "hor" horizon and a "history"
#'
#' @param varbl2plot name of variable to be plotted (string)
#' @param y_h_cond conditional forecast data frame (eg from SimScen) with names c("hor","variable","lower","center","upper"): hor is a Date object
#' @param center (optional, default 0.5) quantile of the mid value
#' @param lower (optional, default 0.16) quantile of lower range
#' @param upper (optional, default 0.84) quantile of upper range
#' @param T.start start date of the forecast
#' @param T.end end of the forecast
#' @param before (integer: optional) periods of data in the plot: default 8 periods
#' @param freq (optional, default 'quarter') frequency of the data (eg 'quarter' or 'month')
#' @param y_data Data used in the estimation eg t(specification$get_data_matrices()$Y) %>% as.data.frame(); true_data$hor=dates
#' @returns list of plot and data used
#' @examples
#' \dontrun{
#' # Example usage with conditional forecast data
#' plot_result <- plot_cond_forc(varbl2plot = "GDP", 
#'                               y_h_cond = forecast_data,
#'                               T.start = as.Date("2023-01-01"),
#'                               T.end = as.Date("2024-01-01"),
#'                               y_data = historical_data)
#' }
#'
#' @export
#'
#' @import dplyr
#' @importFrom lubridate %m-%
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_line ylab geom_ribbon xlab scale_linetype_manual geom_vline theme_minimal theme
#'
plot_cond_forc<-function(varbl2plot=NULL,y_h_cond=NULL,center=0.5,
                         lower=0.16,upper=0.84,T.start=NULL,T.end=NULL,
                         before = 8,freq='quarter',y_data=NULL){
  dates_date<-seq.Date(from=as.Date(y_data$hor[1]),to=as.Date(T.end),by=freq)
  periods<-dates_date[dates_date>=T.start&dates_date<=T.end]
  if(freq=='quarter'){
    multip<-3
  }else{
    multip<-1
  }



    T.before<-periods[1]%m-%months(multip*before)

  y_h_m<-apply(y_h_cond,c(1),FUN=function(x)quantile(x,0.5)) %>% cond2df(.,name = 'center')
  y_h_l<-apply(y_h_cond,c(1),FUN=function(x)quantile(x,0.16))%>% cond2df(.,name = 'lower')
  y_h_u<-apply(y_h_cond,c(1),FUN=function(x)quantile(x,0.84))%>% cond2df(.,name = 'upper')

  cond_for<-full_join(y_h_m,y_h_l,by=c('variable','period'))
  cond_for<-full_join(cond_for,y_h_u,by=c('variable','period'))
  # names(cond_for)
  cond_for<-cond_for %>% group_by(variable) %>%
    mutate(hor=periods)
  cond.for<-cond_for[,c("hor","variable","lower","center","upper")]

  y_data.l<-pivot_longer(y_data,cols =!hor,values_to ='center',names_to = "variable" )
  y_data.l$lower<-NA
  y_data.l$upper<-NA
  y_data.l<-y_data.l[,c("hor","variable","lower","center","upper")]
  cond.for<-rbind(y_data.l[y_data.l$hor<first(periods),],cond.for) %>% as.data.frame()
  cond.for<-cond.for[order(cond.for$hor,cond.for$variable),]
  true_data<-y_data.l[,c('hor','variable','center')]
  names(true_data)[3]<-'data'
  cond.for<-full_join(cond.for,true_data,by=c('hor','variable'))
  cond.for$forc<-'f'
  cond.for[cond.for$hor<T.start,'forc']<-'d'
  cond.for[cond.for$hor<T.start,'center']<-NA


  p<-ggplot(data=cond.for[cond.for$variable==varbl2plot&cond.for$hor>=T.before,],aes(x=hor))+
    geom_line(aes(y=center,linetype='dashed'))+
    geom_line(aes(y=data,linetype='solid'),color='blue')+ylab(varbl2plot)+
    geom_ribbon(aes(ymax=upper,ymin=lower),fill='pink',alpha=0.5)+xlab('')+
    scale_linetype_manual('',values = c("dashed"="dashed","solid"="solid"),
                          labels=c("dashed"="forecast","solid"="data"))+
    geom_vline(xintercept =as.Date(T.start),linetype='dotted',color='red' )+
    theme_minimal()+theme(legend.position = 'bottom')
  return(list(p,cond.for))
}
cond2df<-function(named_vec=NULL,name='Value'){
  df <- data.frame(
    cname = names(named_vec),
    value = as.numeric(named_vec)
  )
  names(df)<-c('cname',name)
  # Split into variable and period using regex
  df$variable <- sub("\\.[0-9]+$", "", df$cname)
  df$period <- as.integer(sub(".*\\.", "", df$cname))

  # Reorder columns
  df <- df[, c("variable", "period", name)]

  # Optional: sort by variable then period
  df <- df[order(df$variable, df$period), ]

  return(df)
}
