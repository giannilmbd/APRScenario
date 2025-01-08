## -----------------------------------------------------------------------------
library(APRScenario)

## ----setup, include=FALSE,echo=FALSE------------------------------------------
knitr::opts_chunk$set(cache = TRUE, echo = FALSE, include = TRUE, warning = FALSE, comment = NA, out.width="70%", out.height="70%",fig.align = "center", dpi=100,units="cm" ,fig.height=20, fig.width=20, out.width="120%")

## ----Set-directory,echo=FALSE,include=FALSE-----------------------------------

# setwd('/home/gianni/Dropbox/SNB/HousePrices/scenarios/')
#following is because of the German way of typing dates etc.
Sys.setlocale("LC_TIME", "English.UTF-8") # Windows
# Sys.getlocale()

system('mkdir -p figures')

## ----Load-necessary-packages--------------------------------------------------
# install.packages("~/Dropbox/SNB/APRScenario_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(APRScenario)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(gridExtra)
library(tidyr)
library(bsvars)

## ----DATA---------------------------------------------------------------------
data(NKdata)
X0<-NKdata[,-1] %>% as.data.frame()
varbls<-names(X0)
p<-autoplot(X0 %>% ts(.,frequency=4,start=round(NKdata$year[1])),facets=T)+ylab('')+xlab('')


ggsave('data.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')


## ----Sign-Restriction---------------------------------------------------------

sr <- matrix(NA, nrow = length(varbls), ncol = length(varbls))

# Fill the matrix with sign restrictions as per the previous setup
sr <- matrix(c(
1,-1,-1,# GDP
1,1,-1, # infl
1,NA,1 # interest rate
), nrow = 3, byrow = TRUE)


## ----estimate-BVAR------------------------------------------------------------


# ############# Subset of variables
# 1. Baseline VAR --------
subset_<-c(1:3) # start with subsets of variables
#/////////////////////////////
n_draws=1500 # increase for final estimation
p=4 # lags
# possible subsampling
X=X0[1:(nrow(X0)-4),]
# specify the model
specification  = bsvarSIGNs::specify_bsvarSIGN$new(data=as.matrix(X,3,3),
                                      p        = p,
                                       sign_irf = sr[1:3,1:3])

# estimate the model
posterior      = bsvars::estimate(specification, S = n_draws)

# compute and plot impulse responses
irf            = bsvars::compute_impulse_responses(posterior, horizon = 40)
# {
#   X11()
# plot(irf, probability = 0.68)
# }
{
  p<-plot_bvars(irf, significance_level = 0.32,
                central_tendency = 'median', variable_names = varbls, shock_names = varbls)
  
  # grid.arrange(p)
}

ggsave('irf_bvars.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')

sr1<-sr %>% as.data.frame()
rownames(sr1)<-varbls 
names(sr1)<-varbls
knitr::kable(sr1,row.names = T,caption = 'Sign Restrictions: IRF 1 (shocks in col.)')
# sr1 %>% cat(file='figures/sr_bvars.tex')


f_bvar<-bsvars::forecast(
  posterior,
  horizon = 3,
  exogenous_forecast = NULL,
  conditional_forecast = NULL
)
# VAR structure Y'=X'*B+u' (see specification$data_matrices) 


## ----dpi=100,units='cm', fig.height=10, fig.width=10, out.width="100%"--------

knitr::include_graphics("figures/irf_bvars.png")


## ----Matrices,echo=TRUE-------------------------------------------------------
# Actual function ---------------
gen_mats()

## ----Run-code-for-unconditional-forecast--------------------------------------
{# run all block 
h=4
y_h_all<-forc_h(h,n_sim=200)
y_h<-y_h_all[[1]]
hist_h<-y_h_all[[3]]
b_h<-y_h_all[[4]]
{# extract quantiles
y_h_m<-apply(y_h,c(1,2),FUN=function(x)quantile(x,0.5))
y_h_l<-apply(y_h,c(1,2),FUN=function(x)quantile(x,0.16))
y_h_u<-apply(y_h,c(1,2),FUN=function(x)quantile(x,0.84))


hist_h_l<-apply(hist_h,c(1,2),FUN=function(x)quantile(x,0.16))
hist_h_u<-apply(hist_h,c(1,2),FUN=function(x)quantile(x,0.84))

}
# convert to data frame for better manipulation
{
y_h_m<-as.data.frame(t(y_h_m))
y_h_u<-as.data.frame(t(y_h_u))
y_h_l<-as.data.frame(t(y_h_l))

hist_h_u<-as.data.frame(t(hist_h_u))
hist_h_l<-as.data.frame(t(hist_h_l))


names(y_h_m)<-varbls
y_h_m$h<-1:nrow(y_h_m)
y_h_tot<-pivot_longer(y_h_m,cols=all_of(varbls),names_to='vars',values_to='Median')
names(y_h_l)<-varbls
y_h_l$h<-1:nrow(y_h_l)
y_h_tot<-left_join(y_h_tot,
                   pivot_longer(y_h_l,cols=all_of(varbls),names_to='vars',values_to='LB'),
                   by=c('h','vars'))

names(y_h_u)<-varbls
y_h_u$h<-1:nrow(y_h_u)
y_h_tot<-left_join(y_h_tot,
                   pivot_longer(y_h_u,cols=all_of(varbls),names_to='vars',values_to='UB'),
by=c('h','vars'))

names(hist_h_l)<-varbls
hist_h_l$h<-1:nrow(hist_h_l)
hist_h_tot<- pivot_longer(hist_h_l,cols=all_of(varbls),names_to='vars',values_to='LB_s')

names(hist_h_u)<-varbls
hist_h_u$h<-1:nrow(hist_h_u)
hist_h_tot<-left_join(hist_h_tot,
                   pivot_longer(hist_h_u,cols=all_of(varbls),names_to='vars',values_to='UB_s'),
by=c('h','vars'))


y_h_tot<-left_join(y_h_tot,hist_h_tot,by=c('h','vars'))

}
# inspect result
dt_t<-as.data.frame(t(Y))
dt_t$h=1:nrow(dt_t)
y_data<-pivot_longer(dt_t,cols =all_of(1:n_var),values_to =names(y_h_tot)[3],names_to = "vars" )
y_data<-y_data[y_data$h>=last(y_data$h)-4,]  
  
y_data$h<-y_data$h-max(y_data$h)
y_data$LB<-NA
y_data$UB<-NA
y_data$LB_s<-NA
y_data$UB_s<-NA

y_h_tot<-rbind(y_data,y_h_tot)
y_h_tot<-y_h_tot[order(y_h_tot$h),]
head(y_h_tot)

uncond.forc<-y_h_tot
# plot 


p <- ggplot(uncond.forc[uncond.forc$vars == varbls[1], ], aes(x = h)) +
  # Median line (solid line)
  geom_line(aes(y = Median, color = "med"), linewidth = 1, show.legend = TRUE) +
  
  # Shaded area for 68% HDI
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = "bnd"), 
              alpha = 0.5, show.legend = TRUE) +
  
  # Dashed lines for 68% high credibility band of history
  geom_line(aes(y = LB_s, color = "hist"), 
            linetype = "dashed", linewidth = 0.8, show.legend = TRUE) +
  geom_line(aes(y = UB_s, color = "hist"), 
            linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  
  # Labels and theme
  labs(title = "Unconditional Forecast", x = "h", y = varbls[1]) +
  theme_minimal() +
  
  # Custom legend for colors
  scale_color_manual(
    name = "Legend",
    labels=c('med'="Median",'hist'="no shock unc."),
    values = c("med" = "blue",  'hist'= "red")
  ) +
  
  # # Custom legend for fill
  # scale_fill_manual(
  #   name = "Legend",
  #   labels=c('bnd'="no shock unc."),
  #   values = c('bnd' = "lightblue")
  # )+
  theme(legend.position = 'bottom')+
  theme_minimal()


} # end run all block

ggsave('uncond_forc.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')

## -----------------------------------------------------------------------------
knitr::include_graphics("figures/uncond_forc.png")


## ----Run-Scenario-------------------------------------------------------------
#NB: Y contans the data n_var x T
h=4 # horizon
n_sim=200 # number of shock draws
obs=c(2) # number of observables
pos_cond_vars=obs 
# given the path of the observables
TT=nrow(X0)
path=X0[(TT-h+1):TT,obs]
bvarSign_path=X0[(TT-h+1):TT,]
bvarSign_path[,-obs]<-NA

# give the shocks that are not driving the scenario: NA if all driving
shocks=NA#c(1,2,3)
tmp<-scenarios(h,path,obs,shocks)
mu_eps<-tmp[[1]]
Sigma_eps<-tmp[[2]]
mu_y<-tmp[[3]]
Sigma_y<-tmp[[4]]
big_b<-tmp[[5]]
big_M<-tmp[[6]]

y_h<-SimScen(mu_eps,Sigma_eps,mu_y,Sigma_y,big_b,big_M,n_sim,h)

# convert to data frames of central and HPD
y_h_m<-apply(y_h,c(1),FUN=function(x)quantile(x,0.5))
y_h_l<-apply(y_h,c(1),FUN=function(x)quantile(x,0.16))
y_h_u<-apply(y_h,c(1),FUN=function(x)quantile(x,0.84))

cond.for<-data.frame(center=y_h_m,variable=rep(varbls,h),hor=sort(rep(1:h,n_var)))
cond.for$lower<-y_h_l
cond.for$upper<-y_h_u
cond.for<-cond.for[,c("hor","variable","lower","center","upper")]




# inspect result
dt_t<-as.data.frame(t(Y))
dt_t$h=1:nrow(dt_t)
y_data<-pivot_longer(dt_t,cols =all_of(1:n_var),values_to =names(cond.for)[4],names_to = "variable" )
y_data<-y_data[y_data$h>=last(y_data$h)-4,]  
  
y_data$hor<-y_data$h-max(y_data$h)
y_data$lower<-NA
y_data$upper<-NA
y_data<-y_data[,names(cond.for)]
cond.for<-rbind(y_data,cond.for)
cond.for<-cond.for[order(cond.for$hor),]

cond.for$hist<-1
cond.for[cond.for$hor<=0,'hist']<-0



p<-plot_cond_forc(varbls[1])
ggsave('cond_forc.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')

p<-plot_cond_histo(variable=varbls[1],horizon=1)
ggsave('cond_histo.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')

## -----------------------------------------------------------------------------

knitr::include_graphics("figures/cond_forc.png")
knitr::include_graphics("figures/cond_histo.png")


## ----Compare-conditions-APR-with-BVARSign-------------------------------------
# As discused by Jarocinski 2010 there can be multiple solutions when k<n*h. The APR Moore-Penrose inverse is just one of these

tmp_frc<-bsvars::forecast(
  posterior,
  horizon = 4,
  exogenous_forecast = NULL,
  conditional_forecast = as.matrix(bvarSign_path)
)

{f_data<-tmp_frc$forecasts
rownames(f_data)<-varbls
f_h_m<-apply(f_data,c(1,2),FUN=function(x)quantile(x,0.5)) %>% t() %>% as.data.frame()
f_h_m$hor<-1:h
f_h_l<-apply(f_data,c(1,2),FUN=function(x)quantile(x,0.16)) %>% t() %>% as.data.frame()
f_h_l$hor<-1:h
f_h_u<-apply(f_data,c(1,2),FUN=function(x)quantile(x,0.84)) %>% t() %>% as.data.frame()
f_h_u$hor<-1:h

f_h_m<-pivot_longer(f_h_m,cols=!"hor",names_to='variable',values_to = 'centF')
f_h_l<-pivot_longer(f_h_l,cols=!"hor",names_to='variable',values_to = 'lowF')
f_h_u<-pivot_longer(f_h_u,cols=!"hor",names_to='variable',values_to = 'higF')
}
f_h<-left_join(f_h_l,f_h_m,by=c('hor','variable'))
f_h<-left_join(f_h,f_h_u,by=c('hor','variable'))

#f_h
names(uncond.forc)[1:2]<-c("hor","variable")

test2<-left_join(cond.for,f_h,by=c('hor','variable'))
test2<-left_join(test2,uncond.forc,by=c('hor','variable'))
p<-ggplot(data=test2, aes(x=hor)) +
  ylab('')+xlab('')+labs(title='Comparison of (un-)conditional forecasts')+
  geom_point(aes(y=center, color='APR', shape='APR'),size=3) +
  geom_point(aes(y=centF, color='BVARSign', shape='BVARSign'),size=3) +
  geom_point(aes(y=Median, color='Uncond', shape='Uncond'),size=3) +
  facet_wrap(~variable, scales = 'free_y',nrow = 3) +
  scale_color_manual('', values=c('APR'='red', 'BVARSign'='blue', 'Uncond'='green')) +
  scale_shape_manual('', values=c('APR'=4, 'BVARSign'=1, 'Uncond'=6))+theme_minimal()
ggsave('cond_forc_comp.png',plot=p,device='png',path='figures',width=18,height=16,units = 'cm')

## -----------------------------------------------------------------------------
knitr::include_graphics("figures/cond_forc_comp.png")


## ----Kullback-Leibler measure-------------------------------------------------

q<-KL(Sigma_eps,mu_eps,h)

hist(q,main='KL measure (ref value 0.5)')


