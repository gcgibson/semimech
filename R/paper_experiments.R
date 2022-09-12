library(tidyverse)
library(cowplot)
library(lubridate)
library(mgcv)
library(splines)
library(MASS)
library(dplyr)
library(rstanarm)


last_date <- as.Date("2022-08-28") - 7*4
# For better plotting
mytext <- element_text(angle=90, hjust=0.5, vjust=0.5)
# mytext <- element_text(angle=45, hjust=0.5, vjust=0.5)

theme_set(theme_minimal_grid() +
            theme(axis.text.x = mytext, legend.position = "top"))

byy1 <- function() scale_x_date(date_labels = "%b '%y", date_breaks = "1 month")
byy3 <- function() scale_x_date(date_labels = "%b '%y", date_breaks = "3 month")

# Get the HHS data --------------------------------------------------------

source("R/get_hhs.R")
source("R/semimech2.R")
hhs <- get_hhs()

hhs_sub <- hhs %>%
  dplyr::select(state, date, fips:covid_per_cap_av7_center) %>%
  filter(!is.na(flu) | !is.na(covid))
hhs_sub <- hhs_sub %>% dplyr::select(state,date,covid)

hhs_sub_us <- hhs_sub %>% dplyr::group_by(date) %>% dplyr::summarise(state="US",covid = sum(covid))
hhs_sub <- rbind(hhs_sub,hhs_sub_us)

### Set OR 2022-08-26 to 52
hhs_sub[hhs_sub$state == "OR" & hhs_sub$date == "2022-08-24", ]$covid <- 38
hhs_sub[hhs_sub$state == "OR" & hhs_sub$date == "2022-08-25", ]$covid <- 40
hhs_sub[hhs_sub$state == "OR" & hhs_sub$date >= "2022-08-26", ]$covid <- 42

forecast_days <- 31
start_date <- "2020-07-01"



statevec <- hhs_sub$state %>% unique() %>% sort()

dflist <- vector("list", length=length(statevec))



for (k in 1:length(statevec)) {

  cat(sprintf("Building model for %s (state %i out of %i)...\n",
              statevec[k], k, length(statevec)))

  mydf <- hhs_sub %>%
    # mutate(flu = flu + 1) %>%
    filter(date >= start_date) %>%
    filter(state == statevec[k]) %>%
    filter(date < last_date) %>%
    mutate(covid = ifelse(is.na(covid), 0, covid)) %>%
    mutate(time = as.numeric(date - min(date)))

  plot(mydf$covid,type='l')
  smoothed_deriv <- c(0,diff(lowess(mydf$covid/max(mydf$covid,na.rm=T),f=.15)$y))

  plot(smoothed_deriv,type='l')
  abline(h=0)
  wave_starts <- c()

  for (t in 2:length(smoothed_deriv)){
    if (!is.na(smoothed_deriv[t-1]) & !is.na(smoothed_deriv[t])){
      if (smoothed_deriv[t-1] < 0 & smoothed_deriv[t] >0){
        wave_starts <-c(wave_starts,t)
        abline(v=t)
      }
    }
  }

  wave_matrix <-matrix(NA,nrow=length(wave_starts)-1,ncol=1000)

  for (wave_start in 2:length(wave_starts)){
    to_replace <- mydf$covid[wave_starts[wave_start-1]:wave_starts[wave_start]]
    while (length(to_replace) < 1000){
      to_replace <- c(to_replace,NA)
    }
    wave_matrix[wave_start-1,] <- to_replace
  }


  NAindex <-function(z){ min(which(is.na(z)))}

  longest_wave_length <- max(apply(wave_matrix,1,NAindex))

  wave_matrix <- wave_matrix[,1:longest_wave_length]



  current_wave <- mydf$covid[wave_starts[length(wave_starts)]:length(mydf$covid)]


  mydf_oos <- mydf %>%
    tail(1) %>%
    slice(rep(1, forecast_days)) %>%
    dplyr::mutate(date = date + row_number(),
                  time = time + row_number())

  mydf2 <- rbind(mydf %>% mutate(period="calibration"),
                 mydf_oos %>% mutate(period="forecast"))


  while(length(current_wave) < dim(wave_matrix)[2]){
    current_wave <- c(current_wave,NA)
  }

  wave_matrix <- rbind(wave_matrix,matrix(current_wave,nrow=1))
  forecast <- run_semi_mech(season_past_matrix_hosp = wave_matrix,season_current_hosp =current_wave )

  plot(colMeans(forecast[[1]]))


  plot(1:length(current_wave),current_wave,type='l',xlim=c(0,length(current_wave)+30))
  #spaghetti <- ggplot(data=data.frame(x=1:length(current_wave),y=current_wave)) + geom_line(aes(x=x,y=y,col="mean"))
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


  col_vals = forecast[[3]][1:10]*10e2
  o_cols = order( col_vals)
  colors_for_plot <-heat.colors(100)[round( col_vals )]

  for (row in 1:1000){
    tmp_row <- forecast[[1]][row,]
    lines(seq(1,length(current_wave)+30),tmp_row,col='red')#,col=colors_for_plot[row])
  }
  legend(x = 10,y=60,forecast[[3]][o_cols]*1e2,bty="n", col=colors_for_plot[o_cols],pch=1,cex=.3)

  plot(forecast[[3]],current_wave-forecast[[1]][,100],type='p')
  plot(forecast[[3]],forecast[[1]][,30])
  cor_res <- c()

  for (j in 1:ncol(forecast[[1]])){
    cor_test <- cor.test(forecast[[3]],forecast[[1]][,j])
    cor_res <- c(cor_res,cor_test$estimate)
  }

  plot(cor_res,type='l')
  lines(colMeans(forecast[[1]]))
  #spaghetti
  ### also not really neccessary


  dflist[[k]] <- forecast[[1]]

}

upper_q <- function(x){quantile(x,probs=c(.975))}
lower_q <- function(x){quantile(x,probs=c(.025))}

dflist_as_df <- list()

for (k in 1:length(statevec)) {
  forecast_dates <- seq(as.Date(tail(mydf$date,1)),as.Date(tail(mydf$date,1))+29,by="day")
  forecast_df <- data.frame(date=forecast_dates,pred=colMeans(dflist[[k]]),
                            upper_95 = apply(dflist[[k]],2,upper_q),
                            lower_95 = apply(dflist[[k]],2,lower_q),state=statevec[k])


  dflist_as_df[[k]] <- forecast_df
}


dfall <- dflist_as_df %>%
  plyr::rbind.fill() %>%
  #dplyr::mutate(
   # covid = ifelse(period=="forecast", NA, covid)) %>%
  # arrange(desc(pop)) %>%
  dplyr::mutate(state = factor(state, levels=unique(state))) %>%
  dplyr::arrange(state, date)

myfun <- function(x, pop) x / pop * 1e5
myfun2 <- function(x, pop) x / pop * 1e5

# Number of admissions
covid_hosps <- dfall %>%
  filter(date >= "2022-01-01") %>%
  # filter(state %in% mystates) %>%
  # filter(state %in% c("NC", "IA", "AR", "KS", "NE", "MT", "ND")) %>%
  # filter(state == "FL") %>%
  ggplot() +
  geom_ribbon(aes(date, ymin= lower_95, ymax=upper_95), alpha=0.5,col='blue') +
  # geom_line(aes(date, covid_av7, lty="covid")) +
  #geom_line(aes(date, flu_av7_center), lty="dashed") +
  geom_line(aes(date, pred)) +
  geom_point(data=hhs_sub[hhs_sub$state == "NY",] %>%  filter(date >= "2022-01-01"),aes(date, covid), alpha=0.25) +
  facet_wrap(~state,scales="free")  +
  # scale_y_log10() +
  byy1() +
  scale_x_date(date_labels = "%b '%y") +
  labs(x="", y="covid admissions") + scale_x_date(date_breaks = '1 month', date_labels = '%b') +
 # geofacet::facet_geo(~state,scales='free') +
  theme(panel.background = element_rect(fill='white', colour='white'),
        plot.background =element_rect(fill='white', colour='white') )

wave_cutoff_dates <- mydf$date[wave_starts]

ggsave(filename = "covid_hosps.png",covid_hosps,device = "png",height = 10,width = 16)


wave_plot_df <- data.frame(val = c(t(wave_matrix)),t =rep(1:ncol(wave_matrix),nrow(wave_matrix)),
                           id = rep(1:nrow(wave_matrix),each=ncol(wave_matrix)))
p1 <- ggplot(wave_plot_df,aes(x=t,y=val,group=id,col=as.factor(id))) + geom_line(alpha=.5) + geom_line( data=data.frame(x=1:length(current_wave),y=current_wave,id=4),aes(x=x,y=y,group=as.factor(id),col=as.factor(id))) +
  geom_line(data=data.frame(y=colMeans(forecast[[1]]),id="Forecast"),aes(x=1:(length(colMeans(forecast[[1]]))),y=y,group=as.factor(id)))

p1 + theme_bw()
phi_post <- forecast[[2]]

p2 <-ggplot(data=data.frame(x=forecast[[2]][,1])) + geom_density(aes(x=x,fill=as.factor(1))) +
  geom_density(data=data.frame(x=forecast[[2]][,2]),aes(x=x,fill=as.factor(2))) +
  geom_density(data=data.frame(x=forecast[[2]][,3]),aes(x=x,fill=as.factor(3))) +

  geom_density(data=data.frame(x=forecast[[3]]),aes(x=x,fill=as.factor(4)) )

print (colMeans(forecast[[3]]))

library(cowplot)

cowplot::plot_grid(p1 + theme_bw() + ylab("Hosp") +xlab("Wave Time"),p2+ theme_bw() + xlab("Theta"),nrow=2)

