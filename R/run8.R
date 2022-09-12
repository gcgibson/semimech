library(rjags)
library(R2jags)
library(splines)
library(ggplot2)
library(rjags)
library(R2jags)
library(splines)








library(covidcast)
library(splines)
library(rjags)
library(R2jags)
library(dplyr)
library(lubridate)
args = commandArgs(trailingOnly=TRUE)

forecast_date <- "2022-01-01"#args[1]
cases <- covidcast_signal(
  data_source = "hhs",
  signal = "confirmed_admissions_covid_1d_prop",
  time_type = "day",
  geo_type = "state",as_of = as.Date(forecast_date) )




cases <- cases[cases$time_value >= "2020-07-01",]
cases <- cases %>% group_by(geo_value) %>% mutate(value = zoo::rollmean(value,7,fill=0)) %>% ungroup()
states <- unique(cases$geo_value)
states <- setdiff(states,c("ut","as","pr","hi","vi"))
fit_res <- list()
wave_matrix_total <- matrix(NA,ncol=1000)
current_waves <- matrix(NA,ncol=1000)
for (geo_val in states){

  print (geo_val)
  mydfcovid <- cases[cases$geo_value ==geo_val ,]$value
  smoothed_deriv <- c(0,diff(lowess(mydfcovid/max(mydfcovid,na.rm=T),f=.2)$y))
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

  if ( length(mydfcovid)-tail(wave_starts,1) < 15){
    wave_starts <- wave_starts[1:(length(wave_starts)-1)]
  }
  wave_matrix <-matrix(NA,nrow=length(wave_starts)-1,ncol=1000)

  if (length(wave_starts) > 1){
      for (wave_start in 2:length(wave_starts)){
        to_replace <- mydfcovid[wave_starts[wave_start-1]:wave_starts[wave_start]]
        while (length(to_replace) < 1000){
          to_replace <- c(to_replace,NA)
        }
        wave_matrix[wave_start-1,] <- to_replace
      }


    NAindex <-function(z){ min(which(is.na(z)))}


    current_wave <- round(mydfcovid[wave_starts[length(wave_starts)]:length(mydfcovid)])
    if (is.null(dim(wave_matrix))){
      wave_matrix <- matrix(wave_matrix,nrow=1)
    }
    longest_wave_length <- max(length(current_wave),max(apply(wave_matrix,1,NAindex))) + 30

    #wave_matrix <- matrix(round(wave_matrix[,1:longest_wave_length]),nrow=nrow(wave_matrix))
  }


  length_of_current_wave <- length(current_wave)
  current_wave[tail(1:length(current_wave),3)] <- NA

  while(length(current_wave) < dim(wave_matrix)[2]){
    current_wave <- c(current_wave,NA)
  }



  wave_matrix <- rbind(wave_matrix,matrix(current_wave,nrow=1))
  wave_matrix_total <- rbind(wave_matrix_total,wave_matrix)
  current_waves <- rbind(current_waves,current_wave)

  ###plotting
  # p0 <- ggplot(data=cases[cases$geo_value ==geo_val & cases$time_value <= as.Date(forecast_date)- 5 ,],aes(x=time_value,y=value)) + geom_line()
  #
  # wave1 <- mydfcovid[1:wave_starts[1]]
  # wave2 <- mydfcovid[wave_starts[1]:wave_starts[2]]
  # wave3 <- mydfcovid[wave_starts[2]:wave_starts[3]]
  # wave4 <- mydfcovid[wave_starts[3]:wave_starts[4]]
  # wave5 <- mydfcovid[wave_starts[4]:length(mydfcovid)]
  # wave5[tail(1:length(wave5),4)] <- NA
  #
  #
  # p1 <- ggplot(data=data.frame(wave2=wave2),aes(x=1:length(wave2),y=wave2,col="Wave 1")) + geom_line()+
  # geom_line(data=data.frame(wave3=wave3),aes(x=1:length(wave3),y=wave3,col="Wave 2")) +
  # geom_line(data=data.frame(wave4=wave4),aes(x=1:length(wave4),y=wave4,col="Wave 3")) +
  # geom_line(data=data.frame(wave5=wave5),aes(x=1:(length(wave5)),y=wave5,col="Wave 4")) + theme_minimal() +
  #   geom_line(data=data.frame(y=colMeans(forecast[[1]]),x=1:length(colMeans(forecast[[1]]))),aes(x=x,y=y,col='forecast'))
  #
  # p1
  # p2 <-ggplot(data=data.frame(x=forecast[[2]][,1])) + geom_density(aes(x=x,fill=as.factor(1))) +
  #   geom_density(data=data.frame(x=forecast[[2]][,2]),aes(x=x,fill=as.factor(2))) +
  #   geom_density(data=data.frame(x=forecast[[2]][,3]),aes(x=x,fill=as.factor(3))) +
  #
  #   geom_density(data=data.frame(x=forecast[[3]]),aes(x=x,fill=as.factor(4)) ) + xlim(c(0,.001))
  #
  #
  # library(cowplot)
  #
  # cowplot::plot_grid(p0+theme_minimal()+ ylab("Hospitalizations") + xlab("Date") +
  #                      theme(legend.title=element_blank()),
  #                    p1+ ylab("Hospitalizations") + xlab("Wave Time") +
  #                      theme(legend.title=element_blank()),
  #                    p2 + xlab(expression(theta)) + theme_minimal() +
  #                      theme(legend.title=element_blank()),nrow=3)
  #



}




plot(NA,xlim=c(0,200),ylim=c(0,20))
for (row in seq(1,nrow(wave_matrix_total),4)){
  lines(wave_matrix_total[row,])
}


index_of_maxes <- c()
maxes <- c()
for (row in 1:nrow(wave_matrix_total)){
  index_of_maxes <- c(index_of_maxes,which.max(wave_matrix_total[row,]))
  maxes <- c(maxes,max(wave_matrix_total[row,],na.rm=T))
}





df_to_submit <- data.frame(matrix(NA,ncol=7,nrow=0))
colnames(df_to_submit) <- c("target","location","forecast_date","target_end_date","quantile","value","type")


for (state in states){
  # extract this states df
  state_forecast_df <- fit_res[[state]]
  # sum by weekly (manually so I don't mess up)
  ncol_ <- ncol(state_forecast_df)
  state_forecast_df <- state_forecast_df[,(ncol_-29):ncol_]

  # define quantiles to populate
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

  # extract quantiles
  state_forecast_df_weekly <- apply(state_forecast_df,2,function(x){
    return (quantile(x,probs=quantiles))
  })

  # make dataframe for this state with correct dimensions, number of quantiles (23)
  # by number of weeks ahead (1-4)
  state_df <- data.frame(h=rep(1:30,23),value =c(t(state_forecast_df_weekly)),
                         quantile = rep(quantiles,each=30))

  # assign correct state
  state_df$state <- state

  # modify week ahead to correct formatting
  state_df$target <- unlist(lapply(state_df$h , function(x){
    return (paste0(x, " day ahead inc hosp"))
  }))

  # assign correct forecast date (latest Monday)

  state_df$forecast_date <- as.Date(forecast_date)

  # assign corresponding next saturday modulo week haead as the target date
  state_df$target_end_date <- Reduce(c,lapply(state_df$h , function(x){
    return (as.Date(forecast_date) + x )
  }))

  # assign the correct state
  state_df$location <- toupper(state)

  # add type column
  state_df$type <- unlist(lapply(state_df$quantile , function(x){

    return ("quantile")

  }))
  # select only columns to be submitted
  state_df <- state_df %>% dplyr::select(target,location,forecast_date,target_end_date,quantile,value,type)



  df_to_submit <- rbind(df_to_submit,state_df)
}


##
library(lubridate)
library(readr)
library(stringr)
df_to_submit$date <- df_to_submit$target_end_date
df_to_submit$covid <- df_to_submit$value



##finally add to fips
df_to_submit <- df_to_submit %>% dplyr::select(target,location,forecast_date,target_end_date,quantile,value,type)

fips <- read_csv("data/fips-codes.csv")
colnames(fips) <- c("fips", "state_full", "state", "alphacount")
glimpse(fips)

fips <- fips %>%
  mutate(fips = str_pad(fips, 2, pad = "0"))
fips$location <- fips$state
df_to_submit <- df_to_submit %>% left_join(fips,by="location")
df_to_submit$location <- df_to_submit$fips
df_to_submit <- df_to_submit %>% dplyr::select(target,location,forecast_date,target_end_date,quantile,value,type)
write.csv(df_to_submit,file = paste0(as.Date(forecast_date),"-UT-Ankh7.csv"),row.names = F)


