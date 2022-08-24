library(covidcast)
library(splines)
library(rjags)
library(R2jags)
library(dplyr)
library(lubridate)
args = commandArgs(trailingOnly=TRUE)

forecast_date <- args[1]
cases <- covidcast_signal(
  data_source = "hhs",
  signal = "confirmed_admissions_covid_1d",
  time_type = "day",
  geo_type = "state",as_of = as.Date(forecast_date) )

cases <- cases[cases$time_value >= "2020-08-01",]
cases <- cases %>% group_by(geo_value) %>% mutate(value = zoo::rollmean(value,7,fill=0)) %>% ungroup()
states <- unique(cases$geo_value)
states <- setdiff(states,c("as","pr","hi","vi"))
fit_res <- list()
for (geo_val in states){

  print (geo_val)
  model_code <- "
  model
  {

     spline[1] <- exp(inprod(basis_2[1,], beta_w2[1:N_knots2]))
     for (i in 2:T_latest){
        latest_wave[i] ~ dnorm(spline[i],1/(spline[i]))
         spline_forecast[i] ~ dnegbin(p[i],r)
          p[i] <- r/(r+spline[i])
        spline[i] <- exp(inprod(basis_2[i,], beta_w2[1:N_knots2]) - phi1*sum(spline[1:(i-1)]))
     }

     for (i in (T_latest+1):(T_latest+30)){
           spline_forecast[i] ~ dnegbin(p[i],r)
          p[i] <- r/(r+spline[i])
          spline[i] <- exp(inprod(basis_oos[i-T_latest,], beta_w2[1:N_knots2]) - phi1*sum(spline[1:(i-1)]))
    }




     spline_w1[1] <- exp(inprod(basis_[1,], beta_w1))
     for (i in 2:T_previous){
        previous_wave[i] ~ dnorm(spline_w1[i],1/(spline_w1[i]))
        spline_w1[i] <- exp(inprod(basis_[i,], beta_w1) - phi2*sum(spline_w1[1:(i-1)]))
    }



    for (i in 1:N_knots2) {
        beta_w2[i] ~ dnorm(beta[i], 10)
      }


    for (i in 1:N_knots) {
      beta_w1[i] ~ dnorm(beta[i], 10000)
    }

    beta[1] ~ dnorm(0, .0001)
    for (i in 2:max(N_knots,N_knots2)) {
      beta[i] ~ dnorm(beta[i-1], .0001)
    }
    sigma_b ~ dt(0, 10^-2, 1)T(0,)
    phi1 ~ dbeta(1,100)
        phi2 ~ dbeta(1,100)

      r <- 10

  }

  "


  ## wave to match

  wave_to_match <- cases[cases$geo_value ==geo_val & cases$time_value %in% seq(as.Date("2020-08-01"),as.Date(forecast_date),by="day"),]$value
  wave_to_match <- pmax(wave_to_match,0)
  if (any(wave_to_match > 30)){
    wave_to_match[round(wave_to_match)==0]<- NA
  }

  diff_wave_to_match <- diff(lowess(wave_to_match,f=.075)$y)
  diff_wave_to_match <- diff_wave_to_match[complete.cases(diff_wave_to_match)]

  zeros <- c()
  for (i in 2:length(diff_wave_to_match)){
    if (diff_wave_to_match[i-1] < 0 && diff_wave_to_match[i] >0 ){
      zeros <- c(zeros,1)
    } else{
      zeros <- c(zeros,0)
    }
  }




  plot(wave_to_match/max(wave_to_match,na.rm=T),type='l',ylim=c(-1,1))

  lines(diff_wave_to_match/max(diff_wave_to_match),type='l',col='red')
  abline(h=0)

  index_of_zeros <- c()
  for (zero_idx in 1:length(zeros)){
    if (zeros[zero_idx]==1){
      abline(v=zero_idx)
      index_of_zeros <- c(index_of_zeros,zero_idx)
    }
  }

  latest_wave_start <- seq(as.Date("2020-08-01"),as.Date(forecast_date),by="day")[tail(index_of_zeros,1)]

  previous_wave_start <- seq(as.Date("2020-08-01"),as.Date(forecast_date),by="day")[tail(index_of_zeros,2)[1]]
  previous_wave_end <- seq(as.Date("2020-08-01"),as.Date(forecast_date),by="day")[tail(index_of_zeros,2)[2]]



  if (is.na(previous_wave_end)){
    previous_wave_end <- previous_wave_start
    previous_wave_start <- as.Date("2020-08-01")
  }

  latest_wave <- cases[cases$geo_value ==geo_val & cases$time_value > latest_wave_start,]$value

  latest_wave[latest_wave == 0 ] <- NA

  previous_wave <- cases[cases$geo_value ==geo_val & cases$time_value %in% seq(as.Date(previous_wave_start),as.Date(previous_wave_end),by="day"),]$value

  basis_ <- ns(1:length(previous_wave),knots=seq(1,length(previous_wave)-14,by=14))

  basis_2 <- ns(1:length(latest_wave),knots=seq(1,length(latest_wave)-14,by=14))

  basis_oos <- splines:::predict.ns(basis_2,newx=seq(length(latest_wave),length(latest_wave)+30))


  model_data <- list( T_previous = length(previous_wave),
                      previous_wave =round(previous_wave)+1,
                      latest_wave =  latest_wave,
                      basis_2 = basis_2,
                      T_latest = length(latest_wave),
                      basis_ = basis_,
                      N_knots2 = dim(basis_2)[2],
                      basis_oos=basis_oos,N_knots=dim(basis_)[2])

  # Choose the parameters to watch
  model_parameters <- c("spline","spline_forecast", "smoother")

  # Run the model
  model_run <- jags(
    data = model_data,
    parameters.to.save = model_parameters,
    model.file = textConnection(model_code),
    n.chains = 4, # Number of different starting positions
    n.iter = 1000, # Number of iterations
    n.burnin = 200, # Number of iterations to remove at start
    n.thin = 2
  ) # Amo


  plot(apply(model_run$BUGSoutput$sims.list$spline,2,mean),type='l')
  lines(latest_wave,col='red')
  fit_res[[geo_val]] <- model_run$BUGSoutput$sims.list$spline_forecast
}
