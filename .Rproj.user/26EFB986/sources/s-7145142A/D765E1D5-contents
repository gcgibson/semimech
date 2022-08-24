library(rjags)
library(R2jags)
library(splines)

run_semi_mech <- function(season_past_hosp,season_current_hosp){

  model_code <- "
  model
  {
     ## Forecast for current season
     for (i in (T_latest+1):(T_latest+30)){
           spline_forecast_current[i] ~ dnegbin(p_current[i],r)
          p_current[i] <- r/(r+season_current_spline[i])
          season_current_spline[i] <- exp(inprod(basis_oos[i-T_latest,], beta_w2[1:N_knots_season_current_basis]) - phi2*sum(season_current_spline[1:(i-1)]))
    }


     ## Spline model for current wave
     season_current_spline[1] <- exp(inprod(season_current_basis[1,], beta_w2[1:N_knots_season_current_basis]))
     for (i in 2:T_latest){
          season_current_hosp[i] ~ dnorm(season_current_spline[i],1/(season_current_spline[i]))
         spline_forecast_current[i] ~ dnegbin(p_current[i],r)
          p_current[i] <- r/(r+season_current_spline[i])
         season_current_spline[i] <- exp(inprod(season_current_basis[i,], beta_w2[1:N_knots_season_current_basis]) - phi2*sum(season_current_spline[1:(i-1)]))
     }





     ## Spline model for previous wave
     season_previous_spline[1] <- exp(inprod(season_past_basis[1,], beta_w1[1:N_knots_season_past_basis]))
     for (i in 2:T_previous){
          season_past_hosp[i] ~ dnorm(season_previous_spline[i],1/(season_previous_spline[i]))
         spline_forecast[i] ~ dnegbin(p_past[i],r)
          p_past[i] <- r/(r+season_previous_spline[i])
         season_previous_spline[i] <- exp(inprod(season_past_basis[i,], beta_w1[1:N_knots_season_past_basis]) - phi1*sum(season_previous_spline[1:(i-1)]))
     }


    for (i in 1:N_knots_season_current_basis) {
        beta_w2[i] ~ dnorm(beta[i], .01)
      }


    for (i in 1:N_knots_season_past_basis) {
      beta_w1[i] ~ dnorm(beta[i], .01)
    }

    beta[1] ~ dnorm(0, .0001)
    for (i in 2:max(N_knots_season_past_basis,N_knots_season_current_basis)) {
      beta[i] ~ dnorm(beta[i-1], .0001)
    }
    sigma_b ~ dt(0, 10^-2, 1)T(0,)
    phi1 ~ dbeta(1,100)
    phi2 ~ dbeta(1,100)

    r <- 10

  }

  "



  ## create past season basis
  season_past_basis <- ns(1:length(season_past_hosp),knots=seq(1,length(season_past_hosp)-14,by=14))


  ## create current season basis
  season_current_basis <- ns(1:length(season_current_hosp),knots=seq(1,length(season_current_hosp)-14,by=14))

  ## get out of sample basis for 30 days ahead
  season_current_basis_oos <- splines:::predict.ns(season_current_basis,newx=seq(length(season_current_hosp),length(season_current_hosp)+30))


  model_data <- list( T_previous = length(season_past_hosp),
                      season_past_hosp =season_past_hosp,
                      season_current_hosp =  season_current_hosp,
                      season_current_basis = season_current_basis,
                      T_latest = length(season_current_hosp),
                      season_past_basis = season_past_basis,
                      N_knots_season_past_basis = dim(season_past_basis)[2],
                      N_knots_season_current_basis = dim(season_current_basis)[2],
                      basis_oos=season_current_basis_oos)

  # Choose the parameters to watch
  model_parameters <- c("season_current_spline")

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

  plot(colMeans(model_run$BUGSoutput$sims.list$season_current_spline[,tail(1:dim(model_run$BUGSoutput$sims.list$season_current_spline)[2],30)]))
  return(model_run$BUGSoutput$sims.list$season_current_spline[,tail(1:dim(model_run$BUGSoutput$sims.list$season_current_spline)[2],30)])
}
