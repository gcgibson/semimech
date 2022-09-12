library(rjags)
library(R2jags)
library(splines)

run_semi_mech <- function(season_past_matrix_hosp,season_current_hosp,length_of_current_wave){

  model_code <- "
  model
  {

      for (i in (length_of_current_wave+1):(length_of_current_wave+30)){
           season_past_hosp[num_past_seasons,i] ~ dnegbin(p_past[i,num_past_seasons],r)
            p_past[i,num_past_seasons] <- r/(r+season_previous_spline[i,num_past_seasons])
            f_t_past[i,num_past_seasons] <- inprod(basis_oos[i-length_of_current_wave,], beta_w1[1:N_knots_current,num_past_seasons])
            exp_f_t_past[i,num_past_seasons] <- exp(f_t_past[i,num_past_seasons])
          season_previous_spline[i,num_past_seasons] <- exp(f_t_past[i,num_past_seasons] - phi[num_past_seasons]*sum(exp_f_t_past[i,num_past_seasons]))
       }

      #####
       phi[num_past_seasons]  ~ dbeta(1,10000)
       ## Spline model for previous wave
       season_previous_spline[1,num_past_seasons] <- inprod(season_current_basis[1,], beta_w1[1:N_knots_current,num_past_seasons])
       for (i in 2:length_of_current_wave){
           season_past_hosp[num_past_seasons,i] ~ dnegbin(p_past[i,num_past_seasons],r)
            p_past[i,num_past_seasons] <- r/(r+season_previous_spline[i,num_past_seasons])
            f_t_past[i,num_past_seasons] <- inprod(season_current_basis[i,], beta_w1[1:N_knots_current,num_past_seasons])
            exp_f_t_past[i,num_past_seasons] <- exp(f_t_past[i,num_past_seasons])
          season_previous_spline[i,num_past_seasons] <- exp(f_t_past[i,num_past_seasons] - phi[num_past_seasons]*sum(exp_f_t_past[i,num_past_seasons]))
       }
     for (j in 1:N_knots_current) {
        beta_w1[j,num_past_seasons] ~ dnorm(beta[j], .001)
     }




    #### Past seasons
    for (j in 1:(num_past_seasons-1)){
        phi[j]  ~ dbeta(1,100)
       ## Spline model for previous wave
       season_previous_spline[1,j] <- inprod(season_past_basis[1,], beta_w1[,j])
       for (i in 2:T_previous){
           season_past_hosp[j,i] ~ dnegbin(p_past[i,j],r)
            p_past[i,j] <- r/(r+season_previous_spline[i,j])
            f_t_past[i,j] <- inprod(season_past_basis[i,], beta_w1[,j])
            exp_f_t_past[i,j] <- exp(f_t_past[i,j])
          season_previous_spline[i,j] <- exp(f_t_past[i,j] - phi[j]*sum(exp_f_t_past[i,j]))
       }
     for (i in 1:N_knots_season_past_basis) {
        beta_w1[i,j] ~ dnorm(beta[i], 1)
     }
      b[j] ~ dnorm(0,.001)
      intercept_b[j] ~ dnorm(0,.001)
    }




   a ~ dnorm(0,.001)


   intercept ~ dnorm(0,.001)


    beta[1] ~ dnorm(0, .0001)
    for (i in 2:max(N_knots_season_past_basis)) {
      beta[i] ~ dnorm(beta[i-1], .0001)
    }
    sigma_b ~ dt(0, 10^-2, 1)T(0,)
    phi2 ~ dbeta(1,100)

    r <-  10

  }

  "



  ## create past season basis
  season_past_basis <- ns(1:dim(season_past_matrix_hosp)[2],knots=seq(1,dim(season_past_matrix_hosp)[2]-14,by=14))


  ## create current season basis
  season_current_basis <- ns(1:(length_of_current_wave),knots=seq(1,(length_of_current_wave)-14,by=14))

  ## get out of sample basis for 30 days ahead
  season_current_basis_oos <- splines:::predict.ns(season_current_basis,newx=seq((length_of_current_wave),(length_of_current_wave)+30))


  model_data <- list( T_previous = dim(season_past_matrix_hosp)[2],
                      num_past_seasons =  dim(season_past_matrix_hosp)[1],
                      season_past_hosp =season_past_matrix_hosp,
                     # season_current_hosp =  season_current_hosp,
                      season_current_basis = season_current_basis,
                     # T_latest = length(season_current_hosp),
                      season_past_basis = season_past_basis,
                     length_of_current_wave=length_of_current_wave,
                      N_knots_season_past_basis = dim(season_past_basis)[2],
                     N_knots_current = dim(season_current_basis)[2],
                     basis_oos=season_current_basis_oos)

  # Choose the parameters to watch
  model_parameters <- c("season_previous_spline","spline_forecast_current","phi","phi2")

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

  dims <-  dim(model_run$BUGSoutput$sims.list$season_previous_spline)

  current_season <- model_run$BUGSoutput$sims.list$season_previous_spline[,1:(length_of_current_wave+30),dims[3]]
  dim(current_season)
  ret_list <- list()
  ret_list[[1]] <- current_season#model_run$BUGSoutput$sims.list$season_current_spline#[,tail(1:dim(model_run$BUGSoutput$sims.list$spline_forecast_current)[2],30)]
  ret_list[[2]] <- model_run$BUGSoutput$sims.list$phi
  ret_list[[3]] <- model_run$BUGSoutput$sims.list$phi2

  return( ret_list)
}
