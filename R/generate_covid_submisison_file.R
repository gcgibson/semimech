
library(tidyverse)
library(cowplot)
library(lubridate)
library(mgcv)
library(splines)

library(MASS)
library(dplyr)
library(rstanarm)
#options(mc.cores = parallel::detectCores())
last_date <- as.Date("2022-08-28")
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
    mutate(covid = ifelse(is.na(covid), 0, covid)) %>%
    mutate(time = as.numeric(date - min(date)))

  smoothed_deriv <- c(0,diff(lowess(mydf$covid/max(mydf$covid,na.rm=T),f=.1)$y))


  wave_starts <- c()

  for (t in 2:length(smoothed_deriv)){
    if (!is.na(smoothed_deriv[t-1]) & !is.na(smoothed_deriv[t])){
      if (smoothed_deriv[t-1] < 0 & smoothed_deriv[t] >0){
        wave_starts <-c(wave_starts,t)
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



  forecast <- run_semi_mech(season_past_matrix_hosp = wave_matrix,season_current_hosp =current_wave )

  ### also not really neccessary


  dflist[[k]] <- forecast

}

upper_q <- function(x){quantile(x,probs=c(.975))}
lower_q <- function(x){quantile(x,probs=c(.025))}

dflist_as_df <- list()

for (k in 1:length(statevec)) {
  forecast_dates <- seq(as.Date(tail(hhs$date,1)),as.Date(tail(hhs$date,1))+29,by="day")
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
  geom_point(data=hhs_sub %>%  filter(date >= "2022-01-01"),aes(date, covid), alpha=0.25) +
  facet_wrap(~state,scales="free")  +
  # scale_y_log10() +
  byy1() +
  scale_x_date(date_labels = "%b '%y") +
  labs(x="", y="covid admissions") + scale_x_date(date_breaks = '1 month', date_labels = '%b') +
 # geofacet::facet_geo(~state,scales='free') +
  theme(panel.background = element_rect(fill='white', colour='white'),
        plot.background =element_rect(fill='white', colour='white') )

ggsave(filename = "covid_hosps.png",covid_hosps,device = "png",height = 10,width = 16)




#####


df_to_submit <- data.frame(matrix(NA,ncol=7,nrow=0))
colnames(df_to_submit) <- c("target","location","forecast_date","target_end_date","quantile","value","type")


for (k in 1:length(statevec)){
  # extract this states df
  state_forecast_df <- dflist[[k]]


  # define quantiles to populate
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)

  # extract quantiles
  state_forecast_df_daily <- apply(state_forecast_df,2,function(x){
    return (quantile(x,probs=quantiles))
  })

  # make dataframe for this state with correct dimensions, number of quantiles (23)
  # by number of days ahead (1-30)
  state_df <- data.frame(h=rep(1:30,23),value =c(t(state_forecast_df_daily)),
                         quantile = rep(quantiles,each=30))

  # assign correct state
  state_df$state <- statevec[k]

  # modify week ahead to correct formatting
  state_df$target <- unlist(lapply(state_df$h , function(x){
    return (paste0(x, " day ahead inc hosp"))
  }))

  # assign correct forecast date (latest Monday)

  state_df$forecast_date <- last_date + 1

  # assign corresponding next saturday modulo week ahea as the target date
  state_df$target_end_date <- Reduce(c,lapply(state_df$h , function(x){
    return (last_date + x + 1)
  }))

  # assign the correct state
  state_df$location <- statevec[k]

  # add type column
  state_df$type <- "quantile"
  # select only columns to be submitted
  state_df <- state_df %>% dplyr::select(target,location,forecast_date,target_end_date,quantile,value,type)



  df_to_submit <- rbind(df_to_submit,state_df)
}


##
library(lubridate)
hhs_sub <- hhs %>%
  dplyr::select(state, date, fips:covid_per_cap_av7) %>%
  filter(!is.na(flu) | !is.na(covid))

df_to_submit$date <- df_to_submit$target_end_date
hhs_sub$location <- hhs_sub$state
df_to_submit$covid <- df_to_submit$value
hhs_sub <- rbind(hhs_sub %>% dplyr::select(location,date,covid),df_to_submit[df_to_submit$quantile == .5,] %>% dplyr::select(location,date,covid))

forecasts <- ggplot(hhs_sub[as.Date(hhs_sub$date) > start_date,],aes(x=as.Date(date),y=covid)) + geom_line()# + facet_wrap(~location,scales='free')
forecasts <- forecasts + facet_wrap(~location,scales='free')
forecasts <- forecasts + geom_point(data=df_to_submit[df_to_submit$quantile == .5,],aes(x=date,y=covid),col='red',size=.2)
ggsave(forecasts + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) ,filename = "covid_forecasts.png",device = "png",height = 10,width = 12)



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
df_to_submit[is.na(df_to_submit$location),]$location <- "US"
write.csv(df_to_submit,file = paste0(last_date + 1,"-UT-Osiris.csv"),row.names = F)




