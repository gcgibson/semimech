# Load the library
setwd("R/")
library(semimech)

# Get the data
library(stringr)
hhs <- get_hhs()

# Define the forecast dates
dates <- c(as.Date("2022-11-20"))

# Get the samples
samples <- semimech::generate_samples("2020-08-01",dates,"states",hhs)

# Plot the samples
library(ggplot2)
library(geofacet)
median_mat <- matrix(NA,nrow=length(names(samples)),ncol=30)
upper_ci_mat <- matrix(NA,nrow=length(names(samples)),ncol=30)
lower_ci_mat <- matrix(NA,nrow=length(names(samples)),ncol=30)
row_idx <- 1
for (key in names(samples)){
  median_mat[row_idx,] <- colMeans(samples[[key]])
  upper_ci_mat[row_idx,] <- apply(samples[[key]],2,function(x){quantile(x,probs=.975)})
  lower_ci_mat[row_idx,] <- apply(samples[[key]],2,function(x){quantile(x,probs=.025)})
  row_idx <- row_idx + 1
}

ci_mat_df <- data.frame(median_ = c(median_mat),upper_ci=c(upper_ci_mat),lower_ci=c(lower_ci_mat),
                        state=rep( names(samples),30),date=rep(dates+1:30,each=length(names(samples))))

p <- ggplot(ci_mat_df,aes(x=date,y=median_)) +
  geom_line(col='blue') +
  geom_point(data=hhs[hhs$date > dates-365 & hhs$date < dates+30,], aes(x=as.Date(date),y=previous_day_admission_adult_covid_confirmed),size=.1)  +
  geom_ribbon(data=ci_mat_df,aes(x=date,ymax=upper_ci,ymin=lower_ci,y=median_),alpha=.5) +
  geofacet::facet_geo(~state,grid="us_state_grid2", scales="free") +
  theme(axis.text.x = element_text(size=6, angle = 90), axis.text.y = element_text(size=6))
  #facet_wrap(~state,scales="free")

ggsave(filename=paste0("../figs/",dates[1],".png"),plot=p,height = 10,width = 12)

# Generate submission file
semimech::generate_submission_file(dates,samples,"UT-Osiris","sub_files")
