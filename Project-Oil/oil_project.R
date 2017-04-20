# Shane Sarnac
# APPM 4660 
# Oil Project 
# April 18, 2017

# Read in file
oil_file = "OilProductionData.csv"
oil_data = read.csv(oil_file, header = TRUE)

# Remove NA values from world oil data set
world_oil_data = oil_data[complete.cases(oil_data),]
year_offset = world_oil_data$Year[which(world_oil_data$World.Oil..10.6.barrels.yr. == world_oil_max)]

world_oil_max = max(world_oil_data$World.Oil..10.6.barrels.yr.)
world_oil_data$World_Oil_Scaled = world_oil_data$World.Oil..10.6.barrels.yr./world_oil_max
world_oil_data$Year_Scaled = world_oil_data$Year - year_offset

f = function(mu_0, sigma_0) {
  Q_inf = 2000000
  f_values = pow(sigma_0*sqrt(2*pi), -1)*Q_inf*exp(-0.5*pow((t - mu_0)/sigma_0,2))
}

plot(oil_data$Year, oil_data$World.Oil..10.6.barrels.yr.)
plot(oil_data$Year, oil_data$US.Oil..10.3.barrels.yr.)

