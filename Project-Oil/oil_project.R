# Shane Sarnac
# APPM 4660 
# Oil Project 
# April 18, 2017

# Install necessary packages
install.packages("pracma")
library("pracma")

# Read in file
oil_file = "OilProductionData.csv"
oil_data = read.csv(oil_file, header = TRUE)

# Define Constants
epsilon = 1.0*(10^-6)
oil_scale_factor = 1.0*(10^4)
Q_inf = 2*(10^6)/oil_scale_factor

# Remove NA values from world oil data set
world_oil_data = oil_data[complete.cases(oil_data),]

world_oil_max = max(world_oil_data$World.Oil..10.6.barrels.yr.)
year_offset = world_oil_data$Year[which(world_oil_data$World.Oil..10.6.barrels.yr. == world_oil_max)]
world_oil_data$World_Oil_Scaled = world_oil_data$World.Oil..10.6.barrels.yr./oil_scale_factor
world_oil_data$Year_Scaled = world_oil_data$Year - year_offset

# Define functions
f = function(mu_0, sigma_0) {
  f_values = Q_inf*exp(-0.5*((world_oil_data$Year_Scaled - mu_0)/sigma_0)^2)/(sigma_0*sqrt(2*pi))
}

f.mu = function(mu_0, sigma_0) {
  f_mu_values = Q_inf*(world_oil_data$Year_Scaled - mu_0)
  f_mu_values = f_mu_values*exp(-0.5*((world_oil_data$Year_Scaled - mu_0)/sigma_0)^2)
  f_mu_values = f_mu_values/((sigma_0^3)*(sqrt(2*pi)))
}

f.sigma = function(mu_0, sigma_0) {
  f_sigma_values1 = Q_inf*(world_oil_data$Year_Scaled - mu_0)^2
  f_sigma_values1 = f_sigma_values1*exp(-0.5*((world_oil_data$Year_Scaled - mu_0)/sigma_0)^2)
  f_sigma_values1 = f_sigma_values1/((sigma_0^4)*sqrt(2*pi))
  
  f_sigma_values2 = Q_inf*exp(-0.5*((world_oil_data$Year_Scaled - mu_0)/sigma_0)^2)
  f_sigma_values2 = f_sigma_values2/((sigma_0^2)*sqrt(2*pi))
  
  f_sigma_values = f_sigma_values1 - f_sigma_values2
}

determine_mu_and_sigma = function(mu_0, sigma_0) {
  mu_guess = mu_0
  sigma_guess = sigma_0
  
  for (i in 1:500) {
    f_values = f(mu_guess, sigma_guess)
    f_mu = f.mu(mu_guess, sigma_guess)
    f_sigma = f.sigma(mu_guess, sigma_guess)
    
    f_mu_squared = f_mu^2
    f_sigma_squared = f_sigma^2
    f_sigma_mu = f_mu*f_sigma
    difference_mu = (world_oil_data$World_Oil_Scaled - f_values)*f_mu
    difference_sigma = (world_oil_data$World_Oil_Scaled - f_values)*f_sigma
    
    sum_f_mu_squared = sum(f_mu_squared)
    sum_f_sigma_squared = sum(f_sigma_squared)
    sum_f_sigma_mu = sum(f_sigma_mu)
    sum_difference_mu = sum(difference_mu)
    sum_difference_sigma = sum(difference_sigma)
    
    solution_list = c(
      sum_f_mu_squared, sum_f_sigma_mu, sum_difference_mu, 
      sum_f_sigma_mu, sum_f_sigma_squared, sum_difference_sigma
    )
    
    solution_matrix = matrix(solution_list, nrow = 2, ncol = 3, byrow = TRUE)
    rref_matrix = rref(solution_matrix)
    
    mu_new = rref_matrix[1,3] + mu_guess
    sigma_new = rref_matrix[2,3] + sigma_guess
    
    if ((abs(mu_new - mu_guess) < epsilon) && (abs(sigma_new - sigma_guess) < epsilon)) {
      mu_guess = mu_new
      sigma_guess = sigma_new
      break
    }
    mu_guess = mu_new
    sigma_guess = sigma_new

    print(i)
  }
  c(mu_new, sigma_new)
}


mu_guess = 0
sigma_guess = 50
solution_data = data.frame(Year_Scaled = -200:200)

# Part 1: Q_inf = 2 trillion
Q_inf = 2.0 * (10^6) / oil_scale_factor
part_1 = determine_mu_and_sigma(mu_guess, sigma_guess)
solution_data$Part_1 = Q_inf*exp(-0.5*((solution_data$Year_Scaled - part_1[1])/part_1[2])^2)
solution_data$Part_1 = solution_data$Part_1 / (part_1[2] * sqrt(2*pi))

# Part 2: Q_inf = 3 trillion
Q_inf = 3.0 * (10^6) / oil_scale_factor
part_2 = determine_mu_and_sigma(mu_guess, sigma_guess)
solution_data$Part_2 = Q_inf*exp(-0.5*((solution_data$Year_Scaled - part_2[1])/part_2[2])^2)
solution_data$Part_2 = solution_data$Part_2 / (part_2[2] * sqrt(2*pi))

# Part 3: Q_inf = 4 trillion
Q_inf = 4.0 * (10^6) / oil_scale_factor
part_3 = determine_mu_and_sigma(mu_guess, sigma_guess)
solution_data$Part_3 = Q_inf*exp(-0.5*((solution_data$Year_Scaled - part_3[1])/part_3[2])^2)
solution_data$Part_3 = solution_data$Part_3 / (part_3[2] * sqrt(2*pi))


#### Plot Results
par(mfrow = c(1,1))
# Part 1
plot(world_oil_data$Year_Scaled, world_oil_data$World_Oil_Scaled, 
     xlim = c(-200,200), ylim = c(0, 4),
     xlab = "Year (2006 is Year 0)", ylab = "Oil Produced", 
     main = "Theoretical Oil Production Curve (Q_inf = 2 Trillion)"
)
lines(solution_data$Year_Scaled, solution_data$Part_1, type = "l")
abline(h = 0.05*max(solution_data$Part_1))

# Part 2
plot(world_oil_data$Year_Scaled, world_oil_data$World_Oil_Scaled, 
     xlim = c(-200,200), ylim = c(0, 4),
     xlab = "Year (2006 is Year 0)", ylab = "Oil Produced", 
     main = "Theoretical Oil Production Curve (Q_inf = 3 Trillion)"
)
lines(solution_data$Year_Scaled, solution_data$Part_2, type = "l")
abline(h = 0.05*max(solution_data$Part_2))

# Part 3
plot(world_oil_data$Year_Scaled, world_oil_data$World_Oil_Scaled, 
     xlim = c(-200,200), ylim = c(0, 4),
     xlab = "Year (2006 is Year 0)", ylab = "Oil Produced", 
     main = "Theoretical Oil Production Curve (Q_inf = 4 Trillion)"
)
lines(solution_data$Year_Scaled, solution_data$Part_3, type = "l")
abline(h = 0.05*max(solution_data$Part_3))

##### Determine year where oil production becomes less than 5% of max
year_five_percent = function(five_percent, dataset) {
  answer = 0
  for (i in 200:length(solution_data$Year_Scaled)) {
    #print(five_percent)
    #print(dataset[i])
    if (dataset[i] <= five_percent) {
      answer = i
      break
    }
  }
  answer
}

part_1_max = max(solution_data$Part_1)
part_2_max = max(solution_data$Part_2)
part_3_max = max(solution_data$Part_3)

part_1_5_percent = 0.05*part_1_max
part_2_5_percent = 0.05*part_2_max
part_3_5_percent = 0.05*part_3_max

part_1_year_under_5_percent = 
  solution_data$Year_Scaled[year_five_percent(part_1_5_percent, solution_data$Part_1)]
part_2_year_under_5_percent = 
  solution_data$Year_Scaled[year_five_percent(part_2_5_percent, solution_data$Part_2)]
part_3_year_under_5_percent = 
  solution_data$Year_Scaled[year_five_percent(part_3_5_percent, solution_data$Part_3)]
