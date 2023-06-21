## generate test data

## phenomeno function
sp_numbers <- function(Time){
  num = 10^6.4 * exp(-0.024 * Time)
  return(num)
}

# ODE formulation for the neutral model
neutraldd <- function(t, y, parms) {
  u <- y[1]  
  v <- y[1]
  
  eps=0.33
  Beta = 1/3.5
  
  FO <- parms[1]
  alpha <- parms[2]
  delta <- parms[3]
  mu <- parms[4]
  beta <- parms[5]
  lambda <- parms[6]
  
  dudt <- psi * (eps) * sp_numbers(t) + rho * (2*v + u) - (delta + Beta) * u
  dvdt <- psi * (1-eps) * sp_numbers(t) + Beta * u - (delta + rho) v
  
  return(list(c(dudt, dvdt)))
}