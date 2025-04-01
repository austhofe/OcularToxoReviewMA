ma.prop.mix <- function(e, n, data, mass0 = TRUE,   
                        n.adapt = 5000,  # Ensures stability in adaptation
                        n.chains = 3, n.burnin = 5000, n.iter = 20000, n.thin = 2,   
                        seed = 1234) {
  
  # If data frame is provided, extract variables
  if (!missing(data)) {
    e <- eval(substitute(e), data, parent.frame())
    n <- eval(substitute(n), data, parent.frame())
  }
  if (length(e) != length(n)) stop("The lengths of event counts and sample sizes do not match.")
  
  # Remove missing values
  idx.na <- which(is.na(e) | is.na(n))
  if (length(idx.na) > 0) {
    e <- e[-idx.na]
    n <- n[-idx.na]
    if(length(idx.na) == 1) print("1 study with NA data entry was removed.")
    if(length(idx.na) > 1) print(paste(length(idx.na), "studies with NA data entry were removed."))
  }
  
  
  # Define the JAGS model as a string
  model.string <- if (mass0) {
    "model{
      for (i in 1:k) {
        e[i] ~ dbinom(theta[i], n[i])
        theta[i] <- (1 - zeta[i]) * beta[i]
        beta[i] <- 1 / (1 + exp(-eta[i]))
        eta[i] ~ dnorm(mu, prec)
        zeta[i] ~ dbern(p.mass)
      }
      ## Priors
      mu ~ dnorm(0, 0.0001)
      prec <- 1 / tau2
      tau2 <- tau^2
      tau ~ dunif(0, 1000)
      p.mass ~ dunif(0, 1)
      prop <- 1 / (1 + exp(-mu))
      ## Prediction
      eta.new ~ dnorm(mu, prec)
      pred.nonmass <- 1 / (1 + exp(-eta.new))
      zeta.new ~ dbern(p.mass)
      pred.overall <- (1 - zeta.new) * pred.nonmass
    }"
  } else {
    "model{
      for (i in 1:k) {
        e[i] ~ dbinom(theta[i], n[i]) 
        theta[i] <- (1 - zeta[i]) * beta[i] + zeta[i]
        beta[i] <- 1 / (1 + exp(-eta[i]))
        eta[i] ~ dnorm(mu, prec)
        zeta[i] ~ dbern(p.mass)
      }
      ## Priors
      mu ~ dnorm(0, 0.0001)
      prec <- 1 / tau2
      tau2 <- tau^2
      tau ~ dunif(0, 1000)
      p.mass ~ dunif(0, 1)
      prop <- 1 / (1 + exp(-mu))
      ## Prediction
      eta.new ~ dnorm(mu, prec)
      pred.nonmass <- 1 / (1 + exp(-eta.new))
      zeta.new ~ dbern(p.mass)
      pred.overall <- (1 - zeta.new) * pred.nonmass + zeta.new
    }"
  }
  model.spec <- textConnection(model.string)
  
  # 🔥 Ensure JAGS gets a consistent RNG seed 🔥
  set.seed(seed)  # Set global seed in R
  #rng.seeds <- seq(seed, by = 1, length.out = n.chains)  # Unique but deterministic seeds per chain
  
  # Initialize JAGS chains with fixed seeds
  #inits <- vector("list", n.chains)
  #for (i in 1:n.chains) {
    #inits[[i]] <- list(mu = 0,  # Fixed value instead of rnorm()
                       #.RNG.name = "base::Wichmann-Hill", 
                       #.RNG.seed = rng.seeds[i])  # Ensures deterministic behavior
  #}
  
  inits<-list(list(mu = 0,  # Fixed value instead of rnorm()
       .RNG.name = "base::Wichmann-Hill", 
       .RNG.seed = 1234),
  list(mu = -1,  # Fixed value instead of rnorm()
       .RNG.name = "base::Wichmann-Hill", 
       .RNG.seed = 3241),
  list(mu = 1,  # Fixed value instead of rnorm()
       .RNG.name = "base::Wichmann-Hill", 
       .RNG.seed = 1424))
  
  # Prepare data
  jags.dat <- list(e = e, n = n, k = length(e))
  params <- c("p.mass", "prop", "tau", "pred.nonmass", "pred.overall")
  
  # 🔥 Force JAGS adaptation to be stable 🔥
  model <- jags.model(model.spec, data = jags.dat, inits = inits, 
                      n.chains = n.chains, n.adapt = n.adapt)
  
  update(model, n.burnin) ## Burn-in phase
  
  # Sample from posterior
  save <- jags.samples(model, params, n.iter = n.iter, thin = n.thin)  
  
  # Extract summary statistics
  p.mass.mean <- mean(save$p.mass)
  prop.median <- quantile(save$prop, probs = 0.5)
  prop.CrI <- quantile(save$prop, probs = c(0.025, 0.975))
  PI.nonmass <- quantile(save$pred.nonmass, probs = c(0.025, 0.975))
  PI.overall <- quantile(save$pred.overall, probs = c(0.025, 0.975))
  
  smry <- list(
    p.mass.mean = p.mass.mean, 
    prop.median = prop.median, 
    prop.CrI = prop.CrI, 
    PI.nonmass = PI.nonmass, 
    PI.overall = PI.overall
  )
  
  out <- list(save = save, smry = smry)
  return(out)
}
