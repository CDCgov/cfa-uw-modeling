# example simulation
# does not use any project-specific modules yet

library(EpiModel)
# load networks
load("mgen/input/fit_temp.rda")

# specify params & simulation controls
params <- param.net(inf.prob = 1, rec.rate = 1 / 60, act.rate = 1)
inits <- init.net(i.num = 10000 * 0.01)
controls <- control.net(type = "SIS", nsims = 1, nsteps = 100)

# simulate
sim <- netsim(est, params, inits, controls)

# plot prevalence and daily new infecteds
plot(sim, y = c("i.num", "si.flow"))
