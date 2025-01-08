# example simulation

# libs and functions
library(EpiModel)
devtools::load_all()

# load networks
load("mgen/input/fit_temp.rda")

# load config file
p <- config::get(file = "mgen/input/model_params.yaml", use_parent = FALSE)

# specify params & simulation controls
params <- param.net(
  recov.state        = p$recov.state,
  act.rate.main      = p$act.rate.main,
  act.rate.casual    = p$act.rate.casual,
  act.rate.inst      = p$act.rate.inst,
  inf.prob           = p$inf.prob,
  inf.prob.g2        = p$inf.prob.g2,
  sympt.prob         = p$sympt.prob,
  sympt.prob.g2      = p$sympt.prob.g2,
  infectious.rate    = p$infectious.rate,
  infectious.rate.g2 = p$infectious.rate.g2,
  rec.rate           = p$rec.rate,
  rec.rate.g2        = p$rec.rate.g2
)

inits <- init.net(
  i.num    = p$i.num,
  i.num.g2 = p$i.num.g2
)

controls <- control.net(
  nsims = p$nsims,
  nsteps = p$nsteps,
  ncores = 1,
  infection.FUN = transmit_mgen,
  progress.FUN = progression_mgen,
  save.other = p$save.other,
  verbose = p$verbose
)

# simulate
sim <- netsim(est, params, inits, controls)

# plot prevalence and daily new infecteds
plot(sim, y = c("i.num", "se.flow", "se.flow.g2"))
plot(sim, y = c("ei.flow", "ei.flow.g2"))
plot(sim, y = c("is.flow", "is.flow.g2"))
