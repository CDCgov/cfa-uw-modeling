# example simulation

# libs and functions
library(EpiModel)
devtools::load_all()

# load networks
load("mgen/input/fit_temp.rda")

# load config file
p <- config::get(file = "mgen/input/model_params.yaml", use_parent = FALSE)

call_with_params <- function(f, p, extra = list()) {
  do.call(f,
  c(p[intersect(names(formals(f)), 
  names(p))],
    extra)
  )
}

params <- call_with_params(param.net, p)

inits <- call_with_params(init.net, p)

controls <- call_with_params(control.net, p,
    list(infection.FUN = transmit_mgen,
         progress.FUN = progression_mgen))

# simulate
sim <- netsim(est, params, inits, controls)

# plot prevalence and daily new infecteds
plot(sim, y = c("i.num", "se.flow", "se.flow.g2"))
plot(sim, y = c("ei.flow", "ei.flow.g2"))
plot(sim, y = c("is.flow", "is.flow.g2"))
