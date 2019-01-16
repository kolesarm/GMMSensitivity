## Load data from AGS replication
d <- R.matlab::readMat("agm_data.mat")
H <- drop(d$H)                          # \partial h(hattheta)

g_init <- drop(d$g.init)                # hatg(hattheta)
h_init <- drop(d$h.init)                # h(hattheta)
perturb <- rep(c(drop(d$demand.perturb),
                 drop(d$supply.perturb)), times = c(13, 18))
sdZ <- drop(d$sd.Z)
ids <- list(iv=unname(unlist(d$names.iv)), th=unname(unlist(d$names.th)))

blp <- list(G=d$G, H=H, W=d$W, g_init=g_init, h_init=h_init, names=ids,
            ZZ=d$Om.ZZ, Sig=d$Omega, sdZ=sdZ, perturb=perturb, n=drop(d$n))
devtools::use_data(blp, internal=FALSE, overwrite=TRUE)

# Instruments, 13 demand (last 8 excluded) and 18 supply (last 12 excluded)
