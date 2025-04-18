-- vim: filetype=lua:

print "Euler equations computing the Rayleigh-Taylor"

nstep = 50
ttyi = 1
cfl = 0.5

part = "rcb"

problem = {
  name = "rayleigh_taylor",
  alpha = 1.0,
    beta = { 1.0, 1.0, 1.0 },
    p0 = 1.0,
    r0 = 1.0,
    kappa = 1.0 -- kappa = 0 -> stationary
}

mat = { spec_heat_ratio = 5/3 }

bc_dir = {
  { 1, 1, 1, 1, 1, 1 },
  { 2, 1, 1, 1, 1, 1 },
  { 3, 1, 1, 1, 1, 1 },
  { 4, 1, 1, 1, 1, 1 },
  { 5, 1, 1, 1, 1, 1 },
  { 6, 1, 1, 1, 1, 1 }
}

fieldout = {
  range = { 0.1, 0.15, 1.0e-3 }
}

diag = {
  iter = 1,
  format = "scientific"
}
