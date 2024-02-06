-- vim: filetype=lua:

print "Euler equations computing the Rayleigh-Taylor"

nstep = 50
ttyi = 10
cfl = 0.5

solver = "kozcg"
fct = false

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
  iter = 10
}

diag = {
  iter = 1,
  format = "scientific"
}
