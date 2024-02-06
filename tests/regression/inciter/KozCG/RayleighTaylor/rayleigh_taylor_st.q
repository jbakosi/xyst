-- vim: filetype=lua:

print "Euler equations computing the stationary Rayleigh-Taylor"

nstep = 10
ttyi = 1
cfl = 0.5
part = "rcb"

solver = "kozcg"
fct = false

mat = { spec_heat_ratio = 5/3 }

problem = {
  name = "rayleigh_taylor",
  alpha = 1.0,
    beta = { 1.0, 1.0, 1.0 },
    p0 = 1.0,
    r0 = 1.0,
    kappa = 0.0 -- kappa = 0 -> stationary
}

bc_dir = {
  { 1, 1, 1, 1, 1, 1 },
  { 2, 1, 1, 1, 1, 1 },
  { 3, 1, 1, 1, 1, 1 },
  { 4, 1, 1, 1, 1, 1 },
  { 5, 1, 1, 1, 1, 1 },
  { 6, 1, 1, 1, 1, 1 }
}

fieldout = {
  iter = 1
}

diag = {
  iter = 1,
  format = "scientific"
}
