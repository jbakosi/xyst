-- vim: filetype=lua:

print "Euler equations computing nonlinear energy growth"

term = 1.0
ttyi = 1
cfl = 0.8

part = "rcb"

problem = {
  name = "nonlinear_energy_growth",
  alpha = 0.25,
  beta = { 1.0, 0.75, 0.5 },
  r0 = 2.0,
  ce = -1.0,
  kappa = 0.8
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
  iter = 5
}

diag = {
  iter = 1,
  format = "scientific"
}
