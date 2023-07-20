-- vim: filetype=lua:

print "Euler equations computing vortical flow, steady state test"

term = 1.0
ttyi = 10
cfl = 0.8

part = "rcb"
reorder = true

problem = {
  name = "vortical_flow",
  alpha = 0.1,
  kappa = 1.0,
  p0 = 10.0,
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
