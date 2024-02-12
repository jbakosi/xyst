-- vim: filetype=lua:

print "Euler equations computing vortical flow"

nstep = 5
ttyi = 1
cfl = 0.8
part = "rcb"

solver = "kozcg"

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
  iter = 1
}

diag = {
  iter = 1,
  format = "scientific"
}
