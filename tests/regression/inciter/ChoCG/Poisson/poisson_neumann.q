-- vim: filetype=lua:

print "Poisson solve with Dirichlet and Neumann BCs"

nstep = 1
ttyi = 1

solver = "chocg"

part = "rcb"

problem = {
  name = "poisson_neumann"
}

pressure = {
  iter = 100,
  tol = 1.0e-6,
  verbose = 2
}

bc_dir = {
  { 1, 1 },     -- x = 0
  { 2, 1 },     -- y = 0
  { 3, 1 }      -- z = 1
}

bc_sym = {
  sideset = {
    4,          -- y = 1
    5,          -- x = pi/4
    6 }         -- z = 0
}

fieldout = {
  iter = 1000
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
