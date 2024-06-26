-- vim: filetype=lua:

print "Poisson solve with a constant rhs and Dirichlet BCs"

nstep = 1
ttyi = 1
dt = 0.5

solver = "chocg"

part = "phg"

problem = {
  name = "poisson_const"
}

pressure = {
  iter = 100,
  tol = 1.0e-6,
  verbose = 2,
  bc_dir = {
    { 1, 1 },
    { 2, 1 },
    { 3, 1 },
    { 4, 1 },
    { 5, 1 },
    { 6, 1 }
  }
}

fieldout = {
  iter = 1000
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
