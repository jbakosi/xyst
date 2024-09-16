-- vim: filetype=lua:

print "Sphere"

nstep = 2
term = 0.2
ttyi = 1
cfl = 0.1

solver = "chocg"

pressure = {
  iter = 300,
  tol = 1.0e-3,
  pc = "jacobi",
  bc_dir = { { 3, 1 } }
}

ic = {
  velocity = { 1.0, 0.0, 0.0 }
}

bc_dir = {
  { 2, 1, 1, 1 }
}

bc_sym = {
  sideset = { 1 }
}

fieldout = {
  iter = 100,
  sideset = { 1, 2, 3 }
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
