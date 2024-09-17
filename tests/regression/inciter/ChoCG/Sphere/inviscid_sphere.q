-- vim: filetype=lua:

print "Inviscid (potential) flow around a sphere"

-- mesh: sphere2_*.exo

nstep = 10
ttyi = 1
cfl = 0.5

solver = "chocg"
part = "phg"
fct = false

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
  { 2, 1, 0, 0 }
}

bc_sym = {
  sideset = { 1, 4 }
}

fieldout = {
  iter = 5
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
