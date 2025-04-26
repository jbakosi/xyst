-- vim: filetype=lua:

print "Sphere with scalar point source "

-- mesh: sphere2_*.exo

nstep = 20
ttyi = 1
cfl = 0.5

solver = "chocg"
flux = "damp2"

part = "phg"

problem = {
  name = "point_src",
  src = {
    location = { -4.95, 0.0, 0.0 },
    radius = 2.0,
    release_time = 0.0
  }
}

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
  { 2, 1, 0, 0, 0 }
}

bc_sym = {
  sideset = { 1, 4 }
}

fieldout = {
  iter = 5,
  sideset = { 1 },
  time = 0.0
}

diag = {
  iter = 5,
  format = "scientific",
  precision = 12
}
