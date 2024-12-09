-- vim: filetype=lua:

print "Viscous full sphere"

-- mesh: sphere2_*.exo

nstep = 20
ttyi = 1

cfl = 0.1

solver = "lohcg"
flux = "damp2"
stab2 = true
stab2coef = 0.01
soundspeed = 100.0
rk = 2

part = "phg"

pressure = {
  iter = 300,
  tol = 1.0e-3,
  pc = "jacobi",
  bc_dir = { { 3, 1 } }
}

Re = 100.0
mat = { dyn_viscosity = 1.0/Re }

ic = {
  velocity = { 1.0, 0.0, 0.0 }
}

bc_dir = {
  { 2, 0, 1, 1, 1 },    -- inflow
  { 3, 1, 0, 0, 0 },    -- outflow
  { 4, 0, 1, 0, 0 }     -- farfield
}

bc_noslip = {
  sideset = { 1 }
}

fieldout = {
  iter = 100,
  sideset = { 1 }
}

integout = {
  iter = 10,
  sideset = { 1 },
  integrals = { "force" }
}

diag = {
  iter = 10,
  format = "scientific",
  precision = 12
}
