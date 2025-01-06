-- vim: filetype=lua:

print "Viscous full sphere"

-- mesh: sphere2_*.exo
-- diam = 2.0
-- plot "out.int" u 2:(-$4*2.0/pi) w l lw 2 t "Cd"

nstep = 20
ttyi = 1
cfl = 0.3

solver = "chocg"
flux = "damp4"
rk = 4
fct = false

part = "phg"
--zoltan_params = {
--  "PHG_COARSENING_METHOD", "AGG",
--  "PHG_COARSEPARTITION_METHOD", "GREEDY",
--  "PHG_REFINEMENT_QUALITY", "5",
--  "PHG_CUT_OBJECTIVE", "CONNECTIVITY"
--}

pressure = {
  iter = 300,
  tol = 1.0e-3,
  pc = "jacobi",
  bc_dir = { { 3, 1 } }
}

Re = 40.0
mat = { dyn_viscosity = 1.0/Re }

ic = {
  velocity = { 1.0, 0.0, 0.0 }
}

bc_dir = {
  { 2, 1, 0, 0 },    -- inflow
  { 4, 1, 1, 1 }     -- farfield
}

bc_noslip = {
  sideset = { 1 }
}

fieldout = {
  iter = 10,
  sideset = { 1 }
}

integout = {
  iter = 1,
  sideset = { 1 },
  integrals = { "force" }
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
