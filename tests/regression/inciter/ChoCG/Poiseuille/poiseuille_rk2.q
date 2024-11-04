-- vim: filetype=lua:

print "Poiseuille flow"

-- mesh: poiseuille?tetz.exo

nstep = 20
ttyi = 1
cfl = 0.5

solver = "chocg"
flux = "damp2"
rk = 2

fct = false

mat = { dyn_viscosity = 0.01 }

pressure = {
  iter = 500,
  tol = 1.0e-3,
  --verbose = 1,
  pc = "jacobi",
  bc_dir = { { 1, 2 },
             { 2, 2 } },
  bc_dirval = { { 1, 2.4 },
                { 2, 0.0 } }
}

ic = {
  velocity = { 0.0, 0.0, 0.0 }
}

bc_noslip = {
  sideset = { 3, 4 }
}

bc_dir = {
  { 1, 0, 1, 1 },
  { 5, 0, 1, 1 }
}

fieldout = {
  iter = 5
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
