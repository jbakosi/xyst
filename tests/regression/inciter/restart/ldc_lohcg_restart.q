-- vim: filetype=lua:

print "Lid-driven cavity"

-- mesh: mms/unitcube_94K.exo
--       mms/unitcube_750K.exo
--       mms/unitcube_6M.exo
--       mms/unitcube_48M.exo

nstep = 10
ttyi = 1

cfl = 0.9

solver = "lohcg"
flux = "damp2"
stab2 = true
stab2coef = 0.001
soundspeed = 100.0
rk = 2

fct = false

part = "rcb"

Re = 100.0
mat = { dyn_viscosity = 1.0/Re }


pressure = {
  iter = 500,
  tol = 1.0e-3,
  verbose = 1,
  pc = "jacobi",
  hydrostat = 0
}

ic = {
  velocity = { 0.0, 0.0, 0.0 }
}

bc_noslip = {
  sideset = { 1, 2, 3, 5, 6 }
}

bc_dir = {
  { 4, 0, 2, 2, 2 }
}

bc_dirval = {
  { 4, 0.0, 1.0, 0.0, 0.0 }
}

fieldout = {
  iter = 1
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
