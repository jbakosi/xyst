-- vim: filetype=lua:

print "Scalar transport: slotted cylinder, cone, hump"

nstep = 20
-- term = math.pi
ttyi = 1
cfl = 0.9

solver = "chocg"
flux = "damp4"
rk = 4

part = "rcb"

freezeflow = 2.0
freezetime= 1.0e-1
problem = {
  name = "slot_cyl"
}

mat = { spec_heat_ratio = 5/3 }

pressure = {
  iter = 300,
  tol = 1.0e-2,
  --verbose = 1,
  pc = "jacobi",
  hydrostat = 0
}

bc_dir = {
 { 1,  1, 1, 1, 1 },
 { 2,  1, 1, 1, 0 }
}

fieldout = {
  iter = 100
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
