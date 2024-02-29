-- vim: filetype=lua:

print "Scalar transport: slotted cylinder, cone, hump"

nstep = 20
-- term 1.57
ttyi = 10
cfl = 0.9

solver = "riecg"
flux = "hllc"

part = "rcb"

problem = {
  name = "slot_cyl"
}

mat = { spec_heat_ratio = 5/3 }

bc_dir = {
 { 1, 1, 1, 1, 1, 1, 1 },
 { 2, 1, 1, 1, 1, 1, 0 }
}

fieldout = {
  iter = 100
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
