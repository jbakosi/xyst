-- vim: filetype=lua:

print "Initial uniform mesh refinement"

nstep = 10
ttyi = 1
cfl = 0.8

part = "rcb"

problem = {
  name = "slot_cyl"
}

mat = { spec_heat_ratio = 5/3 }

bc_dir = {
 { 1, 1, 1, 1, 1, 1, 1 },
 { 2, 1, 1, 1, 1, 1, 1 },
 { 3, 1, 1, 1, 1, 1, 1 },
 { 4, 1, 1, 1, 1, 1, 0 },
 { 5, 1, 1, 1, 1, 1, 1 },
 { 6, 1, 1, 1, 1, 1, 0 }
}

href = {
 t0 = true,
 init = { "uniform" }
}

fieldout = {
  iter = 2
}

diag = {
  iter = 1,
  format = "scientific"
}
