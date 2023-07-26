-- vim: filetype=lua:

print "Sod shocktube"

nstep = 10
term = 0.2
ttyi = 1
cfl = 0.5

part = "rcb"

problem = {
  name = "sod"
}

mat = { spec_heat_ratio = 1.4 }

bc_sym = {
  sideset = { 2, 4, 5, 6 }
}

fieldout = {
  iter = 10000,
  sideset = { 2, 4, 5, 6 }
}

diag = {
  iter = 1,
  format = "scientific"
}
