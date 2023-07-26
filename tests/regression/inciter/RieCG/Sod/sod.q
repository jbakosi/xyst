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
  iter = 10000
}

histout = {
  iter = 1,
  points = {
    { 0.1, 0.05, 0.025 },
    { 0.9, 0.05, 0.025 }
  },
  precision = 6
}

diag = {
  iter = 1,
  format = "scientific"
}
