-- vim: filetype=lua:

print "Sod shocktube"

term = 0.2
ttyi = 1
cfl = 0.1

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
  iter = 20,
  time = 4.0e-3,
  range = {
    { 0.1,  0.15, 1.0e-3 },
    { 0.05, 0.12, 1.0e-5 }
  },
  points = {
    { 0.1, 0.05, 0.025 },
    { 0.9, 0.05, 0.025 },
    { 0.5, 0.05, 0.025 }
  },
  precision = 6
}

diag = {
  iter = 1,
  format = "scientific"
}
