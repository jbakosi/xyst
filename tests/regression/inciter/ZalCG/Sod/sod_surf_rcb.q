-- vim: filetype=lua:

print "Sod shocktube"

nstep = 10
term = 0.2
ttyi = 1
cfl = 0.5

solver = "zalcg"
fctclip = true
fctsys = { 1, 2, 5 }

part = "rcb"

problem = {
  name = "sod"
}

mat = { spec_heat_ratio = 1.4 }

bc_sym = {
  sideset = { 2, 4, 5, 6 }
}

bc_dir = {
 { 1, 1, 1, 1, 1, 1 },
 { 3, 1, 1, 1, 1, 1 }
}

fieldout = {
  iter = 2,
  sideset = { 1, 5 }
}

diag = {
  iter = 1,
  format = "scientific"
}
