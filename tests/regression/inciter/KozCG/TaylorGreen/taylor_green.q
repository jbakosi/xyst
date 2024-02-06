-- vim: filetype=lua:

print "Euler equations computing stationary Taylor-Green"

term = 1.0
ttyi = 10
cfl = 0.8
part = "rcb"

solver = "kozcg"
fct = false

mat = { spec_heat_ratio = 5/3 }

problem = {
  name = "taylor_green"
}

bc_dir = {
  { 1, 1, 1, 1, 1, 1 },
  { 2, 1, 1, 1, 1, 1 },
  { 3, 1, 1, 1, 1, 1 },
  { 4, 1, 1, 1, 1, 1 },
  { 5, 1, 1, 1, 1, 1 },
  { 6, 1, 1, 1, 1, 1 }
}

fieldout = {
  iter = 20
}

diag = {
  iter = 2,
  format = "scientific"
}
