-- vim: filetype=lua:

print "Pipe computing mass flow rate"

-- Percent error between inflow and outflow mass flow rate:
-- plot "out.int" u 2:(($4+$5)/$5*100)

nstep = 20
-- term = 0.01
ttyi = 1
cfl = 0.5

part = "rcb"

mat = { spec_heat_ratio = 1.4 }

ic = {
  density = 1.0,
  velocity = { 0, 0, 0 },
  -- ic pressure = (p_in + p_out) / 2
  pressure = 99950
}

bc_sym = {
  sideset = { 2, 4, 5, 6 }
}

bc_pre = {
  inlet = {
    sideset = { 1 },
    density = 1.0,
    pressure = 100000
  },
  outlet = {
    sideset = { 3 },
    density = 1.0,
    pressure = 99900
  }
}

fieldout = {
  iter = 10
}

integout = {
  iter = 1,
  -- time = 4.0e-3
  -- range = { 0.05, 0.12, 1.0e-5 },
  sideset = { 1, 3 },
  format = "scientific",
  -- precision = 12
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
