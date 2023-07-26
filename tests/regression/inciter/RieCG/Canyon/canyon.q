-- vim: filetype=lua:

print "Dispersion in street canyon from a point source"

nstep = 50       -- Max number of time steps
-- term = 1000.0 -- Max time
cfl = 0.5
ttyi = 10        -- TTY output interval

part = "rcb"

problem = {
  name = "point_src",
  src = {
    location = { 3.0, 0.01, 0.0 },
    radius = 0.2,
    release_time = 0.0
  }
}

ic = {
  density = 1.225,
  pressure = 1.0e+5,
  velocity = { 0.0, 0.0, 0.0 }
}

mat = { spec_heat_ratio = 1.4 }

bc_sym = {
  sideset = { 1, 2, 3, 4, 5 }
}

bc_pre = {
  inlet = {
    density = 1.225,
    pressure = 1.0e+5,
    sideset = { 6 }
  },
  outlet = {
    density = 1.225,
    pressure = 0.9e+5,
    sideset = { 7 }
  }
}

fieldout = {
  iter = 10
}

histout = {
  -- iter = 1,
  -- time = 4.0e-3,
  range = { 15.0, 20.0, 1.0e-5 },
  points = {
    { 0.5,   2.001, 0.0 },
    { 1.0,   2.001, 0.0 },
    { 1.5,   2.001, 0.0 },
    { 2.001, 2.0,   0.0 },
    { 2.001, 1.93,  0.0 },
    { 2.001, 1.5,   0.0 },
    { 2.001, 1.33,  0.0 },
    { 2.001, 1.0,   0.0 },
    { 2.001, 0.67,  0.0 },
    { 2.001, 0.5,   0.0 },
    { 2.001, 0.33,  0.0 },
    { 2.001, 0.17,  0.0 },
    { 3.999, 0.17,  0.0 },
    { 3.999, 0.33,  0.0 },
    { 3.999, 0.5,   0.0 },
    { 3.999, 0.67,  0.0 },
    { 3.999, 1.0,   0.0 },
    { 3.999, 1.33,  0.0 },
    { 3.999, 1.5,   0.0 },
    { 3.999, 1.93,  0.0 },
    { 3.999, 2.0,   0.0 },
    { 4.5,   2.001, 0.0 },
    { 5.0,   2.001, 0.0 },
    { 5.5,   2.001, 0.0 }
  }
}

diag = {
  iter = 10,
  format = "scientific",
  precision = 6
}
