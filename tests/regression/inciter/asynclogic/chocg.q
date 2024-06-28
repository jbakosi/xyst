-- vim: filetype=lua:

print "Asynclogic test"

nstep = 10    -- Max number of time steps
dt = 0.001    -- Time step size
ttyi = 1      -- TTY output interval

solver = "chocg"

ic = {
  velocity = { 0, 0, 0 }
}

bc_sym = {
  sideset = { 1 }
}

mat = {
  spec_heat_ratio = 1.4
}
