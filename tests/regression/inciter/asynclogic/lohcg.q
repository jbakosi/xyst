-- vim: filetype=lua:

print "Asynclogic test"

nstep = 10    -- Max number of time steps
dt = 0.001    -- Time step size
ttyi = 1      -- TTY output interval

solver = "lohcg"
flux = "damp2"

ic = {
  velocity = { 0, 0, 0 }
}

bc_sym = {
  sideset = { 1 }
}

mat = {
  dyn_viscosity = 0.01
}
