-- vim: filetype=lua:

print "Sedov blast wave"

nstep = 20
term = 1.0
ttyi = 1
cfl = 0.5

solver = "zalcg"
fctsys = { 1, 2, 3, 4, 5 }

part = "phg"
zoltan_params = {
   "DEBUG_LEVEL", "2"
 , "CHECK_HYPERGRAPH", "1"
}

problem = {
  name = "sedov",
     p0 = 4.86e+3      -- sedov_coarse.exo,   V(orig) = 8.50791e-06
  -- p0 3.63e+4        -- sedov00.exo,        V(orig) = 1.13809e-06
  -- p0 2.32e+5        -- sedov01.exo,        V(orig) = 1.78336e-07
  -- p0 1.85e+6        -- sedov02.exo,        V(orig) = 2.22920e-08
  -- p0 14773333.33333 -- sedov02.exo+t0ref:u,V(orig) = 2.78650e-09
}

mat = { spec_heat_ratio = 5/3 }

bc_sym = {
  sideset = { 1, 2, 3 }
}

diag = {
  iter = 1,
  format = "scientific"
}

fieldout = {
  iter = 5
}
