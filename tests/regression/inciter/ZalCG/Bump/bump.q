-- vim: filetype=lua:

print "Bump"

nstep = 20
ttyi = 1
cfl = 0.7

solver = "zalcg"

stab2 = true
stab2coef = 0.05

steady = true
residual = 1.0e-9
rescomp = 1

part = "rib"
     
ic = {
  density = 1.0, 
  pressure = 1.0, 
  -- sound speed: a = sqrt(1.4*1.0/1.0) = 1.1832
  -- free stream Mach number: M = 0.675
  -- u = M * a = 0.7987
  velocity = { 0.7987, 0.0, 0.0 }
}

mat = { spec_heat_ratio = 1.4 }
         
bc_sym = {  
  sideset = { 3 }
}
         
bc_far = {
  density = ic.density,
  pressure = ic.pressure,
  velocity = ic.velocity,
  sideset = { 4 }
}

fieldout = {
  iter = 10
}

diag = {
  iter = 1,
  format = "scientific",
  precision = 12
}
