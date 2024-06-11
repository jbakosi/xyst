-- vim: filetype=lua:

print "Bump"

nstep = 20
ttyi = 1
cfl = 0.7

solver = "laxcg"
flux = "hllc"

steady = true
residual = 1.0e-14
rescomp = 1

part = "rcb"

mat = {
  spec_heat_ratio = 1.4
}

-- free stream
rho = 1.225
pre = 101325.0
mach = 0.675
--mach = 0.001
c = math.sqrt(mat.spec_heat_ratio * pre / rho)
print("Free-stream Mach number = ",mach)

ic = {
  density = rho,
  pressure = pre, 
  velocity = { c*mach, 0.0, 0.0 }
}

velinf = ic.velocity

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
