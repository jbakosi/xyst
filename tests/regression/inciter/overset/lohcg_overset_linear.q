-- vim: filetype=lua:

print "Test one-way overset transfer of a scalar"

-- mesh 1, bakground: half_viscous_sphere_bg_11k.exo
-- mesh 2, overset: half_viscous_sphere_overset.exo

nstep = 1

solver = "lohcg"
flux = "damp2"

part_1 = "rcb"
part_2 = "phg"

overset = {
  intergrid_2 = {
    sideset = { 1 },
    layers = { 2, 4, 2 }
  }
}

problem = {
  name = "overset_linear"
}
