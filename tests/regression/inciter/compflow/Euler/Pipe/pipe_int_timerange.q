# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Pipe computing mass flow rate"

# Percent error between inflow and outflow mass flow rate:
# plot "out.int" u 2:(($4+$5)/$5*100)

inciter

  nstep 20
  #term 0.01
  ttyi 1

  cfl 0.5

  partitioning
    algorithm rcb
  end

  compflow
    depvar u
    material
      gamma 1.4 end
    end
    ic
      density 1.0 end
      velocity 0 0 0 end
      # ic pressure = (p_in + p_out) / 2
      pressure 99950 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end
    bc_pressure  # inflow
      sideset 1 end
      density 1.0
      pressure 100000
    end
    bc_pressure  # outflow
      sideset 3 end
      density 1.0
      pressure 99900
    end
  end

  field_output
    interval 10
  end

  integral_output
    #interval 1
    #time_interval 3.0e-5
    time_range 1.3e-4 2.0e-4 1.0e-5 end
    sideset 1 3 end
    format scientific
    #precision 12
  end

  diagnostics
    interval  1
    format scientific
    precision 12
  end

end
