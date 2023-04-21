# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Dispersion in street canyon from a point source"

inciter

  nstep 50    # Max number of time steps
  #term 1000.0     # Max time
  cfl 0.5
  ttyi 10       # TTY output interval

  partitioning
    algorithm rcb
  end

  #steady_state true
  #residual 1.0e-8
  #rescomp 6

  problem point_src

  compflow
    depvar u
    ic
      density 1.225 end
      pressure 1.0e+5 end
      velocity 0.0 0.0 0.0 end
    end
    material
      gamma 1.4 end
    end
    bc_sym
      sideset 1 2 3 4 5 end
    end
    bc_pressure # inlet
      density 1.225
      pressure 1.0e+5
      sideset 6 end
    end
    bc_pressure # outlet
      density 1.225
      pressure 0.9e+5
      sideset 7 end
    end
    source
      3.0 0.01 0.0 # location
      0.2          # radius
      0.0          # release time
    end
  end

  field_output
    interval 10
   end

  history_output
    #interval 1
    #time_interval 4.0e-3
    time_range 15.0 20.0 1.0e-5 end
    point p1 0.5 2.001 end  # trigger error: must be 3 reals
  end

  diagnostics
    interval 10
    format scientific
  end

end
