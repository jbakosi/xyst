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
    point p1 0.5 2.001 0.0 end
    point p2 1.0 2.001 0.0 end
    point p3 1.5 2.001 0.0 end
    point p4 2.001 2.0 0.0 end
    point p5 2.001 1.93 0.0 end
    point p6 2.001 1.5 0.0 end
    point p7 2.001 1.33 0.0 end
    point p8 2.001 1.0 0.0 end
    point p9 2.001 0.67 0.0 end
    point p10 2.001 0.5 0.0 end
    point p11 2.001 0.33 0.0 end
    point p12 2.001 0.17 0.0 end
    point p13 3.999 0.17 0.0 end
    point p14 3.999 0.33 0.0 end
    point p15 3.999 0.5 0.0 end
    point p16 3.999 0.67 0.0 end
    point p17 3.999 1.0 0.0 end
    point p18 3.999 1.33 0.0 end
    point p19 3.999 1.5 0.0 end
    point p20 3.999 1.93 0.0 end
    point p21 3.999 2.0 0.0 end
    point p22 4.5 2.001 0.0 end
    point p23 5.0 2.001 0.0 end
    point p24 5.5 2.001 0.0 end
  end

  diagnostics
    interval 10
    format scientific
  end

end
