# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing the non-stationary Rayleigh-Taylor MMS problem"

inciter

  nstep 50
  ttyi 1
  cfl 0.5

  partitioning
    algorithm mj
  end

  problem rayleigh_taylor
  compflow
    depvar c
    alpha 1.0
    betax 1.0
    betay 1.0
    betaz 1.0
    p0 1.0
    r0 1.0
    kappa 1.0   # kappa = 0 -> stationary
    material
      gamma 1.66666666666667 end
    end
    bc_dirichlet
      sideset 1 1 1 1 1 1 end
      sideset 2 1 1 1 1 1 end 
      sideset 3 1 1 1 1 1 end
      sideset 4 1 1 1 1 1 end
      sideset 5 1 1 1 1 1 end
      sideset 6 1 1 1 1 1 end
    end
  end

  field_output
    time_range 0.1 0.15 1.0e-3 end
  end

  diagnostics
    interval  1
    format scientific
  end

end
