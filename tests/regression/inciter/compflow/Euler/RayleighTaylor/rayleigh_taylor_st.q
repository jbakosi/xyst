# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing the stationary Rayleigh-Taylor MMS problem"

inciter

  nstep 10
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
    kappa 0.0   # kappa = 0 -> stationary
    material
      gamma 1.66666666666667 end
    end
    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  field_output
    interval 1
  end

  diagnostics
    interval  1
    format scientific
  end

end
