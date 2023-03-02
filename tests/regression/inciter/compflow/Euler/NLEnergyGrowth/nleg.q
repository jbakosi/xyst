# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing nonlinear energy growth"

inciter

  term 1.0
  ttyi 1
  cfl 0.8

  partitioning
    algorithm mj
  end

  problem nonlinear_energy_growth
  compflow
    depvar c
    alpha 0.25
    betax 1.0
    betay 0.75
    betaz 0.5
    r0 2.0
    ce -1.0
    kappa 0.8
    material
      gamma 1.66666666666667 end
    end
    bc_dirichlet
      sideset 1 2 3 4 5 6 end
    end
  end

  field_output
    interval 5
  end

  diagnostics
    interval  1
    format scientific
  end

end
