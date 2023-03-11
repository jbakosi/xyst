# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing vortical flow - restarted"

inciter

  nstep 10
  ttyi 1
  cfl 0.8

  partitioning
   algorithm mj
  end

  problem vortical_flow

  compflow
    depvar c
    alpha 0.1
    beta 1.0
    p0 10.0
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
    interval 1
    format scientific
  end

end
