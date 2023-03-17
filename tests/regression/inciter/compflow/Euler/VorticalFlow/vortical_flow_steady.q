# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing vortical flow"

inciter

  nstep 10
  ttyi 10
  cfl 0.8

  partitioning
   algorithm mj
  end

  steady_state true
  residual 1.0e-8
  rescomp 1

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
      sideset 1 1 1 1 1 1 end
      sideset 2 1 1 1 1 1 end 
      sideset 3 1 1 1 1 1 end
      sideset 4 1 1 1 1 1 end
      sideset 5 1 1 1 1 1 end
      sideset 6 1 1 1 1 1 end
    end
  end

  field_output
    interval 10
  end

  diagnostics
    interval 1
    format scientific
  end

end
