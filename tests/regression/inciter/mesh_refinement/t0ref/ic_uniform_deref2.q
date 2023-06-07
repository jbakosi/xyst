# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Initial uniform mesh refinement"

inciter

  nstep 3     # Max number of time steps
  cfl   0.8   # CFL coefficient
  ttyi  1     # TTY output interval

  partitioning
    algorithm rcb
  end

  problem slot_cyl

  compflow
    depvar u
    material
      gamma 1.66666666666667 end
    end
    bc_dirichlet
      sideset 1 1 1 1 1 1 1 end
      sideset 2 1 1 1 1 1 1 end
      sideset 3 1 1 1 1 1 1 end
      sideset 4 1 1 1 1 1 0 end
      sideset 5 1 1 1 1 1 1 end
      sideset 6 1 1 1 1 1 0 end
    end
  end

  amr
    t0ref true
    initial ic
    initial uniform_derefine
    initial uniform_derefine    # no-op
    initial uniform_derefine    # no-op
    initial ic
    initial uniform
    refvar u6 end
    error hessian
  end

  field_output
    interval 1
  end

end
