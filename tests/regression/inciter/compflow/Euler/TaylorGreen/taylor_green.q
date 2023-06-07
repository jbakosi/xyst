# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Euler equations computing stationary Taylor-Green"

inciter

  term 1.0
  ttyi 10
  cfl 0.8

  partitioning
    algorithm rcb
  end

  problem taylor_green
  compflow
    depvar c
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
    interval 20
  end

  diagnostics
    interval 2
    format scientific
  end

end
