# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shocktube"

inciter

  term 0.2
  ttyi 1
  cfl 0.1

  partitioning
    algorithm mj
  end

  problem sod

  compflow
    depvar u
    material
      gamma 1.4 end
    end
    bc_sym
      sideset 2 4 5 6 end
    end
  end

  field_output
    interval 10000
  end

  history_output
    interval 20
    time_interval 4.0e-3
    time_range 0.1 0.15 1.0e-3 end
    time_range 0.2 0.12 0.1 end # error triggered: t0 > t1
    point p1 0.1 0.05 0.025 end
    point p2 0.9 0.05 0.025 end
    point p3 0.5 0.05 0.025 end
  end

  diagnostics
    interval 1
    format scientific
  end

end
