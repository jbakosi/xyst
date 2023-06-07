# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sod shocktube"

inciter

  nstep 10
  term 0.2
  ttyi 1
  cfl 0.5

  partitioning
    algorithm rcb
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
    sideset 2 4 5 6 end
    interval 10000
  end


  diagnostics
    interval 1
    format scientific
  end

end
