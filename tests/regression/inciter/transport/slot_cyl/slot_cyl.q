# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Scalar transport: slotted cylinder, cone, hump"

inciter

  nstep 20
  #term 1.57
  ttyi 10
  cfl 0.9

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
      sideset 2 1 1 1 1 1 0 end 
    end
  end

  field_output
    interval 100
  end

  diagnostics
    interval 1
    format scientific
    precision 12
  end

end
