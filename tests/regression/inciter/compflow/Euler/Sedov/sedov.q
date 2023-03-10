# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Sedov blast wave"

inciter

  nstep 10
  term 1.0

  cfl 0.5
  ttyi 1

  partitioning
    algorithm rcb
  end

  problem sedov

  compflow
    depvar u

     p0 4.86e+3        # sedov_coarse.exo,   V(orig) = 8.50791e-06
    #p0 3.63e+4        # sedov00.exo,        V(orig) = 1.13809e-06
    #p0 2.32e+5        # sedov01.exo,        V(orig) = 1.78336e-07
    #p0 1.85e+6        # sedov02.exo,        V(orig) = 2.22920e-08
    #p0 14773333.33333 # sedov02.exo+t0ref:u,V(orig) = 2.78650e-09

    material
      gamma 1.66666666666667 end
    end
    bc_sym
      sideset 1 2 3 end
    end
  end

  #amr
  # t0ref true
  # initial uniform
  #end

  diagnostics
    interval 1
    format scientific
  end

  field_output
    interval 1000
  end

end
