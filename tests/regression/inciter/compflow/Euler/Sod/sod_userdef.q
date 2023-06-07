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

  #problem sod

  compflow
    depvar u

    ic
      density -1.0 end                  # overwritten by boxes
      velocity 100.0 100.0 100.0 end    # overwritten by boxes
      pressure -1.0 end                 # overwritten by boxes
      box
        xmin -0.5 xmax 0.5
        ymin -0.5 ymax 0.5
        zmin -0.5 zmax 0.5
        density 1.0
        pressure 1.0
      end
      box
        xmin  0.5 xmax 1.5
        ymin -0.5 ymax 0.5
        zmin -0.5 zmax 0.5
        density 0.125
        pressure 0.1
      end
    end

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
    interval  1
    point p1 0.1 0.05 0.025 end
    point p2 0.9 0.05 0.025 end
  end

  diagnostics
    interval 1
    format scientific
  end

end
