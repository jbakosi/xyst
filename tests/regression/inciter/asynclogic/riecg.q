# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Asynclogic test"

inciter

  nstep 10    # Max number of time steps
  dt   0.001  # Time step size
  ttyi 1      # TTY output interval

  compflow
    depvar u
    ic
      density 1.0 end
      pressure 1.0 end
      velocity 0 0 0 end
    end
    material
      gamma 1.4 end
    end
  end

end
