-- vim: filetype=lua:

print "Sod shocktube"

nstep = 10
term = 0.2
ttyi = 1
cfl = 0.5

solver = "zalcg"
fctclip = true
fctsys = { 1, 2, 5 }

part = "rcb"

-- problem = {
--   name = "sod"
-- }

ic = {
  density = -1.0,                     -- overwritten by boxes
  velocity = { 100.0, 100.0, 100.0 }, -- overwritten by boxes
  pressure = -1.0,                    -- overwritten by boxes
  boxes = {
    { x = { -0.5, 0.5 },
      y = { -0.5, 0.5 },
      z = { -0.5, 0.5 },
      density = 1.0,
      pressure = 1.0,
      velocity = { 0, 0, 0 }
    },
    { x = {  0.5, 1.5 },
      y = { -0.5, 0.5 },
      z = { -0.5, 0.5 },
      density = 0.125,
      pressure = 0.1,
      velocity = { 0, 0, 0 }
    }
  }
}

mat = { spec_heat_ratio = 1.4 }

bc_sym = {
  sideset = { 2, 4, 5, 6 }
}

bc_dir = {
 { 1, 1, 1, 1, 1, 1 },
 { 3, 1, 1, 1, 1, 1 }
}

fieldout = {
  iter = 10000
}

histout = {
  iter = 1,
  points = {
    { 0.1, 0.05, 0.025 },
    { 0.9, 0.05, 0.025 }
  },
  precision = 6
}

diag = {
  iter = 1,
  format = "scientific"
}
