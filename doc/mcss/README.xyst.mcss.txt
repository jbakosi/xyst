Updates to <xyst>/doc/mcss

* June 14, 2023 - initial bundle

  This directory is a copy of doxygen/ and pelican-plugins/ folders from MCSS's
  clone at 91ff035 + COPYING + CREDITS + README.

* Jun 3, 2024 - update m.css based on https://github.com/mosra/m.css 5235066

  cp -r <m.css>/documentation <xyst>/doc/mcss/doxygen
  rm -rf <xyst>/doc/mcss/doxygen/test*
  cp -r <m.css>/plugins <xyst>/doc/mcss/plugins
  rm -rf <xyst>/doc/mcss/plugins/m
  cp <m.css>/COPYING <xyst>/doc/mcss/COPYING
  cp <m.css>/CREDITS.rst <xyst>/doc/mcss/CREDITS.rs
  cp <m.css>/README.rst <xyst>/doc/mcss/README.rst
  Update this file
