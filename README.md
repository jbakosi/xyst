[![status-badge](https://ci.codeberg.org/api/badges/xyst/xyst/status.svg)](https://ci.codeberg.org/xyst/xyst)

_Xyst_ is a Navier-Stokes solver for complex domains. Using the
[Charm++](http://charmplusplus.org/) runtime system, it employs _asynchronous_
(or non-blocking) parallel programming and decomposes computational problems
into a large number of work units (that may be more than the available number
of processors) enabling _arbitrary overlap_ of parallel computation,
communication, input, and output. Then the runtime system _dynamically_ and
_automatically_ homogenizes computational load across the simulation
distributed across many computers.

Its ultimate goal is to simulate large and complex engineering problems with a
production-quality code that is extensible and maintainable, using hardware
resources efficiently, even for problems with _a priori_ unknown,
heterogeneous, and dynamic load distribution.

More details at [xyst.cc](https://xyst.cc) or [doc/pages/mainpage.dox](doc/pages/mainpage.dox).
