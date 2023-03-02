[![status-badge](https://ci.codeberg.org/api/badges/xyst/xyst/status.svg)](https://ci.codeberg.org/xyst/xyst)

_Xyst_ is a flow solver for complex domains.

Using the [Charm++](http://charmplusplus.org/) runtime system, it employs
_asynchronous_ (or non-blocking) parallel programming and decomposes
computational problems into a large number of work units (that may be more than
the available number of processors) enabling _arbitrary overlap_ of parallel
computation, communication, input, and output. Then the runtime system
_dynamically_ and _automatically_ homogenizes computational load across the
simulation distributed across many computers.

More details at [xyst.cc](https://xyst.cc) or [doc/pages/mainpage.dox](doc/pages/mainpage.dox).
