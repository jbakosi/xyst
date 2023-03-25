[![status-badge](https://ci.codeberg.org/api/badges/xyst/xyst/status.svg)](https://ci.codeberg.org/xyst/xyst)

_Xyst_ is a Navier-Stokes solver for engineering flows.

Our ultimate goal is to simulate engineering problems with a production-quality
code that is extensible and maintainable, using hardware resources efficiently,
even for problems with a priori unknown, heterogeneous, and dynamic load
distribution.

The software implementation facilitates the effective use of hardware of any
size, from laptops to the largest distributed-memory clusters, by combining
data-, and task-parallelism on top of the Charm++
(\url{http://charmplusplus.org}) runtime system. Charm++'s execution model is
asynchronous by default, allowing arbitrary overlap of computation and
communication. Built-in automatic load balancing enables redistribution of
heterogeneous computational load based on real-time CPU load measurement at
negligible cost. The runtime system also features automatic checkpointing,
fault tolerance, resilience against hardware failure, and supports power-, and
energy-aware computation.

More details at [xyst.cc](https://xyst.cc).
