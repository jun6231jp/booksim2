BookSim Interconnection Network Simulator
=========================================

BookSim is a cycle-accurate interconnection network simulator.
Originally developed for and introduced with the [Principles and Practices of Interconnection Networks](http://cva.stanford.edu/books/ppin/) book, its functionality has since been continuously extended.
The current major release, BookSim 2.0, supports a wide range of topologies such as mesh, torus and flattened butterfly networks, provides diverse routing algorithms and includes numerous options for customizing the network's router microarchitecture.

---
New Features
+ topology PolarFly+
  - PolarFly+ is a novel topology that combining hypercube and PolarFly.
  - F2-F7 PolarFly and F0 (No PolarFly, just a hypercube) are supported
  - 1-7D hypercube and 0D (No hypercube, just a PolarFly) are supported
+ routing algorithm
  - basic algorithm dim_order_polarflyplus
  - fault adoivance 
+ collective traffics
  - pairwise exchange
  - ring

How to use
+ random simulation
  
