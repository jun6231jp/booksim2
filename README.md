Customized BookSim Interconnection Network Simulator
=========================================

BookSim is a cycle-accurate interconnection network simulator.
Originally developed for and introduced with the [Principles and Practices of Interconnection Networks](http://cva.stanford.edu/books/ppin/) book, its functionality has since been continuously extended.

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
For a simlulation for 4DxF7 PolarFly+ with random traffic, set the config file as follows:
```

// Topology
topology = polarflyplus;
hypercubeport = 4; //hypercube port
polarflyport = 8; //polarfly port F7:8, F5:6, F3:4, F2:3
nic = 1; 
seq = 3;
//chunk = 114;
// Routing
routing_function = dim_order;
credit_delay   = 1; //TX/RX XB (1GHz 1cycle=1ns) 
routing_delay  = 1; //TX/RX XB 
vc_alloc_delay = 1; //TX/RX XB
sw_alloc_delay = 1; //TX/RX XB
st_final_delay = 10; //XB
st_prepare_delay = 10;//XB
input_speedup     = 1; //TX NIC inject rate 1 flit per 1 cycle
output_speedup    = 1; //RX NIC eject rate 1 flit per 1 cycle
internal_speedup  = 1.0; //intermediate router speed 1cycle 1flit
num_vcs = 6;
vc_allocator = separable_input_first;
sw_allocator = separable_input_first;
vc_buf_size = 800;

// Traffic
traffic = uniform;
//latency: drains all packet, throughput:no drain?
sim_type = latency;
warmup_periods = 300;
sim_count      = 1;
sample_period  = 300;
max_samples    = 1000;
//collective_threshold = 2;//1;
//delay_threshold = 1;

//1: batch mode, 0: injection mode
use_read_write = 0;

//for injection mode
packet_size = 4;
injection_rate = 0.25;

//for batch mode
read_request_size = 1;
write_request_size = 1;
read_reply_size = 1;
write_reply_size = 1;

//link_failures = 1;
//fail_seed=time;
```
  
