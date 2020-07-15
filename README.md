# CyanFlow

My power network calculate project,mainly using slow,ugly,lazy C++ code.
Using <code>Egien</code> algbra library. 

# input file format

<code>

4 3 0.00010

1 3 2 0.0600 0.2000

2 2 1 0.0300 0.0900

3 1 4 0.0300 0.0900

1 PQ -0.5000 -01500

2 PQ 0.5000 -0.150

3 PV 0.4000 0.9500

4 BL 1.0000 0.0000

</code>

## explain:

* [number of nodes] [number of lines] [precision]
* following [number of lines] lines:
* [node A number] [node B number] [p] [q] 
* indicate the line power flow from A to B.
* following [number of nodes] lines:
* [node number] [node type] [data 1] [data 2]
* data 1 &data 2 are depend on the node type:

 |node type |data 1  | data 2|
|:----:|:----:|:----:|
 |PQ|value of p|value of q|
 |PV|value of P|value of v|
|BL|value of v|angle|

# Known problem
* Only support blance node is the last node.(course limit)
* multi-thread may not be useable.
* The result may not be right.
* Potential memory leakage when free Network class.

# preformence enhance ways:
* reduce new complex.
# Linsece
GPL v3 
