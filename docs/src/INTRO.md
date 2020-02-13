# a julia implementation of Random Projection Tree Classifier

interesting distances are (using package Distances) 

* Euclidean     L2
* Cityblock     L1
* Jaccard       fuzzy nucleotide : 1 - sum min (a,b)/ sum max(a,b).
* CosineDist

## Method and Algorithm

### A rough description of the algorithm


The data are first stored in a root node of a tree.

Then the data is split in 2 descendant nodes.
There are 2 modes of splitting:

* By random sampling a random vector on the unit sphere and deciding
  to send each data vector in the left or right node depending on the sign of the cosine between the sampled unit vector and the data.

* By computing a gravity center of the points in a node and taking half the node near the
 gravity center in the one node (the left for example) and the other half in the other descendant node. This mode of splitting helps getting some sphericity of data in the final leaves.

We thus expect to get data in leaves defined by hyperplanes and having minimal sphericity.
The balance between the 2 modes of node splitting is adjusted by the **_threshold_** parameters in type RPTreeArg.  

### Implementation

The implementation relies on a minimal binary Tree type implemented in file Tree.jl.

The construction of the random projection tree **_(RPTree)_** is somewhat parallelised: 

* If ``nworkers() > 1`` the 2 nodes of the first generation are affected to 2 differents tasks.

* If ``Threads.nthreads() > 1`` the nodes of the second generations are affected to 2 different threads, thus achieving a moderate level of parallelism. (Some work to be done to better adapt this part of the code).

## Testing

The code was tested in Mass Spectrometry Imaging problems mainly with distances:

* Euclidean     L2
* Cityblock     L1
* Jaccard       fuzzy nucleotide :

    $d(a,b) = 1 - \frac{\sum_{i}{\min(a_i,b_i)}}{\sum_i{\max(a_i, b_i)}}$
    
* CosineDist

(Cf package Distances)

## References

1. **DasGupta S. Freund Y. Random projection trees for vector quantization 2009**.
2. **DasGupta S. Freund Y. Random projection trees and low dimensional manifolds 2007**.
