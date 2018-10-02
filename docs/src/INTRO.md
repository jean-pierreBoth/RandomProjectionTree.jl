# a julia implementation of Random Projection Tree Classifier

interesting distances are (using package Distances)
* Euclidean     L2
* Cityblock     L1
* Jaccard       fuzzy nucleotide : 1 - sum min (a,b)/ sum max(a,b) 
* CosineDist


## References

1. [ DasGupta S. Freund Y. Random projection trees for vector quantization 2009 ]
2. [ DasGupta S. Freund Y. Random projection trees and low dimensional manifolds 2007 ]
