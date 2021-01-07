# a julia implementation of Random Projection Tree Classifier

The implementation follows the papers of Das Gupta and Freund.

## Algorithm

The data are first stored in a root node of a tree.
Then a binary tree is created by propagating data down to leaves according to random projection and a sphericity constrains inside nodes.
We thus get a collection of $2^{d}$ leaves if d is the depth of the tree asked for.
  
## License

Licensed under either of

* Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>

* MIT license  <http://opensource.org/licenses/MIT>

at your option.

## References

1. **DasGupta S. Freund Y. Random projection trees for vector quantization 2009**.
2. **DasGupta S. Freund Y. Random projection trees and low dimensional manifolds 2007**.
