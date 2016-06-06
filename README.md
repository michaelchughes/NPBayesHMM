NPBayesHMM : Nonparametric Bayesian HMM toolbox, for Matlab

Website: http://michaelchughes.github.com/NPBayesHMM/

Author:  Mike Hughes (www.michaelchughes.com, mike AT michaelchughes.com)

This toolbox provides code for running Markov chain Monte Carlo (MCMC) posterior inference for the Beta Process (BP) hidden markov model (HMM).

This software is released under the BSD 3-clause license. Please see the LICENSE file for details.

*Update 2016:* This repo is not under active development, but I'm still happy to handle support requests.

## Organization

The repository is organized as follows:

* code/

contains relevant Matlab code. This should be the working dir in Matlab. within code/, you can find a fast intro script in code/demo/EasyDemo.m

* [doc/QuickStartGuide.md](doc/QuickStartGuide.md)

This plain-text documentation file is your one-stop resource for getting started. Any further questions, please contact Mike.

* data/ 

contains one example dataset (6 Mocap sequences of various exercises)
See the demos for how to run posterior inference on this data. 
Other example datasets (from our NIPS 2012 paper) are available by contacting Mike via email.
      
## Academic Citation

If you find this toolbox useful, please cite one of our papers:


#### AOAS 2014 journal article

> This long-form article presents a coherent reference for understanding the BP-HMM and BP-AR-HMM models, as well as an updated split-merge MCMC inference algorithm.

* "Joint Modeling of Multiple Time Series via the Beta Process with Application to Motion Capture Segmentation."
Emily Fox, Michael C. Hughes, Erik B. Sudderth, Michael I. Jordan
Annals of Applied Statistics, Vol. 8(3), 2014.
[[paper]](http://michaelchughes.com/papers/FoxHughesSudderthJordan_AOAS_2014.pdf)
[[supplement]](http://michaelchughes.com/papers/FoxHughesSudderthJordan_AOAS_2014_supplement.pdf)


#### NIPS 2012 conference paper

> Our NIPS 2012 paper introduced improved split-merge inference algorithm for the BP-HMM.

* "Effective Split-Merge Monte Carlo Methods for Nonparametric Models of Sequential Data". 
Michael C. Hughes, Emily B. Fox, and Erik B. Sudderth.
NIPS 2012.
[[paper]](http://michaelchughes.com/papers/HughesFoxSudderth_NIPS_2012.pdf)
[[supplement]](http://michaelchughes.com/papers/HughesFoxSudderth_NIPS_2012_supplement.pdf)


## Acknowledgements

This code is inspired by (and heavily based upon) the
BP-AR-HMM toolbox, released by Emily Fox. Most functions have been completely 
re-written for speed, readability, and extensibility, but Emily deserves most
credit for coming up with the original solid inference algorithms.

I also thank 

* Tom Minka for his excellent Lightspeed toolbox
* The development team at Eigen for a blazingly-fast matrix library
Eigen made the HMM dynamic programming routines much much faster. 
* William Allen for providing baseline code for these efficient routines.
