NPBayesHMM : Nonparametric Bayesian HMM toolbox, for Matlab

Website: http://michaelchughes.github.com/NPBayesHMM/

Author:  Mike Hughes (www.michaelchughes.com, mike AT michaelchughes.com)

This toolbox provides code for running Markov chain Monte Carlo (MCMC) posterior inference for the Beta Process (BP) hidden markov model (HMM).

This software is released under the BSD 3-clause license. Please see the LICENSE file for details.

*Update 2016:* This repo is not under active development, but I'm still happy to handle support requests.

# Organization

The repository is organized as follows:  
  code/ contains relevant Matlab code. This should be the working dir in Matlab.
        within code/, you can find a fast intro script in code/demo/EasyDemo.m
  doc/  contains human-readable documentation.
        QuickStartGuide.md should be your one-stop resource for getting started.
  data/ contains one example dataset (6 Mocap sequences of various exercises)
        See the demos for how to run posterior inference on this data. 
        Other example datasets (from our NIPS 2012 paper) are available by contacting Mike via email.
      
# Academic Citation

If you find this toolbox useful, please cite one of our papers:

NIPS 2012 paper, which 
M. Hughes, E. Fox, and E. Sudderth. 
"Effective Split-Merge Monte Carlo Methods for Nonparametric Models of Sequential Data". 
Advances of Neural Information Processing Systems (NIPS) 2012.
http://web.michaelchughes.com/research/bp-hmm-split-merge-inference/HughesFoxSudderth_NIPS2012.pdf

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
