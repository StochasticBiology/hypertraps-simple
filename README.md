# hypertraps-simple
Simple, low-level HyperTraPS code for projects

Summary 
---

HyperTraPS (hypercubic transition path sampling) is an algorithm for inferring the likely dynamic pathways that give rise to a series of observations. These observations consist of binary strings describing the presence (1) or absence (0) of each of L traits of interest. HyperTraPS can broadly be thought of as a Hidden Markov Model where the state space is the L-hypercube where each vertex is one of these strings and an edge exists between vertices whose strings differ in exactly one place. 

If you use HyperTraPS please cite the original paper Johnston & Williams, Cell Systems 2, 101-111 (2016), and/or the methods paper Greenbury et al., Cell Systems 10, 39-51 (2020).

Getting Started 
---

If you're running a Linux-like command line with gcc, the quickest way to get HyperTraPS running is to use the `demo.sh` script. You may need to mark it as executable:
`chmod +x demo.sh`

then run it with
`./demo.sh`

This compiles and runs the source code, and runs a quick instance of HyperTraPS on an example synthetic data file. If R with ggplot2 is installed it will also produce a summary plot.

IMPORTANT NOTE: this simple demo script performs NO checking about the performance of HyperTraPS. This is a really simple test file and the parameters used will not in general be sufficient to get converged results for bigger/more complex data. The script `demo-parallel.sh` runs several chains in parallel and visualises the corresponding output as a first check for convergence; better approaches like Gelman-Rubin are recommended in more formal settings.

Running Windows? You'll need to compile the source code below and run it via the command line; unfortunately this Bash script won't work in Windows.

Code and I/O 
---

`hypertraps.c` is C source code for HyperTraPS inference.
`posteriors.c` is C source code for the subsequent analysis of HyperTraPS posteriors.

In Linux, at the command line, compile both with gcc and link the math library:
`gcc -o3 hypertraps.c -lm -o hypertraps.ce`
`gcc -o3 posteriors.c -lm -o posteriors.ce`

Invoke HyperTraPS with
`./hypertraps.ce [datafile] [random number seed] [length index] [kernel index] [direction index]`

for example
`./hypertraps.ce fake-cross-samples-1.txt 1 1 5 0`

This will produce an output file `[datafile]-posterior-0-[random number seed]-[length index]-[kernel index].txt`. This file contains samples from the posterior distribution that HyperTraPS has learned from the data. Use the analysis code on this with
`./posteriors.ce 0 [output file]`

for example
`./posteriors.ce 0 fake-cross-samples-1.txt-posterior-0-1-1-5.txt`

This will produce a summary of the dynamics in file `[output file].process` (among others).

For Windows, use a C compiler and the command prompt to issue these commands. I'm not sure if `./` works in Windows; best to navigate in the command prompt to the location with the executable and run locally?

Command-line arguments 
---

In addition to the data file name, HyperTraPS takes these command-line parameters:

* [random number seed] -- an integer initialising the random number generator. Use different values to run parallel MCMC chains and check for convergence.
* [length index] -- the length of the MCMC simulation. HyperTraPS will run for 10^([length index]+2) steps, so 1 = 1000, 2 = 10000 etc.
* [kernel index] -- gives the kernel used to obtain new parameter values in MCMC. Index 1 adds values drawn from N(0, 0.005) to 10% of parameters. Indices 2-7 add normal variates to all parameters, with mean 0 and s.d 0.05, 0.05, 0.1, 0.25, 0.5, 0.75 respectively.
* [direction index] -- whether the system acquires traits starting from 0000... (index 0), or loses traits starting from 1111... (index 1).

Data and input files 
---

HyperTraPS observations can be independent (taken from independent instances of some dynamic system) or linked (where observations of the same system are made over time). This latter category encompasses observations that are longitudinal (a single system is observed at different timepoints) and phylogenetic (a system is observed over time, during which it may branch and lead to more systems that subsequently evolve independently).

Example independent data:
* Patient A: 1001010
* Patient B: 1000101
* Patient C: 0110101

Example longitudinal data:
* Patient A: 0000100 -> 0110100 -> 1110101
* Patient B: 0010001 -> 0110001 -> 0110001

Example phylogenetic data:
* Leaf A: 1001000
* Leaf B: 1001001
* Leaf C: 0000001
* Phylogeny ((A,B),C)

HyperTraPS requires data in the form of pairs of observations -- a "before" and an "after" state. For independent data, the "before" state is taken to be 000... for each pair, as each instance of the system independently starts from the same initial state. For longitudinal data, the pairs are the ith and (i+1)th observations from a longitudinal chain. For phylogenetic data, the pairs are ancestor and descendant vertices on the phylogeny. An independent algorithm must be used to infer the states associated with ancestral nodes. A convenient choice for this is the "maximum parsimony" approach, where traits present in all an ancestor's descendants are assumed to also have been present in the ancestor.

Hence, for our examples -- example independent data:
* 0000000 -> 1001010
* 0000000 -> 1000101
* 0000000 -> 0110101`

Example longitudinal data:
* 0000000 -> 0000100
* 0000100 -> 0110100
* 0110100 -> 1110101
* 0000000 -> 0010001
* 0010001 -> 0110001
* 0110001 -> 0110001 (NB)

Example phylogenetic data:
* 0000000 -> 0000001
* 0000000 -> 1001000
* 1001000 -> 1001000 (NB)
* 1001000 -> 1001001 

(NB) note that "before" and "after" states here are identical, indicating an absence of any transition. Such observations are not used by HyperTraPS in discrete time (but are used in continuous-time HyperTraPS).

The data file input to HyperTraPS has two options:

_Unlabelled_ -- contains these transitions arranged in pairs of rows, so that each odd-numbered row contains a "before" state and the subsequent even-number row contains the corresponding "after" state. Presence/absence markers are separated by spaces. So our first example would be 
```
0 0 0 0 0 0 0  
1 0 0 1 0 1 0  
0 0 0 0 0 0 0  
1 0 0 0 1 0 1  
0 0 0 0 0 0 0  
0 1 1 0 1 0 1
```

_CSV format_ -- comma-separated value file with a header line and exactly these columns: Label of "before" observation | Label of "after" observation | the L traits of the "before" observation | the L traits of the "after observation". For example, with a simpler example than above for brevity:
```
Ancestor,Descendant,A.trait1,A.trait2,A.trait3,D.trait1,D.trait2,D.trait3
Node_1,Node_2,0,0,0,0,0,1
Node_2,Node_3,0,0,1,0,1,1
```

Output files 
---

The output of `posteriors.ce` includes a simple summary of the dynamics inferred by HyperTraPS: the probability that each feature is acquired/lost at each possible ordering. There are 5 columns in the `[output file].process` file:
`[ordering index] [feature index when sorted by mean ordering] [original feature index] [feature label] [probability]`

So, for example, plotting column 3 horizontally, column 1 vertically, and column 5 as colour or point size will give a heatmap-style summary plot of the inferred dynamics.

If you call `posteriors.c` with 1 as the first argument (as opposed to 0 above), there is also a `[output file]-routes.txt` file which stores samples of individual trajectories supported by the inferred parameters. Here, each row is an ordered list of the features that are acquired/lost over one simulation of the system.

Please note this is a quick-and-dirty version -- a more official release can be found here https://github.com/StochasticBiology/HyperTraPS .
