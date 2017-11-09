Overview

This repository contains 

(a) Implementation of the MPF-BML framework (an algorithm which infers the parameters of the Maximum Entropy potts distribution)  
(b) Preprocessed MSA of gp160 and
(c) Gp160 landscape

as described in 

RHY Louie, KJ Kaczorowski, JP Barton, A Chakraborty, MR McKay, "The fitness landscape of the Human Immunodeficiency Virus envelope protein that is targeted by antibodies", Proc. Natl. Acad. Sci., 2017

Installation

To run the MPF and BML components of the framework, there are two C MEX files in the "Helper Functions" which need to be built: K_dK_MPF.c and gibbs_potts_mex.c. Typically

mex K_dK_MPF.c
mex gibbs_potts_mex.c

should work, however, there are many tutorials online which describe how to build MEX files.

Details and usage of the MPF-BML implementation

The computational framework is an algorithm to find the field and coupling parameters of the Maximum Entropy distribution.  
The code to get started is main_MPF_BML.m, which runs the complete framework, and plots various statistics to confirm the inferred parameters. The code has been deliberately left as a script, not a function, to allow users  to explore the different steps of the framework. Example data is provided. The framework comprises of three steps, each of which can be run independently of the other:

(1) Mutant Combining

The purpose of this step is to reduce the number of states (resulting in a decrease in the number of couplings)  to achieve a balance between bias and variance. Usage:

phi_opt = mutantCombining(msa_aa, weight_seq);

where the inputs are:

msa_aa - a matrix of characters, with each row representing a sequence of observed states (or in the context of the PNAS paper, the aminoi acid multiple-sequence-alignment (MSA)) 
weight_seq - the weighting of each sequence,

and the outputs:

phi_opt -  optimal combining factor  which represents the fraction of the entropy obtained by "coarse-graining" or combining the least-frequent states to one state, compared to the entropy without combining. Note that phi_opt=0 corresponds to the pure Ising case, while phi_opt=1 corresponds to the pure Potts case.

(2) MPF

This step runs a regularized Potts and mex-function extension of the Minimum-Probability-Flow (MPF) algorithm, as originally proposed in 

Sohl-Dickstein J, Battaglino P, DeWeese MR (2009) Minimum Probability Flow learning. Proc 28th ICML 107(Ml):12.

Usage:

J_MPF = MPF_run(msa_bin_unique,weight_seq_unique,num_mutants_combine_array,phi_opt,options_MPF)

where msa_bin_unique is a binary potts extension of the original MSA (which is produced by the provided function binMatAfterComb), weight_seq_unique is the weighting of each sequence in msa_bin_unique, num_mutants_combine_array is the number of mutants at each residue after mutant combining, phi_opt is the mutant combining factor obtained from step 1 (this can be manually set by the user, though this should also be manually set in the binMatAfterComb function)

Any questions or comments, please email raylouie@hotmail.com
