## Overview

This repository contains 

1. Implementation (in MATLAB) of the MPF-BML framework, an algorithm which infers the parameters of the Maximum Entropy distribution with the Potts model.  This is used to infer the fitness landscape of gp160 based on its sequence data. gp160 is a protein in HIV which is the primary target of antibodies.  
2. Preprocessed mulitple sequence alignment (MSA) of gp160 (fasta format) and
3. Gp160 landscape (MATLAB .mat format)

as described in 

RHY Louie, KJ Kaczorowski, JP Barton, AK Chakraborty, MR McKay, "The fitness landscape of the Human Immunodeficiency Virus envelope protein that is targeted by antibodies", 2017

## Installation

To run the MPF and BML components of the framework, there are two C MEX files in the "Helper Functions" folder which need to be built. To install

1. Open MATLAB
2. Change to the "Helper Functions" folder.
3. Type `mex K_dK_MPF.c` and `mex gibbs_potts_mex.c`

## Details and usage of the MPF-BML implementation

The MPF-BML computational framework is an algorithm to infer the field and coupling parameters of the Maximum Entropy distribution.  The code to get started is 

`main_MPF_BML.m`

which runs the complete framework, and plots various statistics to confirm the inferred parameters. The code has been deliberately left as a script, not a function, to allow users  to explore the different steps of the framework. Example data is provided. The framework comprises of three steps, each of which can be run independently of the other:

#### (1) Mutant Combining

The purpose of this step is to reduce the number of states (resulting in a decrease in the number of couplings)  to achieve a balance between bias and variance. 

##### Usage

`phi_opt = mutantCombining(msa_aa, weight_seq);`

where the inputs are:

msa_aa - a matrix of characters, with each row representing a sequence of observed states (or in the context of the paper, the aminoi acid multiple-sequence-alignment (MSA)) 

weight_seq - the weighting of each sequence,

and the outputs:

phi_opt -  optimal combining factor  which represents the fraction of the entropy obtained by "coarse-graining" or combining the least-frequent states to one state, compared to the entropy without combining. Note that phi_opt=0 corresponds to the pure Ising case, while phi_opt=1 corresponds to the pure Potts case.

#### (2) MPF

This step runs a regularized Potts and mex-function extension of the Minimum-Probability-Flow (MPF) algorithm, as originally proposed in 

Sohl-Dickstein J, Battaglino P, DeWeese MR (2009) Minimum Probability Flow learning. Proc 28th ICML 107(Ml):12.

##### Usage

`J_MPF = MPF_run(msa_bin_unique,weight_seq_unique,num_mutants_combine_array,phi_opt,options_MPF);`

where the inputs are:

msa_bin_unique  - a binary potts extension of the original MSA (which is produced by the provided function `binMatAfterComb`)

weight_seq_unique  -  the weighting of each sequence in msa_bin_unique (which is produced by the provided function `binMatAfterComb`)

num_mutants_combine_array  -  the number of mutants at each residue after mutant combining (which is produced by the provided function `binMatAfterComb`)

phi_opt  - the mutant combining factor obtained from step 1. This can be manually set by the user, though this should also be manually set in the `binMatAfterComb` function.

The output is fields/couplings matrix. Note that the non-diagonal elements are a factor of 1/2 the true couplings, thus the energy of sequence x is calculated as x' J_MPF x.

#### (3) BML

This step implements the RPROP algorithm to  refine the parameters inferred from MPF.

##### Usage

`J_MPF_BML =BML_run(J_MPF(:),msa_bin_unique,weight_seq_unique,num_mutants_combine_array,options_BML);`

where the inputs are as described in Step 2, with the exception of the first argument J_MPF(:), which is a flattened fields/couplings vector which initalizes the BML algorithm.

## Details of the MSA

The processed MSA (as described in the paper) in fasta format `hivgp160_processed_MSA.fasta` is in the folder "MSA and Landscape". The weighting of each sequence is in `hivgp160_patient_weighting.mat`.

## Details of the Landscape

The gp160 field and coupling (landscape) parameters are in `hivgp160_landscape.mat`in the folder "MSA and Landscape", where J_MPF_BML is the field/coupling matrix (NB the off-diagonal entries are half the value of the true couplings, thus the energy of sequence x is calculated as x' J_MPF_BML x), amino_acid_after_combining is the amino acids in decreasing order of frequency at each of the 815 residues, and mut_mat is the mutant probaility matrix (after mutant combining).


Any questions or comments, please email raylouie@hotmail.com
