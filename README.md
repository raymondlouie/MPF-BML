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
2. Change directory to the "Helper Functions" folder.
3. In the command prompt, enter` mex K_dK_MPF.c`
4. In the command prompt, enter  `mex gibbs_potts_mex.c`

## Toolboxes

The following MATLAB toolboxes are required:

1. Bioinformatics
2. Communications System
3. Parallel Computing

## Details and usage of the MPF-BML implementation

The MPF-BML computational framework is an algorithm to infer the field and coupling parameters of the Maximum Entropy distribution.  An example working code is the script

`main_MPF_BML.m`

which runs the complete framework, and plots various statistics to confirm the inferred parameters. The code has been deliberately left as a script, not a function, to allow users  to explore the different steps of the framework. Example data is provided. The framework comprises of three main functions corresponding to the three key steps of the algorithm, each of which can be run independently of the other. From now on, we assume all example commands are to be run in the MATLAB command prompt.

First note that the input to each function is a sample character matrix `msa_aa`, which can be formed from a fasta file with name `fasta_name` by

```
[Header_fasta, Sequence_fasta] = fastaread(fasta_name);
msa_aa = cell2mat(Sequence_fasta');
```

#### (1) Mutant Combining

The purpose of this step is to reduce the number of states (resulting in a decrease in the number of couplings)  to achieve a balance between bias and variance. The function which implements this is `mutantCombining` and the output `phi_opt`  is the optimal combining factor  which represents the fraction of the entropy obtained by "coarse-graining" or combining the least-frequent states to one state, compared to the entropy without combining. Note that `phi_opt=0` corresponds to the pure Ising case, while `phi_opt=1` corresponds to the pure Potts case.

##### Example usage

Choose the optimal combining factor from `phi_array`, a vector of possible values, and `weight_seq`, the weighting per sequence.

```
phi_array = [0:0.1:1]; 
weight_seq = ones(size(msa_aa,1),1) ; % equal weighting per patient
phi_opt = mutantCombining(msa_aa, 'weight_seq',weight_seq,'phi_array',phi_array);
```

The default values of weight_seq is set to equal weighting per patient, i.e.,

`weight_seq = ones(size(msa_aa,1),1) ; % equal weighting per patient `

while the default value of phi_array is

`phi_array=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9:0.1:1]; `

#### Intermediate step: helper variables

The MPF (Step 2) and BML (Step 3) functions both require  helper variables, produced by the function `binMatAfterComb.m`. These helper variables are 

`msa_bin` - binary extended matrix after combining with factor phi_opt
`msa_bin_unique` - unique rows of msa_bin
`weight_seq_unique` - weight of each sequence in msa_bin_unique
`freq_single_combine_array` - frequency of each amino acid after combining with factor phi_opt.
`amino_single_combin_array` - amino acid sorted in decreasing order of frequency after combining with factor phi_opt
`num_mutants_combine_array` - number of mutants at each residue
`phi_opt` - optimal combining factor

##### Example usage

Calculate the helper variables using `msa_aa` only, in which case the default

```
phi_opt = 0; % Ising case
weight_seq = ones(size(msa_aa,1),1) ; % equal weighting per patient
[...]  = binMatAfterComb(msa_aa,'weight_seq','phi_opt');
```

The default values of weight_seq is set to equal weighting per patient, i.e.,

`weight_seq = ones(size(msa_aa,1),1) ; % equal weighting per patient `

while the default value of phi_opt is the Potts case, i.e.,

`phi_opt=1; `

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

options_MPF - an options struct file which controls various paramters of the algorithm. The most relevant parameters to tune are the regularization parameters, which can be manually set, e.g., by

```
options_MPF.lambda_J = 0.01; % L1 regularization parameter for the couplings
options_MPF.gamma_J = 0.02; % L2 regularization parameter for the couplings
```

The output is fields/couplings matrix. Note that the non-diagonal elements are a factor of 1/2 the true couplings, thus the energy of sequence `x` can be calculated as 

`x'*J_MPF*x`

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
