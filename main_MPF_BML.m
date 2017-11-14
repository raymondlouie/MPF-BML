close all
clear all;

% Inputs:
% msa_aa: matrix of characters (aka amino acid MSA). Rows correspond to 
%         sequences, and colums to observed states at a particular residue
%
% Optional input:
% weight_seq: weight of each sequence, length is equal to the the number of
%             rows in msa_aa. Default=equal weighting.
% phi_opt   : Mutant combining factor as described in paper.
%             phi_opt=0 is Ising, phi_opt=1 is Potts. 
%             Default=output of combining algorithm as described in paper
% 
% Outputs:
% J_MPF_BML_mat - field and coupling matrix (max entropy parameters). 
% freq_single_combine_array - frequency of each amino acid after
%                         combining with factor phi_opt.
% amino_single_combine_array - amino acid sorted in decreasing order
%                         of frequency after combining with factor phi_opt
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of this framework should cite:
%
%     RHY Louie, KJ Kaczorowski, JP Barton, A Chakraborty, MR McKay
%     (2017), The fitness landscape of the Human Immunodeficiency Virus 
%     envelope protein that is targeted by antibodies, 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load example data. This data conatins the processed MSA and sequence 
% weighting (by patients) as described in the PNAS paper. Comment the below
% out if there is user-provided "msa_aa" and "weight_seq"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_example=1; % load the test example
if test_example==1
    fasta_name = 'hivgp160_processed_MSA.fasta';
    load hivgp160_patient_weighting
    [Header_fasta, Sequence_fasta] = fastaread(fasta_name);
    msa_aa = cell2mat(Sequence_fasta');

    % For testing purposes, reduce the number of sequences and residues
    % for faster computation
    
    num_seq_test=1000; % number of sequences used for testing
    num_residue_test = 70; % number of residues used in testing
    
    msa_aa = msa_aa(1:num_seq_test,1:num_residue_test); % amino acid MSA
    weight_seq = weight_seq(1:num_seq_test); % weight of each sequence
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default weight and remove 100% conserved sites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Set default weight_seq if not specified by user
% num_seq = size(msa_aa,1); % number of sequences
% if ~exist('weight_seq')
%     % set equal weighting if weighting vector not provided
%     weight_seq = ones(num_seq,1);
% end

% num_patients = sum(weight_seq); % number of patients

% Remove and find location of 100% conserved residues
num_residue = size(msa_aa,2);
ind_conserve=[];
ind_non_conserve=[];
for ind_residue = 1:num_residue
    if length(unique(msa_aa(:,ind_residue)))==1
        ind_conserve=[ind_conserve ind_residue];
    else
        ind_non_conserve = [ind_non_conserve ind_residue];
    end
end

msa_aa(:,ind_conserve) = [];

if length(ind_conserve)==1
    disp(['There is one 100% conserved residue detected and removed, resulting in ' ...
        num2str(length(ind_non_conserve)) ' number of residues.  The location of residues are in the variable ind_conserve']);
elseif length(ind_conserve) > 1
    disp(['There are ' num2str(length(ind_conserve)) ' 100% conserved residues detected and removed, resulting in ' ...
        num2str(length(ind_non_conserve)) ' number of residues. The location of residues are in the variable ind_conserve']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Mutant combining
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_step1_mutantcombine = tic();

% phi_opt = mutantCombining(msa_aa, 'phi_array',[0:0.01:1]);
phi_opt = mutantCombining(msa_aa, 'weight_seq',weight_seq);
% phi_opt = mutantCombining(msa_aa);

time_step1_mutantcombine = toc(time_step1_mutantcombine);

disp(['Step 1: Mutant combining, Time: ' num2str(time_step1_mutantcombine) ' seconds']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate step - calculate extended binary matrix, and other statistics
%                     required for MPF and RPROP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_extendBin = tic();

[msa_bin,msa_bin_unique,weight_seq_unique, freq_single_combine_array,amino_single_combine_array,num_mutants_combine_array] = ...
    binMatAfterComb(msa_aa,'weight_seq',weight_seq,'phi_opt',phi_opt);
% [msa_bin,msa_bin_unique,weight_seq_unique, freq_single_combine_array,amino_single_combine_array,num_mutants_combine_array] = binMatAfterComb(msa_aa,weight_seq,phi_opt);

% mut_mat = full(((msa_bin_unique')*diag(weight_seq_unique)*msa_bin_unique))/num_patients; % mutant probability matrix
% [delta_cij delta_cij_bound] = calculateStd(mut_mat,num_patients);

time_extendBin = toc(time_extendBin);

disp(['Intermediate step - Extended binary matrix and other stats, Time: ' num2str(time_extendBin) ' seconds']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: MPF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_step2_MPF = tic();

options_MPF.lambda_J = 0.01; % L1 regularization parameter
options_MPF.gamma_J = 0.02; % L2 regularization parameter

J_MPF = MPF_run(msa_bin_unique,weight_seq_unique,num_mutants_combine_array,phi_opt,options_MPF);

time_step2_MPF = toc(time_step2_MPF);

disp(['Step 2: MPF, Time: ' num2str(time_step2_MPF) ' seconds']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: BML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_step2_BML = tic();

options_BML.no_iterations=50;
options_BML.epsMax = 1.25;

J_MPF_BML =BML_run(J_MPF(:),msa_bin_unique,weight_seq_unique,num_mutants_combine_array,options_BML);

time_step2_BML = toc(time_step2_BML);

disp(['Step 3: BML, Time: ' num2str(time_step2_BML) ' seconds']);

num_residues_binary = size(msa_bin_unique,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verification of the landscape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = verifyParam(J_MPF(:),msa_bin_unique,weight_seq_unique,num_mutants_combine_array);

out = verifyParam(J_MPF_BML(:),msa_bin_unique,weight_seq_unique,num_mutants_combine_array);
