function out = var_sample_covar(p1,p2,p12,n)
% calculate the variance of the sample covariance
% variance given by p12 - p1*p2 with n samples
% derivation in dropbox/reports/variance of covariance in johns sce paper
% main calling function in tes851_bernoulliTestv4

firstMu_ana = (n*p12 + n*(n-1)*(p12^2))/(n^2);

secondMu_ana = (n*(n-1)*(n-2)*p12*p1*p2 + n*(n-1)*p12*(p1+p2) + n*(n-1)*(p12^2)  + n*p12)/(n^3);

thirdMu1 = (n-3)*(n-2)*(n-1)*n*(p1^2)*(p2^2);
thirdMu2 =  (n-2)*(n-1)*n*p1*(p2^2) + (n-2)*(n-1)*n*(p12)*p1*p2;
thirdMu3 = (n-2)*(n-1)*n*p12*p1*p2 + (n-2)*(n-1)*n*p1*p2*p12;
thirdMu4 = (n-2)*(n-1)*n*(p1^2)*p2 + (n-2)*(n-1)*n*p12*p1*p2;
thirdMu5 = 2*(n-1)*n*p12*p2 + 2*(n-1)*n*p12*p1;
thirdMu6 = n*p12;
thirdMu7 = n*(n-1)*p1*p2 + n*(n-1)*(p12^2) + n*(n-1)*(p12^2);
thirdMu_ana= (thirdMu1+(thirdMu2+thirdMu3+thirdMu4)+thirdMu5 + thirdMu6+thirdMu7)/(n^4);

out = firstMu_ana-2*secondMu_ana+ thirdMu_ana - (p12-p1*p2)^2;


%%%%%%%%%
% 
% pij=p12;
% pi=p1;
% pj=p2;
% zeroOrder = 6*pij - 4*pi*pij - 4*pij^2 + 2*pi*pj + 4 pi^2*pj - 4 pij*pj + 12*pi*pij*pj + 4*pi*pj^2 - 12*pi^2*pj^2; 
% firstOrder = (-10*pij + 8*pi*pij + 8*pij^2 - 4*pi*pj - 6*pi^2*pj + 8*pij*pj - 26*pi*pij*pj - 6*pi*pj^2 + 22*pi^2*pj^2)*n; 
% secondOrder = (2 pij - 4*pi*pij - 4*pij^2 + 2*pi*pj + 2*pi^2*pj - 4*pij*pj + 18*pi*pij*pj + 2*pi*pj^2 - 12*pi^2*pj^2)*n^2; 
% thirdOrder = - 2*pij^2*n^3; 
% fourthOrder = (pij - pij^2)*n^4; 
% fifthOrder = pij^2*n^5;