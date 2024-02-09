%  The goal is to compute the in-plane bistatic radar cross
%  section for a boxed torus geometry
%
%
%
addpath(genpath('~/git/fmm3dbie/matlab'))

% Load boxed torus mesh
% norder is the order of discretization
% iref is the additional number of refinements


% Plot boxed torus

zk = 1.1;

%% Precompute boundary correction
%

% 
%
nrhs = 100;
%% Solve with precomputed quadrature corrections
for i=1:nrhs
   % Compute plane wave boundary data with 
   % in-plane angle (\theta), i.e. 
   % incident wave coming from (\cos(theta), \sin(theta), 0)

   % Solve boundary value problem


   % Evaluate bistatic radar cross section back at reciever
   % given by
   % \int_{\Gamma} \sigma (signature)

end

% plot the bi-static radar cross-section


