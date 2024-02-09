% 
%  This file illustrates the computation of Helmholtz layer potentials
%  on spheres in three dimensions and compares them to analytically 
%  known solutions.
%
%  Then solve a boundary value problem due to an interior point source
%
%  Todo: convert the interior point source test to a plane wave test
%        and demo exact solution, refer to same referene pdf
%
% Additional dependencies: chebfun
%
addpath(genpath('~/git/fmm3dbie/matlab'))
addpath(genpath('~/git/chebfun'))

% Load boxed torus mesh
% norder is the order of discretization
% iref is the additional number of refinements
norder = 6;
iref = 0;

fname = ['../geometries_go3/simple_torus_o0' int2str(norder) '_r0' ...
       int2str(iref) '.go3'];
S = surfer.load_from_file(fname);

[srcvals,~,~,~,~,wts] = extract_arrays(S); 
tol = 1e-8;

% Plot sphere
plot(S); 

zk = 1.1;

%% Test boundary value solve
xyz_in = [0.51;-0.49;0.52];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk, src_info, S, 's');


zpars = [zk,-1j*zk,1];
tic, sig = helm3d.dirichlet.solver(S, zpars, rhs, tol); toc;

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk, S, targ_info, 'c', zpars(2), zpars(3));

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk, src_info, targ_info, 's');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));

%% Solve with precomputed quadrature corrections
tic, Q = helm3d.dirichlet.get_quadrature_correction(S, zpars, tol); toc
opts = [];
opts.quadrature_correction = Q;
sig_precomp = helm3d.dirichlet.solver(S, zpars, rhs, tol, opts);
fprintf('Error between densities computed with and without precomp=%d\n',...
           norm((sig - sig_precomp).*sqrt(wts(:))));

