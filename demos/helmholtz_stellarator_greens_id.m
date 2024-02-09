%
%  This code illustrates loading of surfaces, plotting them
%  and testing green's identity at far away points on a sphere
% 
%  Todo: make it into a stellarator
addpath(genpath('~/git/fmm3dbie/matlab'))

% Load sphere mesh
% norder is the order of discretization
% nu is the number of initial triangles/quads on each face of the cube
% iref is the additional number of refinements
% iptype is the patch type = 1 for triangles and 11/12 for quads

norder = 6;
nu = 1; 
iref = 2;
iptype = 11;
S = surfer.sphere(norder, nu, iref, iptype);

tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;

% Plot sphere
tic, plot(S); toc;


%% Test greens identity off surface
%
%  u = D[u] - S[dudn]  
%
zk = 1.1;
xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
u = helm3d.kern(zk, src_info, S, 's');
dudn = helm3d.kern(zk, src_info, S, 'sprime');


targ_info = [];
targ_info.r = xyz_out;
g_bdry_to_targ = helm3d.kern(zk, S, targ_info, 's');
dgdn_bdry_to_targ = helm3d.kern(zk, S, targ_info, 'd');



uex = helm3d.kern(zk, src_info, targ_info, 's');
utest = dgdn_bdry_to_targ*(u.*wts(:)) - g_bdry_to_targ*(dudn.*wts(:));

fprintf('Error in Greens identity=%d\n',abs(utest - uex)/abs(uex));


%%  Test greens identity on surface
%
d_u = helm3d.dirichlet.eval(S, [zk, 0, 1], u, eps) + u/2;
s_dudn = helm3d.dirichlet.eval(S, [zk, 1, 0], dudn, eps);

fprintf('Error in greens identity on surface=%d\n', norm(u - d_u + s_dudn));



