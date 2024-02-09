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

zk = 1.1;
zpars = complex([zk, 1.0, 0.0]);
ndeg = 2;

jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);

% Evaluate spherical harmonic on grid
f = spherefun.sphharm(ndeg,1);

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = f(S.r(1,:)./rr,S.r(2,:)./rr,S.r(3,:)./rr);


% Compute layer potential with spherical harmonic as data
rhs = rhs(:);
tol = 1e-7;
p = helm3d.dirichlet.eval(S,zpars,rhs,tol);

% Compare to exact solution
p_ex = rhs*jn*hn*1j*zk;

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));

fprintf('Error in single layer potential=%d\n',err1);

%% Test boundary value solve
xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
rhs = helm3d.kern(zk,src_info,S,'s');


zpars = [zk,-1j*zk,1];
sig = helm3d.dirichlet.solver(S, zpars, rhs, tol);

targ_info = [];
targ_info.r = xyz_out;

dat = helm3d.kern(zk,S,targ_info,'c',zpars(2),zpars(3));

pot = dat*(sig.*wts);
pot_ex = helm3d.kern(zk,src_info,targ_info,'s');
fprintf('Error in iterative solver=%d\n',abs(pot-pot_ex)/abs(pot_ex));


