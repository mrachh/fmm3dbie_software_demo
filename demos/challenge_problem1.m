%% Challenge problem
%
%  The challenge problem is split into 3 parts:
%
%  Part 0: Construct a geometry consisting of 6 ellipsoids on 
%          the vertices of the unit hexagon 
%          (cos(2*pi*i*j/6), sin(2*pi*i*j/6), 0), j = 0,1,2,3,4,5.
%          The ellipsoids should have semi major axes (0.3, 0.25, 0.35),
%          (0.25, 0.35, 0.3), (0.35, 0.3, 0.25), (0.3, 0.25, 0.35), 
%          (0.25, 0.35, 0.3), (0.35, 0.3, 0.25). 
%
%  Part 1: Demonstrate accuracy of your transmission solver, either
%          using a self-convergence test or using an analytic solution
%          test
%
%  Part 2: For the given configuration of spheres in section 2, 
%          find the index of refraction, or alternatively the interior
%          wave number such that the L2 norm of the potential in box
%          of size () centered at (). 
%          The correct interior wave number lies between
%          k_{exterior} [1 \sqrt{3}]; i.e. the index of refraction lies
%          between [1,3]. You may use inbuild matlab optimization
%          functions to assist you or just plain old bisection can be
%          your friend. The grid and quadrature weights for computing
%          smooth functions on boxes is provided in this file.
%          When choosing the number of points in the box, think
%          about the spacing needed in terms of the wavelength of the
%          problem. Before proceeding to part 3, you must compute
%          the correct index of refraction to 3 significant digits.
%
%  Part 3: Keeping the radii of the spheres fixed, using the
%          index of refraction you found in part 2, 
%          and maintaining at least a distance of () between them, 
%          find the configuration of spheres which minimizes the L2 norm
%          over the union of the following two boxes.
%

%% Part 1: Verifying order of convergence
clear
run ~/git/fmm3dbie/matlab/startup.m
addpath ~/git/fmm3d/matlab

close('all')

zexp = exp(2*pi*1i/6);
zuse = 1;
na = 3*ones(3,1);
S1 = geometries.ellipsoid([0.3, 0.25, 0.35], na, ...
      [real(zuse);imag(zuse);0]);

zuse = zuse*zexp;
S2 = geometries.ellipsoid([0.25, 0.35, 0.3], na, ...
      [real(zuse);imag(zuse);0]);

zuse = zuse*zexp;
S3 = geometries.ellipsoid([0.35, 0.3, 0.25], na, ...
      [real(zuse);imag(zuse);0]);

zuse = zuse*zexp;
S4 = geometries.ellipsoid([0.3, 0.25, 0.35], na, ...
      [real(zuse);imag(zuse);0]);

zuse = zuse*zexp;
S5 = geometries.ellipsoid([0.25, 0.35, 0.3], na, ...
      [real(zuse);imag(zuse);0]);

zuse = zuse*zexp;
S6 = geometries.ellipsoid([0.35, 0.3, 0.25], na, ...
      [real(zuse);imag(zuse);0]);


S = merge([S1, S2, S3, S4, S5, S6]);
plot(S);

%%

 zk = 2*pi/0.15;

 ri = 1.7;
 zks = zk*[1, sqrt(ri)];
 rep_params = ones(4,1);

 dir = [0,0,1];
 [uin, uingrad] = helm3d.planewave(zks(1),dir,S);
 
 rhs = complex(zeros(2,S.npts));
 rhs(1,:) = -uin;
 rhs(2,:) = -(uingrad(1,:).*S.n(1,:) + uingrad(2,:).*S.n(2,:) + ...
     uingrad(3,:).*S.n(3,:));
 eps = 1e-4;
 t1 = tic; [sig, errs] = solver(S, 'helm', 'trans', rhs, eps, zks, rep_params);
 tend = toc(t1);
 fprintf('Time taken in solver=%d\n',tend);

 figure
 semilogy(errs);
%%%%%%%%%%%%
%%%%%%%%%%%%
%%

umin = -5;
umax = 5;

vmin = 1.2;
vmax = 10;


xsec = [0,1,0];
xsec = xsec.'/norm(xsec);

vnulls = null(xsec.');
uaxi = vnulls(:,1);
vaxi = vnulls(:,2);

us = umin:0.05:umax;
vs = vmin:0.05:vmax;
[U,V] = meshgrid(us,vs);

X = U*uaxi(1) + V*vaxi(1);
Y = U*uaxi(2) + V*vaxi(2);
Z = U*uaxi(3) + V*vaxi(3);


sz = size(X);
xyz_out = [X(:).';Y(:).';Z(:).'];

targ_info = [];
targ_info.r = xyz_out;

t2 = tic;
eps = 1E-4;
pot = eval_fields(S, 'h', 'trans', sig, targ_info, eps, zks, rep_params); 
tend = toc(t2);
fprintf('time taken in eval=%d\n', tend)
%%

%pot_in = helm3d.kern(zk,src_info,targ_info,'s');
pot_in = helm3d.planewave(zk, dir, targ_info);
utot = pot_in + pot;
utot_plt = reshape(utot, sz);

figure;
h = surf(X,Y,Z,abs(utot_plt),'EdgeAlpha',0.1,'FaceAlpha',0.9);
colorbar;
lims = clim;

hold on

h = plot(S);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.3);
set(h,'FaceColor','none');
caxis(lims);

figure;
h = surf(X,Y,Z,real(utot_plt),'EdgeAlpha',0.1,'FaceAlpha',0.9);
colorbar;
lims = clim;

hold on

h = plot(S);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.3);
set(h,'FaceColor','none');
caxis(lims);

figure;
h = surf(X,Y,Z,imag(utot_plt),'EdgeAlpha',0.1,'FaceAlpha',0.9);
colorbar;
lims = clim;

hold on

h = plot(S);
set(h,'EdgeColor','black');
set(h,'EdgeAlpha',0.3);
set(h,'FaceColor','none');
caxis(lims);