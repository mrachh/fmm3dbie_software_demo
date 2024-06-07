%
%  Write a code that illustrates the order of convergence
%  of the solution to the Helmholtz dirichlet problem
%  in the exterior of the stellarator and compare the
%  it to the expected $O(h^{p-1})
%  
%
%  todo: convert sphere code to stellarator later

%  Base number of quad patches in each direction
nu_base = 5;
nv_base = 15;

%  Set max number of refinements (more than 3 might take too long) 
nref_max = 3;

%  Set order of discretization, must be less than 10
norder  = 4;

nu = 1; 
iptype = 11;

S = surfer.sphere(norder, nu, iref, iptype);

%% Plot the geometry
%  Create a geometry with base number of patches, and plot
%  the geometry and interior and exterior points

plot(S);

%% Identify and test points in the interior and exterior

%% Load points in the interior and exterior of stellarator
nin = 5;
[xyz_in] = get_interior_points_stellarator(nin);
strengths = rand(nin,1) + 1j*rand(nin,1);

nout = 5;
[xyz_out] = get_exterior_points_stellarator(nout);

plot(S);
hold on;

scatter3(xyz_in(1,:), xyz_in(2,:), xyz_in(3,:), 'r.');
scatter3(xyz_out(1,:), xyz_out(2,:), xyz_out(3,:), 'b.');
%  Compute exact potential at the exterior targets

srcinfo = [];
targinfo = [];
srcinfo.r = xyz_in;
targinfo.r = xyz_out; 
potex = helm3d.kern(zk, srcinfo, targinfo, 's')*strengths;

errs= zeros(size(0:nref_max));

tol = 1e-8;
zpars = [zk, -1j*zk ,1];
%%  For each refinement parameter, load 
for iref = 0:nref_max
    
%  Set number of patches
    


%  Load the geometry with given number of patches
    S = surfer.sphere(norder, nu, iref, iptype);
    wts = S.wts;
    wts = wts(:);


%  Generate boundary data for given discretization
    rhs = helm3d.kern(zk, srcinfo, S, 's')*strengths;
   

%  Solve the boundary value problem
    sig = helm3d.dirichlet.solver(S, zpars, rhs, tol);

%  Evaluate the potential at the target locations
    dat = helm3d.kern(zk,S,targinfo,'c',zpars(2),zpars(3));
    pot = dat*(sig.*wts);


%  Compute and store the error for given refinement
    errs(iref+1) = norm(pot-potex);

end

%  Plot the errors and verify order of convergence


function xyz_in = get_interior_points_stellarator(n)
    thet = rand(n,1)*pi;
    phi = rand(n,1)*2*pi;
    rr = rand(n,1)*0.7;
    
    xyz_in = zeros(3,n);
    xyz_in(1,:) = rr.*sin(thet).*cos(phi);
    xyz_in(2,:) = rr.*sin(thet).*sin(phi);
    xyz_in(3,:) = rr.*cos(thet);
    

end


function xyz_out = get_exterior_points_stellarator(n)
    thet = rand(n,1)*pi;
    phi = rand(n,1)*2*pi;
    rr = 1.2 + rand(n,1)*0.7;
    
    xyz_out = zeros(3,n);
    xyz_out(1,:) = rr.*sin(thet).*cos(phi);
    xyz_out(2,:) = rr.*sin(thet).*sin(phi);
    xyz_out(3,:) = rr.*cos(thet);
end
