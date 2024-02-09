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
nref_max = ;

%  Set order of discretization, must be less than 10
norder = ;


%  Load points in the interior and exterior of stellarator
nin = ;
[xyz_in] = get_interior_points_stellarator(nin);
strengths = rand(nin,1) + 1j*rand(nin,1);

nout = ;
[xyz_out] = get_exterior_points_stellarator(nout);


%  Compute exact potential at the exterior targets
potex = ;

%  Create a geometry with base number of patches, and plot
%  the geometry and interior and exterior points



%  For each refinement parameter, load 
for i = 1:nref_max

%  Set number of patches


%  Load the geometry with given number of patches


%  Generate boundary data for given discretization


%  Solve the boundary value problem


%  Evaluate the potential at the target locations


%  Compute and store the error for given refinement

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
