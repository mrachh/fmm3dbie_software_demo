%%  Learning Goals:
%  
% In this demo, you will learn how to manipulate three dimensional
% surfaces, defined via analytic parametrizations.
%
% You will also learn how to use 
%  Load simple geometry (step 3)
%
%
R = 1;
S = geometries.sphere(R, 4);

% Plot the domain
figure
plot(S);

% plot the wireframe

% scatter plot
figure
scatter(S);
axis equal


% Look at the normals
% quiver(S);

% Look inside S (don't change alpha, really look inside, search your
% feelings)


%% Compute the surface area, centroid, and moment of inertia
a = sum(S.wts(:));
erra = abs(a - 4*pi*R^2);


%% User exercise: 
% Use this to compute the centroid, and moment of interia 


%% Harder integral (Gauss' Law)
% \int_{\Gamma} \nabla_{x} \frac{1}{4 \pi |x-y|} \cdot n(x) dx = \delta_{S}

xyz = -[0.3; 0.2; 0.1];
dx = S.r(1,:) - xyz(1);
dy = S.r(2,:) - xyz(2);
dz = S.r(3,:) - xyz(3);
r = sqrt(dx.^2 + dy.^2 + dz.^2);

rdotn = dx.*S.n(1,:) + dy.*S.n(2,:) + dz.*S.n(3,:);
fint = (rdotn./r.^3/4/pi);

% Plot surface with function
plot(S, fint)

% Show quality of function
% ferr = estimate_error(S, fint);
% plot(S, ferr);

err = (fint*S.wts - 1);
fprintf('Error at easy point = %d\n', err);

%% Harder test, using anonymous functions
% Move point closer to boundary -> make a function

xyz = -[0.95; 0.01; 0.03];
dx = @(S, xyz) S.r(1,:) - xyz(1);
dy = @(S, xyz) S.r(2,:) - xyz(2);
dz = @(S, xyz) S.r(3,:) - xyz(3);
r = @(S, xyz) sqrt(dx(S, xyz).^2 + dy(S, xyz).^2 + dz(S, xyz).^2);

rdotn = @(S, xyz) dx(S, xyz).*S.n(1,:) + ...
                   dy(S, xyz).*S.n(2,:) + dz(S, xyz).*S.n(3,:);
fint = @(S, xyz) (rdotn(S, xyz)./r(S, xyz).^3/4/pi);

fint1 = fint(S, xyz);

plot(S, fint1);
err = (fint1*S.wts - 1);
fprintf('Error at difficult point = %d\n', err);

% 
S2 = oversample(S, 20);

fint2 = fint(S2, xyz);
err = (fint2*S2.wts - 1);
fprintf('Error at difficult point post oversampling= %d\n', err);

%% User exercise:
abc = [2;3.1; 1.7];
S_ellip = geometries.ellipsoid(abc);
plot(S_ellip)


%% Repeat above experiments to verify area using mathematica/wikipedia or
% elliptic integrals if you are brave, and the repeat the gauss test

% Rotate

% translate

% scale


%% Surface smoother
% Surface smoother (supports many different low order meshes)

%% User exercise: Play around with smoothing parameter, and a different
% geometry, plot the errors in refined surfaces, see convergence


