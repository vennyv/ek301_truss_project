%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input C matrix
connectivity_matrix = [1 1 0 0 0 0 0 0 0 0 0 0 0;
                       1 0 1 1 0 0 0 0 0 0 0 0 0;
                       0 1 1 0 1 1 0 1 0 0 0 0 0;
                       0 0 0 1 1 0 1 0 0 0 0 0 0;
                       0 0 0 0 0 1 1 0 1 1 0 0 0;
                       0 0 0 0 0 0 0 1 1 0 1 1 0;
                       0 0 0 0 0 0 0 0 0 1 1 0 1;
                       0 0 0 0 0 0 0 0 0 0 0 1 1];

sizes = size(connectivity_matrix);
joints = sizes(1,1);
members = sizes(1,2);

% Input Sx and Sy matrices
Sx = zeros(8,3);
Sx(1,1) = 1;

Sy = zeros(8,3);
Sy(1,2) = 1;
Sy(8,3) = 1;

% Input X and Y vectors
X = [0; 6.25; 12.5; 12.5; 19.5; 19.5; 25.25; 31];
Y = [0; 6; 0; 12; 12; 0; 6; 0];

% Input L vector (external loads)
L = zeros(16,1);
L(11) = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up A matrix
A = zeros(joints*2, members+3);
A(1:joints, 1:members) = connectivity_matrix;
A(joints+1:end, 1:members) = connectivity_matrix;

[a, b] = find(connectivity_matrix == 1);

% Calculate total distance and lengths of members
distvec = zeros(members, 1);
for x = 1:members
    member_indices = find(connectivity_matrix(:,x));
    xdist = X(member_indices(2)) - X(member_indices(1));
    ydist = Y(member_indices(2)) - Y(member_indices(1));
    distvec(x) = sqrt(xdist^2 + ydist^2);
end

tdist = sum(distvec);

% Calculate total cost
cost = 10 * joints + tdist;

% Calculating A matrix
for x = 1:members
    member_indices = find(connectivity_matrix(:,x));
    xdist = X(member_indices(2)) - X(member_indices(1));
    ydist = Y(member_indices(2)) - Y(member_indices(1));
    dist = sqrt(xdist^2 + ydist^2);
    A(member_indices(1), x) = xdist/dist;
    A(member_indices(2), x) = -xdist/dist;
    A(member_indices(1)+joints, x) = ydist/dist;
    A(member_indices(2)+joints, x) = -ydist/dist;
end

A(1:joints, members+1:end) = Sx;
A(joints+1:end, members+1:end) = Sy;

% Calculate T vector using backslash operator for stability
T = A \ L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print Results
fprintf('EK301, Section A2, Group 20: Darren S., Venessa M., Vikram B. 4/2/2024.\n');
fprintf('Load: %.1f oz\n', sum(L));
fprintf('Member forces in oz:\n');
for z = 1:members
    if T(z) > 0
        fprintf('m%d: %.3f (T)\n', z, T(z));
    elseif T(z) < 0
        fprintf('m%d: %.3f (C)\n', z, -T(z));
    else
        fprintf('m%d: %.3f\n', z, T(z));
    end
end

fprintf('Reaction forces in oz:\n');
fprintf('Sx1: %.2f\n', T(members+1));
fprintf('Sy1: %.2f\n', T(members+2));
fprintf('Sy2: %.2f\n', T(members+3));

fprintf('Cost of truss: $%.2f\n',cost)
fprintf('Theoretical max load/cost ratio in oz/$: 0.241\n')

% Constants for the fit formula
fit_coefficient = 3055; % oz/in
alpha = 2.009;
U_fit = 1.36; % oz

% Calculate nominal buckling load for each member
Pnom = fit_coefficient * distvec.^(-alpha);

% Initialize vectors for failure loads
failure_load_nominal = zeros(members, 1);
failure_load_strong = zeros(members, 1);
failure_load_weak = zeros(members, 1);

% Calculate failure loads for compressive forces
for i = 1:members
    if T(i) < 0 % Compressive force
        failure_load_nominal(i) = Pnom(i);
        failure_load_strong(i) = Pnom(i) + U_fit;
        failure_load_weak(i) = Pnom(i) - U_fit;
    end
end

% Find the minimum failure load across all compressive members for each case
W_failure_nominal = min(failure_load_nominal(failure_load_nominal > 0));
W_failure_strong = min(failure_load_strong(failure_load_strong > 0));
W_failure_weak = min(failure_load_weak(failure_load_weak > 0));

% Estimate the uncertainty in the failure load
Delta_W = (W_failure_strong - W_failure_weak) / 2;
W_failure = W_failure_nominal; % Your predicted failure load is the nominal one

% Printing the results for buckling analysis
fprintf('Failure load (Nominal): %.2f oz\n', W_failure_nominal);
fprintf('Failure load (Strong): %.2f oz\n', W_failure_strong);
fprintf('Failure load (Weak): %.2f oz\n', W_failure_weak);
fprintf('Estimated failure load: %.2f oz ± %.2f oz\n', W_failure, Delta_W);