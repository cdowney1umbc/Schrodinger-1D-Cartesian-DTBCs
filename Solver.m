% This script calculates the values of Psi_j^n
% To create an animation using those values, run the script Animator.m

% This program uses a tridiagonal matrix solver written by Tamas Kis
% Tamas Kis (2024). Tridiagonal Matrix Algorithm
% (https://github.com/tamaskis/tridiagonal-MATLAB/releases/tag/v6.0.1),
% GitHub. Retrieved May 30, 2024.
import tridiagonal_MATLAB_main.*;

% Physical constants
hbar = 1;
m = 1;

%% Step 1
% Continuous bounds
L = 10;
T = 60;

% Discrete bounds
J = 100;
N = 600;

% Step sizes
h_x = L/J;
h_t = T/N;

% Spatial and temporal arrays
% x_j is 2J+1 by 1 with x_j(j) = x_{j-1-J}
% And t_n is 1 by N+1 with t_n(n) = t^{n-1}
x_j = linspace(-L, L, (2 * J) + 1).';
t_n = linspace(0, T, N + 1);

% Spacetime grid
% x_jn and t_jn are 2J+1 by N+1
% With (x_jn(j,n), t_jn(j,n)) = (x_{j-1-J}, t^{n-1})
x_jn = repmat(x_j, 1, N + 1);
t_jn = repmat(t_n, (2 * J) + 1, 1);

%% Step 2
% Potential function, selected from the '+Potential_Functions' subdirectory
import Potential_Functions.*;
Potential_Function = @(X,T) Harmonic_Oscillator(X,T,0.5,2.5);

% Initial condition, selected from the '+Initial_Conditions' subdirectory
import Initial_Conditions.*;
Initial_Condition = @(X) Wave_Packet_Min_Sigma(X, hbar, -6, 10, 1);

% Potential array, initialized using the selected potential
% V_jn is 2J+1 by N+1 with V_jn(j,n) = V(x_{j-1-J}, t^{n-1} + h_t/2)
V_jn = Potential_Function(x_jn, t_jn + (0.5 * h_t));
% We enforce the support restriction even if the potential function
% does not obey it
V_jn([1,2,2*J,(2*J)+1],:) = 0;

% Wave function / solution
% Psi_jn is 2J+1 by N+1 with Psi_jn(j,n) = Psi_{j-1-J}^{n-1}
Psi_jn = zeros((2 * J) + 1, N + 1);

% Initialize t=0 values of Psi_jn using the selected initial condition
Psi_jn(:,1) = Initial_Condition(x_j);
% We enforce the support restriction even if the initial condition
% does not obey it
Psi_jn([1,2,2*J,(2*J)+1],:) = 0;

%% Step 3
% Shorthand constants
c_V = (2 * m * h_x^2) / hbar^2;
c_t = (4 * m * h_x^2) / (hbar * h_t);
phi = atan(4 / c_t);
mu = c_t / sqrt(c_t^2 + 16);
chi0 = 1 - (1i * c_t / 2) - (1i * exp(1i * phi / 2) * c_t / (2 * sqrt(mu)));

%% Step 4
% Calculate P_n(mu) for n in [0,N]
% P_n_mu is 1 by N+1 with P_n_mu(n) = P_{n-1}(mu)
P_n_mu = ones(1, N + 1);
P_n_mu(2) = mu;
for n = 2:N
    P_n_mu(n+1) = ((((2 * n) - 1) * mu * P_n_mu(n)) ...
        - (n - 1) * P_n_mu(n-1)) / n;
end

% Calculate chi_n for n in [1,N]
% chi_n is 1 by N with chi_n(n) = chi_n
% This prioritizes speed over readability. You are not required to
% understand this
chi_n = (-1i * c_t * (-1).^(1:N)) ...
    - (1i * exp(-1i * phi / 2) * c_t / (2 * sqrt(mu))) ...
    * ((exp(-1i * phi * ((1:N) - 1)) .* P_n_mu(2:N+1)) ...
    + (exp(-1i * phi * (1:N)) .* P_n_mu(1:N)) ...
    + 4 * (-1).^(1:N) * mu .* cumsum((-1).^(0:N-1) .* exp(-1i * phi * (0:N-1)) .* P_n_mu(1:N)));

%% Step 5
% Calculate all values that stay the same at all time steps
subdiag = ones(2 * J, 1);
diag_const = [-chi0; (-2 + (1i * c_t)) * ones((2 * J) - 1, 1); -chi0];
superdiag = ones(2 * J, 1);

% For each time step but the last
for n = 1:N
    diag = diag_const - c_V * [0; V_jn(2:2*J, n); 0];

    d_n = [sum(Psi_jn(1, 1:n) .* chi_n(n:-1:1), "all"); ...
        - Psi_jn(1:(2*J)-1, n) ...
        + ((2 + (1i * c_t)) * Psi_jn(2:2*J, n)) ...
        + (c_V * V_jn(2:2*J, n) .* Psi_jn(2:2*J, n)) ...
        - Psi_jn(3:(2*J)+1, n); ...
        sum(Psi_jn((2*J)+1, 1:n) .* chi_n(n:-1:1), "all")];

    % Solve for the values of the wave function at the next time step
    Psi_jn(:, n+1) = tridiagonal_vector(subdiag, diag, superdiag, d_n);
    % Tamas Kis (2024). Tridiagonal Matrix Algorithm
    % (https://github.com/tamaskis/tridiagonal-MATLAB/releases/tag/v6.0.1),
    % GitHub. Retrieved May 30, 2024.
end

%% Post processing
% The real and imaginary parts of the wave function
RePsi_jn = real(Psi_jn);
ImPsi_jn = imag(Psi_jn);

% The magnitude of the wave function
MagPsi_jn = abs(Psi_jn);

% The maximum value that the magnitude takes
MagPsi_max = max(MagPsi_jn, [], "all");

% When the solution is displayed we want the potential to be displayed as
% well, so we calculate the potential at every grid point
% Valigned_jn is 2J+1 by N+1
% with Valigned_jn(j,n) = V(x_{j-1-J}, t^{n-1})
Valigned_jn = Potential_Function(x_jn, t_jn);
% We enforce the support restriction even if the potential function
% does not obey it
Valigned_jn([1,2,2*J,(2*J)+1],:) = 0;

% The minimum and maximum potentials
V_max = max(Valigned_jn, [], "all");
V_min = min(Valigned_jn, [], "all");
