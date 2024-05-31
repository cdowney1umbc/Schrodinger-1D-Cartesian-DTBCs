% The potential of a harmonic oscillator with constant k
% To abide by the restriction that the potential must vanish far from the
% origin, we implement a "breaking displacement" xbreak beyond which the
% potential is 0
function V = Harmonic_Oscillator(X,T,k,xbreak)
    V = zeros(size(X));

    % Check whether each X(i) is less than or equal to the breaking
    % displacement
    nonbreaking_X = (abs(X) <= xbreak);

    % For the X(i) which are less than or equal to the breaking
    % displacement, calculate the potential
    % The formula ensures continuity at +- xbreak
    V(nonbreaking_X) = 0.5 * k * (X(nonbreaking_X).^2 - xbreak^2);
end
