% A wave packet with a given mean position, mean momentum,
% and standard deviation of position,
% and with the minimum standard deviation of momentum allowed
% by the Heisenberg uncertainty principle
function F = Wave_Packet_Min_Sigma(X, hbar, mu_x, mu_p, sigma_x)
    F = exp(- (X - mu_x).^2 / (4 * sigma_x^2)) ...
        .* exp(1i * mu_p * X / hbar) / sqrt(sqrt(2*pi) * sigma_x);
end
