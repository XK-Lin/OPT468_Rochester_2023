clear; close all; clc
nlambdas = 256;
lambdas = linspace(2090.5, 2092, nlambdas) * 1e-9;
ifsave = false;
T = zeros(nlambdas, 2);
i = 1;
for lambda = lambdas
    wvgonce
    L = 2 * pi * R;
    a = exp(-alpha * L / 2);
    t = cos(K * l);
    T(i, :) = (t.^2 + a^2 - 2 * a * t .* cos(beta * L)) ./ (1 + a^2 * t.^2 - 2 * a * t .* cos(beta * L));
    i = i + 1;
end
%%
close all
figure(1)
set(gcf, 'Position', [100, 50, 700, 500])
set(0, 'defaultTextInterpreter', 'latex')

plot(lambdas * 1e9, T(:, 1), lambdas * 1e9, T(:, 2), 'LineWidth', 1.2)
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
xlabel('Wavelength [nm]'); ylabel('Transmittance $T$'); title('Ring Resonator Transmission Spectrum')
legend({'$\mathrm{H_2O}$', '$\mathrm{H_2O}$ with Ethanol'}, 'FontSize', 20, 'Box', 'off', 'Interpreter', 'latex', 'Location', 'southeast')
xlim([lambdas(1), lambdas(end)] * 1e9); grid on
if ifsave
    print(figure(1), 'Transmission.png', '-dpng', '-r400')
end