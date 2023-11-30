close all; clc

nwavelengths = 50;
wavelengths = linspace(2075, 2085, nwavelengths) * 1e-9;
waveguides = ["bus", "ring"];
subs = ["b", "e"];
couplings = ["b", "e"];
Ts = zeros(2, nwavelengths);

i = 1;
for lambda = wavelengths
    clear C beta xfields yfields
    for waveguide = waveguides
        if strcmp(waveguide, 'ring')
            for sub = subs
                wvgonce
            end
        else
            wvgonce
        end
    end
    j = 1;
    for coupling = couplings
        coup
        Ts(j, i) = T;
        j = j + 1;
    end
    i = i + 1;
end

close all
figure(1)
plot(wavelengths, Ts(1, :), wavelengths, Ts(2, :))
legend({'b', 'be'})

% mode1 = 'TE';           % starting mode
% waveguide = 'ring';     % bus or ring
% sub = 'e';              % b for blood or e for ethanol (for ring only)
% concentration = 0.0008; % concentration of ethanol (for ring only)
% height = 0.5e-6;        % height of waveguide
% width = 1e-6;           % width of waveguide
% ifshow = false;         % whether show plots
% ifsave = false;         % whether save figures
% d = 1e-6;               % distance between waveguides (side to side, not center)
% textp = [-0.9, -0.2, 0.6, 0.6, 0, -0.6, -0.2, -1];
% innerR = 50e-6;         % inner radius
% R = innerR + width / 2; % average of inner R and outer R
% alpha = 10e-2;          % attenuation
% l = 14e-6;              % coupling distance
% coupling = 'e';         % bus coupled with blood only, or blood + ethanol
% 
% % materials
% nblood = 1.2966 + 0.00023294i;  % Rowe et al. 2017
% nC2H5OH = 1.3436 + 0.00017782i; % Myers et al. 2018
% nSi3N4 = 1.9810;                % Luke et al. 2015
% nSiO2 = 1.4369;                 % Malitson 1965
% nSi = 3.4699 + 4.8100e-8i;      % Shkondin et al. 2017
% nBaF2 = 1.4644;                 % Li 1980
% nCaF2 = 1.4235;                 % Li 1980
% 
% ns = nCaF2; nw = nSi3N4;
% nsname = '$\mathrm{CaF_2}$'; nwname = '$\mathrm{Si_3N_4}$';
% 
% % constants
% c = 299792458;
% epsilon0 = 8.8541878128e-12;
% mu = 4 * pi * 1e-7;
% k0 = 2 * pi / lambda;
% omega = k0 * c;
% 
% % Plot1
% mlist = [0, 1]; endsep = 0.2; nb = 2^7;
% 
% figure(1)
% set(gcf, 'Position', [100, 50, 1200, 750])
% set(0, 'defaultTextInterpreter', 'latex')
% 
% subplot(2, 3, 1)
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'YDir', 'normal')
% rectangle('Position', [-width, -height, 2 * width, 2 * height] / 2 * 1e6)
% hold on
% xline([-width, width] / 2 * 1e6, 'LineStyle', '--')
% yline([-height, height] / 2 * 1e6, 'LineStyle', '--')
% xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]')
% axis([-width, width, -height, height] * 1e6)
% text(textp(1), textp(5), nsname, 'FontSize', 18)
% text(textp(3), textp(5), nsname, 'FontSize', 18)
% text(textp(2), textp(6), nsname, 'FontSize', 18)
% text(textp(2), textp(5), nwname, 'FontSize', 18)
% if strcmp(waveguide, 'bus')
%     text(textp(2), textp(4), nsname, 'FontSize', 18)
% elseif strcmp(waveguide, 'ring')
%     if strcmp(sub, 'b')
%         text(textp(7), textp(4), '$\mathrm{Blood}$', 'FontSize', 18)
%     elseif strcmp(sub, 'e')
%         text(textp(8), textp(4), "$" + num2str((1 - concentration) * 100, '%.2f') + "\%\ \mathrm{Blood}+" +...
%             num2str(concentration * 100, '%.2f') + "\%\ \mathrm{Ethanol}$", 'FontSize', 18)
%     else
%         error('Undefined substrate. (blood or ethanol)')
%     end
% else
%     error('Undefined waveguide.')
% end
% axis equal; hold off
% 
% if strcmp(waveguide, 'bus')
%     struct_pars1 = [ns, nw, ns, height, lambda];
% elseif strcmp(waveguide, 'ring')
%     if strcmp(sub, 'b')
%         struct_pars1 = [real(nblood), nw, ns, height, lambda];
%     else
%         mixedn = concentration * real(nC2H5OH) + (1 - concentration) * real(nblood);
%         struct_pars1 = [mixedn, nw, ns, height, lambda];
%     end
% end
% 
% subplot(2, 3, 2)
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
% hold on
% for m = mlist
%     if m == 0
%         [b1, V1, blist1, Vlist1] = nee(mode1, struct_pars1, m, nb, endsep);
%     else
%         [~, V1, blist1, Vlist1] = nee(mode1, struct_pars1, m, nb, endsep);
%     end
%     plot(Vlist1, blist1, 'LineWidth', 1.2)
% end
% neff1 = sqrt(b1 * (nw^2 - ns^2) + ns^2); % field along y
% xline(V1); yline(b1)
% xlabel('$V$'); ylabel('$b$'); title('Step 1')
% legend({'$m=0$', '$m=1$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off')
% hold off
% 
% if strcmp(mode1, 'TM')
%     mode2 = 'TE';
% else
%     mode2 = 'TM';
% end
% struct_pars2 = [ns, neff1, ns, width, lambda];
% 
% subplot(2, 3 ,3)
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
% hold on
% for m = mlist
%     if m == 0
%         [b2, V2, blist2, Vlist2] = nee(mode2, struct_pars2, m, nb, endsep);
%     else
%         [~, V2, blist2, Vlist2] = nee(mode2, struct_pars2, m, nb, endsep);
%     end
%     plot(Vlist2, blist2, 'LineWidth', 1.2)
% end
% neff2 = sqrt(b2 * (neff1^2 - ns^2) + ns^2); % field along x
% xline(V2); yline(b2)
% xlabel('$V$'); ylabel('$b$'); title('Step 2')
% legend({'$m=0$', '$m=1$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off')
% hold off
% 
% beta = neff2 * k0;
% kappa_y = k0 * sqrt(nw^2 - neff1^2);
% gammaup_y = k0 * sqrt(neff1^2 - struct_pars1(1)^2);
% gammadown_y = k0 * sqrt(neff1^2 - ns^2);
% gammaleft_x = k0 * sqrt(neff2^2 - ns^2);
% gammaright_x = gammaleft_x;
% kappa_x = k0 * sqrt(neff1^2 - neff2^2);
% fprintf('neff1 = %.4f, neff2 = %.4f, beta = %.4f\n', neff1, neff2, beta)
% fprintf('kappa_y = %.4f, gammaup_y = %.4f, gammadown_y = %.4f\n', kappa_y, gammaup_y, gammadown_y)
% fprintf('kappa_x = %.4f, gammaleft_x = %.4f, gammaright_x = %.4f\n', kappa_x, gammaleft_x, gammaright_x)
% 
% % if input is TM, then y is TM, x is TE; if input is TE, then y is TE, x is
% % TM
% dx = width / 2^6; dy = height / 2^5;
% x = -2.88 * width:dx:2.88 * width; y = -12.5 / 2 * height:dy:10.5 / 2 * height;
% xseps = {x < -width / 2, x >= -width / 2 & x <= width / 2, x > width / 2};
% xs = {x(xseps{1}), x(xseps{2}), x(xseps{3})};
% yseps = {y < -height, y >= -height & y <= 0, y > 0};
% ys = {y(yseps{1}), y(yseps{2}), y(yseps{3})};
% nc = struct_pars1(1);
% shifty = y + height / 2;
% 
% if ~exist('xfields', 'var')
%     xfields = cell(3);
% end
% if ~exist('yfields', 'var')
%     yfields = cell(3);
% end
% if strcmp(waveguide, 'bus')
%     n = 1;
% elseif strcmp(waveguide, 'ring')
%     if strcmp(sub, 'b')
%         n = 2;
%     elseif strcmp(sub, 'e')
%         n = 3;
%     end
% end
% if strcmp(mode1, 'TM')
%     xfields(n, :) = {@(x) exp(gammaleft_x * (x + width / 2)), @(x) cos(kappa_x * x) / cos(kappa_x * width / 2),...
%         @(x) exp(-gammaright_x * (x - width / 2))};
%     yfields(n, :) = {@(y) 1 / ns^2 * (cos(kappa_y * height) + nw^2 * gammaup_y / (nc^2 * kappa_y) * sin(kappa_y * height)) * exp(gammadown_y * (y + height)),...
%         @(y) 1 / nw^2 * (cos(kappa_y * y) - nw^2 * gammaup_y / (nc^2 * kappa_y) * sin(kappa_y * y)),...
%         @(y) 1 / nc^2 * exp(-gammaup_y * y)};
% elseif strcmp(mode1, 'TE')
%     xfields(n, :) = {@(x) 1 / ns^2 * exp(gammaleft_x * (x + width / 2)),...
%         @(x) 1 / nw^2 * cos(kappa_x * x) / cos(kappa_x * width / 2),...
%         @(x) 1 / ns^2 * exp(-gammaright_x * (x - width / 2))};
%     yfields(n, :) = {@(y) (cos(kappa_y * height) + gammaup_y / kappa_y * sin(kappa_y * height)) * exp(gammadown_y * (y + height)),...
%         @(y) cos(kappa_y * y) - gammaup_y / kappa_y * sin(kappa_y * y),...
%         @(y) exp(-gammaup_y * y)};
% end
% xfield = [xfields{n, 1}(xs{1}), xfields{n, 2}(xs{2}), xfields{n, 3}(xs{3})];
% yfield = [yfields{n, 1}(ys{1}), yfields{n, 2}(ys{2}), yfields{n, 3}(ys{3})];
% 
% subplot(2, 3, 4)
% plot(x * 1e6, xfield, 'LineWidth', 1.2)
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
% hold on
% xline([-width, width] / 2 * 1e6, 'LineWidth', 1, 'LineStyle', '--')
% xlabel('$x$ [$\mu$m]'); ylabel('$E_x$'); title('$x$ Field')
% xlim([min(x), max(x)] * 1e6)
% grid on; hold off
% 
% subplot(2, 3, 5)
% plot(shifty * 1e6, yfield, 'LineWidth', 1.2)
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
% hold on
% xline([-height, height] / 2 * 1e6, 'LineWidth', 1, 'LineStyle', '--')
% xlabel('$y$ [$\mu$m]'); ylabel('$E_y$'); title('$y$ Field')
% xlim([min(shifty), max(shifty)] * 1e6)
% grid on; hold off
% 
% subplot(2, 3, 6)
% [Xfield, Yfield] = meshgrid(xfield, yfield);
% imagesc(x * 1e6, shifty * 1e6, Xfield .* Yfield)
% hold on
% rectangle('Position', [-width, -height, 2 * width, 2 * height] / 2 * 1e6, 'LineStyle', '--')
% set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'YDir', 'normal')
% xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]'); title("Total Field (" + mode1 + ")")
% colormap gray; axis equal; hold off
% 
% if strcmp(waveguide, 'ring')
%     sgtitle('Ring', 'FontSize', 24)
% else
%     sgtitle('Bus Waveguide', 'FontSize', 24)
% end
% 
% if ifsave
%     if strcmp(waveguide, 'ring')
%         if strcmp(sub, 'b')
%             print(figure(1), "ei_blood.png", '-dpng', '-r400');
%         elseif strcmp(sub, 'e')
%             print(figure(1), "ei_blood_ethanol.png", '-dpng', '-r400');
%         end
%     elseif strcmp(waveguide, 'bus')
%         print(figure(1), "ei_bus.png", '-dpng', '-r400');
%     end
% end
% 
% %% Normalization
% xfieldleft2 = @(x) xfields{n, 1}(x).^2;
% xfieldcore2 = @(x) xfields{n, 2}(x).^2;
% xfieldright2 = @(x) xfields{n, 3}(x).^2;
% yfielddown2 = @(y) yfields{n, 1}(y).^2;
% yfieldcore2 = @(y) yfields{n, 2}(y).^2;
% yfieldup2 = @(y) yfields{n, 3}(y).^2;
% xint = integral(xfieldleft2, -Inf, -width / 2) + integral(xfieldcore2, -width / 2, width / 2) +...
%     integral(xfieldright2, width / 2, Inf);
% yint = integral(yfielddown2, -Inf, -height) + integral(yfieldcore2, -height, 0) +...
%     integral(yfieldup2, 0, Inf);
% int = xint * yint;
% rhs = 2 * omega * mu / beta;
% if ~exist('C', 'var')
%     C = zeros(1, 3);
% end
% C(n) = sqrt(rhs / int);         % C = Ex0 * Ey0
% 
% %% coupling
% if ~all(C)
%     error('Rerun simulations for all cases before calculating K!')
% end
% cons = omega * epsilon0 * (nw^2 - ns^2) / 4;
% if ~exist('K', 'var')
%     K = zeros(1, 2);
% end
% if strcmp(coupling, 'b')
%     xpartfun = @(x) C(2) * C(1) * xfields{2, 2}(x) .* xfields{1, 1}(x);
%     ypartfun = @(y) yfields{2, 2}(y) .* yfields{1, 2}(y);
%     m = 1;
% elseif strcmp(coupling, 'e')
%     xpartfun = @(x) C(3) * C(1) * xfields{3, 2}(x) .* xfields{1, 1}(x);
%     ypartfun = @(y) yfields{3, 2}(y) .* yfields{1, 2}(y);
%     m = 2;
% else
%     error('Undefined coupling.')
% end
% xpart = integral(xpartfun, -3 / 2 * width - d, -width - d / 2);
% ypart = integral(ypartfun, -height, 0);
% K(m) = cons * xpart * ypart;
% 
% L = 2 * pi * R;
% a = exp(-alpha * L / 2);
% t = cos(K(m) * l);
% T = (t^2 + a^2 - 2 * a * t * cos(beta * L)) / (1 + a^2 * t^2 - 2 * a * t * cos(beta * L));
