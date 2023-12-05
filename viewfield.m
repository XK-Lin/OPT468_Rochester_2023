clear; close all; clc
lambda = 2.09e-6;
ifsave = true;
wvgonce
textp = [-0.9, -0.2, 0.6, 0.6, 0, -0.6, -0.2, -1];
s = 3;

mlist = [0, 1]; endsep = 0.2; nb = 2^7;

figure(1)
set(gcf, 'Position', [100, 50, 1200, 750])
set(0, 'defaultTextInterpreter', 'latex')

subplot(2, 3, 1)
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'YDir', 'normal')
rectangle('Position', [-bus.width, -bus.height, 2 * bus.width, 2 * bus.height] / 2 * 1e6)
hold on
xline([-bus.width, bus.width] / 2 * 1e6, 'LineStyle', '--')
yline([-bus.height, bus.height] / 2 * 1e6, 'LineStyle', '--')
xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]')
axis([-bus.width, bus.width, -bus.height, bus.height] * 1e6)
text(textp(1), textp(5), nsname, 'FontSize', 18)
text(textp(3), textp(5), nsname, 'FontSize', 18)
text(textp(2), textp(6), nsname, 'FontSize', 18)
text(textp(2), textp(5), nwname, 'FontSize', 18)
if s == 1;
    text(textp(2), textp(4), nsname, 'FontSize', 18)
elseif s == 2
    text(textp(7), textp(4), '$\mathrm{Water}$', 'FontSize', 18)
elseif s ==3
    text(textp(8), textp(4), "$" + num2str((1 - concentration) * 100, '%.2f') + "\%\ \mathrm{Water}+" +...
            num2str(concentration * 100, '%.2f') + "\%\ \mathrm{Ethanol}$", 'FontSize', 18)
else
        error('Undefined substrate number. (1, 2, 3)')
end
axis equal; hold off

% change index according to the chose material
if s == 1
    struct_pars1 = [bus.ntop, bus.ncore, bus.nbottom, bus.height, lambda];
elseif s == 2
    struct_pars1 = [ringb.ntop, ringb.ncore, ringb.nbottom, ringb.height, lambda];
elseif s == 3
    struct_pars1 = [ringbe.ntop, ringbe.ncore, ringbe.nbottom, ringbe.height, lambda];
end

subplot(2, 3, 2)
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
hold on
for m = mlist
    if m == 0
        [b1, V1, blist1, Vlist1] = nee(mode, struct_pars1, m, nb, endsep);
    else
        [~, V1, blist1, Vlist1] = nee(mode, struct_pars1, m, nb, endsep);
    end
    plot(Vlist1, blist1, 'LineWidth', 1.2)
end
neff1 = sqrt(b1 * (bus.ncore^2 - bus.nbottom^2) + bus.nbottom^2); % field along y
xline(V1); yline(b1)
xlabel('$V$'); ylabel('$b$'); title('Step 1')
legend({'$m=0$', '$m=1$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off')
hold off

if strcmp(mode, 'TM')
    modee = 'TE';
else
    modee = 'TM';
end
struct_pars2 = [bus.nright, neff1, bus.nleft, bus.width, lambda];

subplot(2, 3 ,3)
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
hold on
for m = mlist
    if m == 0
        [b2, V2, blist2, Vlist2] = nee(modee, struct_pars2, m, nb, endsep);
    else
        [~, V2, blist2, Vlist2] = nee(modee, struct_pars2, m, nb, endsep);
    end
    plot(Vlist2, blist2, 'LineWidth', 1.2)
end
xline(V2); yline(b2)
xlabel('$V$'); ylabel('$b$'); title('Step 2')
legend({'$m=0$', '$m=1$'}, 'Interpreter', 'latex', 'Location', 'southeast', 'Box', 'off')
hold off

dx = bus.width / 2^6; dy = bus.height / 2^5;
x = -2.88 * bus.width:dx:2.88 * bus.width; y = -12.5 / 2 * bus.height:dy:10.5 / 2 * bus.height;
xseps = {x < -bus.width / 2, x >= -bus.width / 2 & x <= bus.width / 2, x > bus.width / 2};
xs = {x(xseps{1}), x(xseps{2}), x(xseps{3})};
yseps = {y < -bus.height, y >= -bus.height & y <= 0, y > 0};
ys = {y(yseps{1}), y(yseps{2}), y(yseps{3})};
shifty = y + bus.height / 2;
if s == 1
    xfield = [busfields.leftfield(xs{1}), busfields.xcorefield(xs{2}), busfields.rightfield(xs{3})];
    yfield = [busfields.bottomfield(ys{1}), busfields.ycorefield(ys{2}), busfields.topfield(ys{3})];
elseif s == 2
    xfield = [ringbfields.leftfield(xs{1}), ringbfields.xcorefield(xs{2}), ringbfields.rightfield(xs{3})];
    yfield = [ringbfields.bottomfield(ys{1}), ringbfields.ycorefield(ys{2}), ringbfields.topfield(ys{3})];
elseif s == 3
    xfield = [ringbefields.leftfield(xs{1}), ringbefields.xcorefield(xs{2}), ringbefields.rightfield(xs{3})];
    yfield = [ringbefields.bottomfield(ys{1}), ringbefields.ycorefield(ys{2}), ringbefields.topfield(ys{3})];
end

subplot(2, 3, 4)
plot(x * 1e6, xfield, 'LineWidth', 1.2)
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
hold on
xline([-width, width] / 2 * 1e6, 'LineWidth', 1, 'LineStyle', '--')
xlabel('$x$ [$\mu$m]'); ylabel('$E_x$'); title('$x$ Field')
xlim([min(x), max(x)] * 1e6)
grid on; hold off

subplot(2, 3, 5)
plot(shifty * 1e6, yfield, 'LineWidth', 1.2)
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
hold on
xline([-height, height] / 2 * 1e6, 'LineWidth', 1, 'LineStyle', '--')
xlabel('$y$ [$\mu$m]'); ylabel('$E_y$'); title('$y$ Field')
xlim([min(shifty), max(shifty)] * 1e6)
grid on; hold off

subplot(2, 3, 6)
[Xfield, Yfield] = meshgrid(xfield, yfield);
imagesc(x * 1e6, shifty * 1e6, Xfield .* Yfield)
hold on
rectangle('Position', [-width, -height, 2 * width, 2 * height] / 2 * 1e6, 'LineStyle', '--')
set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex', 'YDir', 'normal')
xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]'); title("Total Field (" + mode + ")")
colormap gray; axis equal; hold off

if s == 1
    sgtitle('Bus Waveguide', 'FontSize', 24)
else
    sgtitle('Ring', 'FontSize', 24)
end

if ifsave
    if s == 1
        print(figure(1), 'ei_bus.png', '-dpng', '-r400')
    elseif s == 2
        print(figure(1), 'ei_H2O.png', '-dpng', '-r400')
    elseif s ==3
        print(figure(1), 'ei_H2O_ethanol.png', '-dpng', '-r400')
    end
end