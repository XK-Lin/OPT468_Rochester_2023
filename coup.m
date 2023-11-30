if ~all(C)
    error('Rerun simulations for all cases before calculating K!')
end
cons = omega * epsilon0 * (nw^2 - ns^2) / 4;
if ~exist('K', 'var')
    K = zeros(1, 2);
end
if strcmp(coupling, 'b')
    xpartfun = @(x) C(2) * C(1) * xfields{2, 2}(x) .* xfields{1, 1}(x);
    ypartfun = @(y) yfields{2, 2}(y) .* yfields{1, 2}(y);
    betaa = beta(2);
    m = 1;
elseif strcmp(coupling, 'e')
    xpartfun = @(x) C(3) * C(1) * xfields{3, 2}(x) .* xfields{1, 1}(x);
    ypartfun = @(y) yfields{3, 2}(y) .* yfields{1, 2}(y);
    betaa = beta(3);
    m = 2;
else
    error('Undefined coupling.')
end
xpart = integral(xpartfun, -3 / 2 * width - d, -width - d / 2);
ypart = integral(ypartfun, -height, 0);
K(m) = cons * xpart * ypart;

L = 2 * pi * R;
a = exp(-alpha * L / 2);
t = cos(K(m) * l);
T = (t^2 + a^2 - 2 * a * t * cos(betaa * L)) / (1 + a^2 * t^2 - 2 * a * t * cos(betaa * L));
