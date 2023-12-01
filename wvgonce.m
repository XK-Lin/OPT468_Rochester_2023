close all; clc
defpars

% initialize C and K
C = zeros(1, 3); K = zeros(1, 2); beta = zeros(1, 2);

% effective index
[neff1, neff2, kghcws, kgvcws] = effi(ns, nw, ns, ns, ns, width, height, lambda, mode);
bus = wvg(ns, nw, ns, ns, ns, width, height, neff1, neff2, kghcws, kgvcws);
[neff1, neff2, kghcws, kgvcws] = effi(ns, nw, ns, nb, ns, width, height, lambda, mode);
ringb = wvg(ns, nw, ns, nb, ns, width, height, neff1, neff2, kghcws, kgvcws);
[neff1, neff2, kghcws, kgvcws] = effi(ns, nw, ns, nbe, ns, width, height, lambda, mode);
ringbe = wvg(ns, nw, ns, nbe, ns, width, height, neff1, neff2, kghcws, kgvcws);

% find fields of all waveguides
busfields = wvgfields(bus, mode);
ringbfields = wvgfields(ringb, mode);
ringbefields = wvgfields(ringbe, mode);

% normalize fields
C(1) = getC(bus, busfields, omega, mu, k0);
C(2) = getC(ringb, ringbfields, omega, mu, k0);
C(3) = getC(ringbe, ringbefields, omega, mu, k0);

% coupling coefficient
[K(1), beta(1)] = getK(ringb, ringbe, busfields, ringbfields, ringbefields, C, d, k0, omega, epsilon0, 'b');
[K(2), beta(2)] = getK(ringb, ringbe, busfields, ringbfields, ringbefields, C, d, k0, omega, epsilon0, 'be');

% coupling coefficient

function C = getC(wvgobj, wvgfieldobj, omega, mu, k0)
    xfields2 = {@(x) wvgfieldobj.leftfield(x).^2, @(x) wvgfieldobj.xcorefield(x).^2, @(x) wvgfieldobj.rightfield(x).^2};
    yfields2 = {@(y) wvgfieldobj.bottomfield(y).^2, @(y) wvgfieldobj.ycorefield(y).^2, @(y) wvgfieldobj.topfield(y).^2};
    xint = integral(xfields2{1}, -Inf, -wvgobj.width / 2) + integral(xfields2{2}, -wvgobj.width / 2, wvgobj.width / 2) +...
        integral(xfields2{3}, wvgobj.width / 2, Inf);
    yint = integral(yfields2{1}, -Inf, -wvgobj.height) + integral(yfields2{2}, -wvgobj.height, 0) +...
        integral(yfields2{3}, 0, Inf);
    int = xint * yint; rhs = 2 * omega * mu / (wvgobj.neff2 * k0);
    C = sqrt(rhs / int);
end

function [K, beta] = getK(ringbobj, ringbeobj, busfieldobj, ringbfieldobj, ringbefieldobj, C, d, k0, omega, epsilon0, coup)
    cons = omega * epsilon0 * (ringbobj.ncore^2 - ringbobj.nbottom^2) / 4;
    if strcmp(coup, 'b')
        xpartfun = @(x) C(2) * C(1) * busfieldobj.leftfield(x) .* ringbfieldobj.xcorefield(x);
        ypartfun = @(y) busfieldobj.ycorefield(y) .* ringbfieldobj.ycorefield(y);
        beta = ringbobj.neff2 * k0;
    elseif strcmp(coup, 'be')
        xpartfun = @(x) C(3) * C(1) * busfieldobj.leftfield(x) .* ringbefieldobj.xcorefield(x);
        ypartfun = @(y) busfieldobj.ycorefield(y) .* ringbefieldobj.ycorefield(y);
        beta = ringbeobj.neff2 * k0;
    else
        error('Undefined coupling')
    end
    xpart = integral(xpartfun, -3 / 2 * ringbobj.width - d, -ringbobj.width - d / 2);
    ypart = integral(ypartfun, -ringbobj.height, 0);
    K = cons * xpart * ypart;
end
