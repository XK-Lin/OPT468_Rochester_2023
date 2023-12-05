close all; clc
defpars
lambda = 2090e-9;
wvgonce
L = 2 * pi * R;
a = exp(-alpha * L / 2);
t = cos(K * l);
K
FSR = [lambda^2 ./ (ringb.neff2 * L), lambda^2 ./ (ringbe.neff2 * L)] * 1e9
FWHM = [lambda^2 * (1 - a * t(1)) ./ (pi * ringb.neff2 * L * sqrt(a * t(1))),...
    lambda^2 * (1 - a * t(2)) ./ (pi * ringbe.neff2 * L * sqrt(a * t(2)))] * 1e9
Q = [pi * ringb.neff2 * L * sqrt(a * t(1)) / (lambda * (1 - a * t(1))),...
    pi * ringbe.neff2 * L * sqrt(a * t(2)) / (lambda * (1 - a * t(2)))]