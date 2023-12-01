function [b, V, blist, Vlist] = nee(mode, struct_pars, m, nb, endsep)
%
% Calculate normalized eigenvalue equation for slab waveguide.
%
% Inputs are mode, struct_pars, m, nb, endsep. nb is the array size of b.
%
% struct_pars is an array [nc, nw, ns, h, lambda].
% 
% endsep is how close blist should approach 1 (1e-2 would be fine).
%
% Outputs are the crossing's b and V, and arrays for V and b for plot.
%

% check whether the input has the correct number of entries
if length(struct_pars) ~= 5
    error('Please give an array [nc, nw, ns, h, lambda] as the structure input.')
end

% define indices, waveguide size, and wavelength
nc = struct_pars(1); nw = struct_pars(2); ns = struct_pars(3);
h = struct_pars(4); lambda = struct_pars(5);

% define free-space wavevector
k0 = 2 * pi / lambda;

% calculate normalized parameter a and V
a = (ns^2 - nc^2) / (nw^2 - ns^2);
V = h * k0 * sqrt(nw^2 - ns^2);

% solve the normalized-parameter eigenvalue equation
syms bs

% for TE mode
if strcmp(mode, 'TE')
    sol = vpasolve(V == 1 ./ sqrt(1 - bs) .* (atan(sqrt(bs ./ (1 - bs))) ...
        + atan(sqrt((a + bs) ./ (1 - bs))) + m * pi), bs, 0.5);
% for TM mode
elseif strcmp(mode, 'TM')
    sol = vpasolve(V == 1 ./ sqrt(1 - bs) .* (atan(nw^2 / ns^2 * ...
        sqrt(bs ./ (1 - bs))) + atan(nw^2 / nc^2 * sqrt((a + bs) ./ ...
        (1 - bs))) + m * pi), bs, 0.5);
% display error if the mode does not match TE or TM
else
    error('Undefined mode.')
end

% convert MATLAB symbol to double
b = double(sol);

% define a list of b for plot
blist = linspace(1e-6, 1 - endsep, nb);

% for TE mode, calculate a list of V
if strcmp(mode, 'TE')
    Vlist = 1 ./ sqrt(1 - blist) .* (atan(sqrt(blist ./ (1 - blist))) +...
        atan(sqrt((a + blist) ./ (1 - blist))) + m * pi);
% for TM mode, calculate a list of V
elseif strcmp(mode, 'TM')
    Vlist = 1 ./ sqrt(1 - blist) .* (atan(nw^2 / ns^2 * sqrt(blist ./ (1 - blist))) + ...
        atan(nw^2 / nc^2 * sqrt((a + blist) ./ (1 - blist))) + m * pi);
end

