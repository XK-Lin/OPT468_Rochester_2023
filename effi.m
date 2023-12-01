function [neff1, neff2, kghcws, kgvcws] = effi(nl, nc, nr, nt, nb, w, h, lambda, mode)
%
% Implement the effective index method to solve for the
% fundamental mode of a rectangular waveguide.
%
% Mode is either 'TE' or 'TM'.
%
% isleft and isright asks whether we need nee for the left and right
% region. If they are false, then the ori
%
% The outputs are neff1, neff2, and the final beta.
%

n1left = nl; n1right = nr;

% for TE mode
if strcmp(mode, 'TE')
    % Step 1
    % calculations for the core region
    b1center = nee('TE', [nt, nc, nb, h, lambda], 0, 1, 1e-4);
    n1center = sqrt(b1center * (nc^2 - nb^2) + nb^2);

    % Step 2
    b = nee('TM', [n1right, n1center, n1left, w, lambda], 0, 1, 1e-4);
    n = sqrt(b * (n1center^2 - n1left^2) + n1left^2);

% for TM mode
elseif strcmp(mode, 'TM')
    % Step 1
    % calculations for the core region
    b1center = nee('TM', [nt, nc, nb, h, lambda], 0, 1, 1e-4);
    n1center = sqrt(b1center * (nc^2 - nb^2) + nb^2);

    % Step 2
    b = nee('TE', [n1right, n1center, n1left, w, lambda], 0, 1, 1e-4);
    n = sqrt(b * (n1center^2 - n1left^2) + n1left^2);
else
    error('What is the polarization? (horizontal or vertical)')
end

% define neff1 and neff2 (the outputs of the function)
neff1 = n1center;
neff2 = n;

% calculate the wavevector
k0 = 2 * pi / lambda;

% using the expressions we derived during class, calculate gamma and kappa
gammacy = k0 * sqrt(neff1^2 - nt^2);
kappay = k0 * sqrt(nc^2 - neff1^2);
gammasy = k0 * sqrt(neff1^2 - nb^2);
gammacx = k0 * sqrt(neff2^2 - n1right^2);
kappax = k0 * sqrt(neff1^2 - neff2^2);
gammasx = k0 * sqrt(neff2^2 - n1left^2);

% define output for vertical direction (kgv) and horizontal direction (kgh)
kgvcws = [gammacy, kappay, gammasy];
kghcws = [gammacx, kappax, gammasx];

% warn the user about the definition of cladding and substrate
% fprintf('Warning: c is up and right, s is down and left (influences gammac and gammas)!\n')

