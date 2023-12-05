mode = 'TE';            % starting mode
concentration = 0.05;   % concentration of ethanol (for ring only)
height = 0.5e-6;        % height of waveguide
width = 1e-6;           % width of waveguide
innerR = 15e-6;         % inner radius
R = innerR + width / 2; % average of inner R and outer R
alpha = 100;            % attenuation
l = 2e-6;               % coupling distance
d = 0.5e-6;               % distance between waveguides (side to side, not center)

% materials
nblood = 1.2966 + 0.00023294i;  % Rowe et al. 2017
nC2H5OH = 1.3436 + 0.00017782i; % Myers et al. 2018
nSi3N4 = 1.9810;                % Luke et al. 2015
nSiO2 = 1.4369;                 % Malitson 1965
nSi = 3.4699 + 4.8100e-8i;      % Shkondin et al. 2017
nBaF2 = 1.4644;                 % Li 1980
nCaF2 = 1.4235;                 % Li 1980
nH2O = 1.3020;                  % Hale and Querry, 1973

% waveguide materials
ns = nCaF2; nw = nSi3N4;
nsname = '$\mathrm{CaF_2}$'; nwname = '$\mathrm{Si_3N_4}$';
nb = real(nH2O);
nbe = concentration * real(nC2H5OH) + (1 - concentration) * nb;

% constants
c = 299792458;
epsilon0 = 8.8541878128e-12;
mu = 4 * pi * 1e-7;
k0 = 2 * pi / lambda;
omega = k0 * c;
