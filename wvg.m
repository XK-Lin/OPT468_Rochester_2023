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
