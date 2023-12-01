classdef wvgfields
    properties
        leftfield
        xcorefield
        rightfield
        topfield
        ycorefield
        bottomfield
    end
    methods
        function obj = wvgfields(wvgobj, mode)
            if strcmp(mode, 'TM')
                obj.leftfield = @(x) exp(wvgobj.gammaleft * (x + wvgobj.width / 2));
                obj.xcorefield = @(x) cos(wvgobj.kappax * x) / cos(wvgobj.kappax * wvgobj.width / 2);
                obj.rightfield = @(x) exp(-wvgobj.gammaright * (x - wvgobj.width / 2));
                obj.bottomfield = @(y) 1 / wvgobj.nbottom^2 * (cos(wvgobj.kappay * wvgobj.height) +...
                    wvgobj.ncore^2 * wvgobj.gammatop / (wvgobj.ntop^2 * wvgobj.kappay) *...
                    sin(wvgobj.kappay * wvgobj.height)) * exp(wvgobj.gammabottom * (y + wvgobj.height));
                obj.ycorefield = @(y) 1 / wvgobj.ncore^2 * (cos(wvgobj.kappay * y) - wvgobj.ncore^2 *...
                    wvgobj.gammatop / (wvgobj.ntop^2 * wvgobj.kappay) * sin(wvgobj.kappay * y));
                obj.topfield = @(y) 1 / wvgobj.ntop^2 * exp(-wvgobj.gammatop * y);
            elseif strcmp(mode, 'TE')
                obj.leftfield = @(x) 1 / wvgobj.nleft^2 * exp(wvgobj.gammaleft * (x + wvgobj.width / 2));
                obj.xcorefield = @(x) 1 / wvgobj.ncore^2 * cos(wvgobj.kappax * x) / cos(wvgobj.kappax * wvgobj.width / 2);
                obj.rightfield = @(x) 1 / wvgobj.nright^2 * exp(-wvgobj.gammaright * (x - wvgobj.width / 2));
                obj.bottomfield = @(y) (cos(wvgobj.kappay * wvgobj.height) + wvgobj.gammatop / wvgobj.kappay *...
                    sin(wvgobj.kappay * wvgobj.height)) * exp(wvgobj.gammabottom * (y + wvgobj.height));
                obj.ycorefield = @(y) cos(wvgobj.kappay * y) - wvgobj.gammatop / wvgobj.kappay * sin(wvgobj.kappay * y);
                obj.topfield = @(y) exp(-wvgobj.gammatop * y);
            end
        end
    end
end
