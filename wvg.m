classdef wvg
    properties
        nleft {mustBeNonnegative}
        nright {mustBeNonnegative}
        ncore {mustBeNonnegative}
        ntop {mustBeNonnegative}
        nbottom {mustBeNonnegative}
        width {mustBeNonnegative}
        height {mustBeNonnegative}
        gammaleft {mustBeNonnegative}
        kappax {mustBeNonnegative}
        gammaright {mustBeNonnegative}
        gammatop {mustBeNonnegative}
        kappay {mustBeNonnegative}
        gammabottom {mustBeNonnegative}
        neff1 {mustBeNonnegative}
        neff2 {mustBeNonnegative}
    end
    methods
        function obj = wvg(nl, nc, nr, nt, nb, w, h, n1, n2, kghcws, kgvcws)
            obj.nleft = nl;
            obj.ncore = nc;
            obj.nright = nr;
            obj.ntop = nt;
            obj.nbottom = nb;
            obj.width = w;
            obj.height = h;
            obj.gammaleft = kghcws(3);
            obj.kappax = kghcws(2);
            obj.gammaright = kghcws(1);
            obj.gammatop = kgvcws(1);
            obj.kappay = kgvcws(2);
            obj.gammabottom = kgvcws(3);
            obj.neff1 = n1;
            obj.neff2 = n2;
        end
    end
end
