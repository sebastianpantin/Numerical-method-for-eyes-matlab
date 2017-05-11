classdef Lenses
%   Creates objects for the cornea and lens.
%   Also calculates focal length.
%   Can be used with values for radius of curvature or optical power

% Properties:
%   P1          - Frontside optical power
%   P2          - Backside optical power
%   RefIndx1    - Refractive index prior to the lens
%   RefIndx2    - Refractive index in the lens
%   d           - Thickness of the lens
%   focal       - Focal length of the lens


% Methods:
%   this_lens   - Defines the properties of the lens
%   get_power   - Calculates the power of the lens


    properties
%        FrontRad       % Frontside radius of curvature
%        BackRad        % Backside radius of curvature
        P1
        P2
        RefIndx1
        RefIndx2
        d
        focal
        power
        TF
    end

    methods
        function this_lens = Lenses(P1,P2,RI1,RI2,d)
            this_lens.P1 = P1;
            this_lens.P2 = P2;
            this_lens.RefIndx1 = RI1;
            this_lens.RefIndx2 = RI2;
            this_lens.d = d;
        end
        
%        function f = get_P(obj)
%            n=obj.RefIndx2/obj.RefIndx1;
%            f=(1/((n-1)*(1/obj.FrontRad-1/obj.BackRad+((n-1)*obj.d)/(n*obj.FrontRad*obj.BackRad))));
%            %f=(1/((obj.RefIndx1/obj.RefIndx2-1)*(1/obj.FrontRad+1/obj.BackRad)));
%        end

        function power = get_power(obj)
            power = obj.P1+obj.P2-obj.P1*obj.P2*obj.d/obj.RefIndx2;
        end        
    end
end
