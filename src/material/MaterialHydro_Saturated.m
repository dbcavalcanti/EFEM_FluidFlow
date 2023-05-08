%% Material_Elastic class
%
% This class defines an linear elastic stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MaterialHydro_Saturated < MaterialHydro  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialHydro_Saturated(parameters, anm)
            this = this@MaterialHydro('saturated',parameters, anm);
            this.nStVar = 0;
        end
    end

    %% Public methods
    methods
        
        %------------------------------------------------------------------
        % Compute the permeability matrix
        function K = permeabilityMtrx(this, ~, ~, ~, ~)

            % Conductivity
            k = this.parameters(1);

            % Water specific weight
            gw = this.parameters(2);

            % Permeability matrix: K/gw = k/mu
            K = [k 0.0; 0.0 k]/gw;

        end


    end
end