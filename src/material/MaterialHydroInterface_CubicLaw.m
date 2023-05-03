%% MaterialInterface_Elastic class
%
% This class defines an elastic traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2023
%%%
%
%% Class definition
classdef MaterialHydroInterface_CubicLaw < MaterialHydroInterface      
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialHydroInterface_CubicLaw(parameters)
            this = this@MaterialHydroInterface('interfaceFlow', parameters, '');
            this.nStVar = 0;
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the permeability matrix
        function Te = permeabilityMtrx(this,~,~)

            % Material parameters
            w  = this.parameters(1);    % Aperture
            mu = this.parameters(2);    % Fluid dynamic viscosity
            ct = this.parameters(3);
            cb = this.parameters(4);

            % Longitudinal permeability
            kl = w*w*w/12/mu;
            
            % Permeability 
            Te = [kl, ct, cb];

        end

    end
end