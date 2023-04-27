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
classdef MaterialInterface_CubicLaw < MaterialInterface      
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_CubicLaw(parameters)
            this = this@MaterialInterface('interfaceFlow', parameters, '');
            this.nStVar = 0;
        end
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the elastic traction vector
        %
        % Input:
        %   w0: current jump displacement vector in the local system [w0s, w0n]
        %   dw: increment of the jump displacement vector in the local
        %       system
        %
        % Output:
        %   t : traction vector in the local system [ts, tn]
        %
        function t = stressVct(this, dw, pt)

            Te = this.constitutiveMtrx(dw, pt);
            t  = Te*(pt.strainOld + dw);
            
        end

        %------------------------------------------------------------------
        % Compute the constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : elastic constitutive matrix
        %
        function Te = constitutiveMtrx(this, dw, pt)

            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

        end

        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : elastic constitutive matrix
        %
        function Te = elasticConstitutiveMtrx(this,w)

            % Material parameters
            w  = this.parameters(1);    % Aperture
            mu = this.parameters(2);    % Fluid dynamic viscosity
            
            % Permeability 
            Te = w*w*w/12/mu;

        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [t,Te] = evalConstitutiveModel(this,dw,pt)

            % Constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Stress vector
            t  = Te*(pt.strainOld + dw);

        end

    end
end