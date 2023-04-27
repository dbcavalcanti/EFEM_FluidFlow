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
classdef Material_Saturated < Material  
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Material_Saturated(parameters, anm)
            this = this@Material('saturated',parameters, anm);
            this.nStVar = 0;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector
        function stress = stressVct(this, dStrain, pt)
            De      = this.constitutiveMtrx();
            stress  = De*(pt.strainOld + dStrain);
        end
        
        %------------------------------------------------------------------
        % Compute the elastic constitutive matrix
        function K = constitutiveMtrx(this, ~, ~, ~, ~)

            % Conductivity
            k = this.parameters(1);

            % Fluid dynamic viscosity
            mu = this.parameters(2);

            % Constitutive matrix
            K = [k 0.0; 0.0 k]/mu;

        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [vm, K] = evalConstitutiveModel(this,dStrain,pt)

            % Conductivity
            k = this.parameters(1);        

            % Constitutive matrix
            K = [k 0.0; 0.0 k];

            % Stress vector
            vm  = K*(pt.strainOld + dStrain);

        end


    end
end