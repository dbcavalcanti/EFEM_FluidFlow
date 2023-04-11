%% MaterialInterface_IsotropicDamage class
%
% This class defines an traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: February 2023
%%%
%
%% Class definition
classdef MaterialInterface_ElastoplasticModMC < MaterialInterface_Elastic       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_ElastoplasticModMC(parameters, penal)
            this = this@MaterialInterface_Elastic(parameters, penal);

            % The elastoplastic model requires the storage of the shear
            % and normal components of the plastic part of the displacement
            % jump.
            this.nStVar = 2; 

        end
 
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [t,T] = evalConstitutiveModel(this,dw,pt)

            % Current total jump
            w = pt.strainOld + dw;

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Trial elastic test
            ttrial = Te * w;

            % Evaluate yield criteria
            f = yieldFunction(this, pt, ttrial);

            % Check the yield criteria
            if f < 1.0e-9
                t = ttr;
                T = Te;
            else
                % Return mapping algorithm based on the Closest Point
                % projection

                % Compute the Algorithmic Constitutive Matrix
                
            end


        end

        %------------------------------------------------------------------
        % Evaluate the yield function
        % The yield surface is a Modified Mohr-Coulomb law.
        % Softening laws are considered for the tensile strength and the
        % cohesion.
        function f = yieldFunction(this, pt, t)

            % Get material parameters
            ft0 = this.parameters(3);
            c0  = this.parameters(4);
            phi = this.parameters(5);
            Gf  = this.parameters(6);
            Gfs = this.parameters(7);

            % Auxiliary variables
            beta   = c0*Gf / (ft0*Gfs);
            tanphi = tan(phi);

            % Components of the traction vector
            ts = t(1); tn = t(2);

            % Get the shear and normal components of the plastic part of
            % the jump displacement vector
            ups = pt.statevarOld(1);
            upn = pt.statevarOld(2);

            % Internal variables
            kappa  = upn + beta * ups;
            kappaS = kappa / beta;

            % Compute the current tensile strength and cohesion
            ft = ft0 * exp(-ft0*kappa/Gf);
            c  = c0  * exp(-c0*kappaS/Gfs);
            
            % Yield function
            f = ts*ts - tn*tn*(ft*ft + 2.0*c*ft*tanphi - c*c)/(ft*ft) - c*c*(1.0 + tanphi*tanphi) + (tn + c*tanphi)^2;

        end

    end

end