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
classdef MaterialInterface_IsotropicDamage < MaterialInterface_Elastic       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_IsotropicDamage(parameters, penal)
            this = this@MaterialInterface_Elastic(parameters, penal);

            % This isotropic scalar damage model requires two state
            % variables, the maximum shear jump and the maximum (positive)
            % normal jump
            this.nStVar = 2; 

        end
 
    end
    
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Compute the traction vector
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

            % Secant constitutive matrix
            Tsec = this.secantConstitutiveMtrx(pt.strainOld + dw);

            % Stress vector
            t  = Tsec*(pt.strainOld + dw);
            
        end

        %------------------------------------------------------------------
        % Compute the tangent constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : tangent constitutive matrix
        %
        function T = constitutiveMtrx(this, dw, pt)

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            d = this.scalarDamage(dw, pt);

            % Tangent constitutive matrix
            T = (1.0 - d)*Te;

        end

        %------------------------------------------------------------------
        % Compute the secant constitutive matrix
        %
        % Input:
        %   w : updated jump displacement vector in the local system [ws, wn]
        %
        % Output:
        %   Te : secant constitutive matrix
        %
        function Tsec = secantConstitutiveMtrx(this, dw, pt)

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            d  = this.scalarDamage(dw, pt);

            % Secant constitutive matrix
            Tsec = (1.0 - d)*Te;

        end

        %------------------------------------------------------------------
        % Compute the stress vector and the constitutive
        % matrix
        function [t,T] = evalConstitutiveModel(this,dw,pt)

            % Current total jump
            w = pt.strainOld + dw;

            % Elastic constitutive matrix
            Te = this.elasticConstitutiveMtrx(pt.strainOld + dw);

            % Scalar damage
            [d,DdDweq]  = this.scalarDamage(w, pt);

            % Secant constitutive matrix
            Tsec = (1.0 - d)*Te;

            % Stress vector
            t  = Tsec * w;

            % Part of the tangent stiffness matrix associated with the
            % damage
            Td = this.damageConstitutiveMatrix(w,DdDweq);

            % Tangent constitutive matrix
            T = Tsec - Td;

        end

    end

end