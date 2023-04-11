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
classdef MaterialInterface_IsotropicDamage_Linear < MaterialInterface_IsotropicDamage       
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialInterface_IsotropicDamage_Linear(parameters, penal)
            this = this@MaterialInterface_IsotropicDamage(parameters, penal);

            % This isotropic scalar damage model requires two state
            % variables, the maximum shear jump and the maximum (positive)
            % normal jump
            this.nStVar = 2; 

        end
 
    end
    
    %% Public methods
    methods


        %------------------------------------------------------------------
        % Damage scalar value. 
        % Considers a exponential evolution law
        function [d,DdDweq] = scalarDamage(this, w, pt)

            % Get the material parameters
            kn = this.parameters(2);
            ft = this.parameters(3);
            Gf = this.parameters(4);

            % Damage threshold
            w0 = ft / kn;

            % Maximum admissible equivalent jump
            wmax = w0 + 2.0 * Gf / ft;

            % Compute the equivalent displacement jump
            weq0 = this.evaluateStateVariable([0.0,0.0],pt);
            weq  = this.evaluateStateVariable(w,pt);
            weq  = max(weq,weq0);

            % Before damage initiated
            if (weq - w0) <= 0.0

                % Scalar damage
                d = 0.0;

                % Derivative of the scalar damage wrt the variable weq
                DdDweq = 0.0;

            % Damage initiated, but did not evolved
            elseif (weq - weq0) <= 0.0

                % Scalar damage
                d = (wmax / (wmax - w0)) * (1.0 - w0/weq0); 

                if d > 0.99999
                    d = 1.0;
                end

                % Derivative of the scalar damage wrt the variable weq
                DdDweq = 0.0;

            % Damage initiated and evolved
            else

                % Scalar damage
                d = (wmax / (wmax - w0)) * (1.0 - w0/weq);

                if d > 0.99999
                    d = 1.0;
                end
    
                % Derivative of the scalar damage wrt the variable weq
                DdDweq =  (w0*wmax)/(weq*weq*(wmax - w0)); 

            end

        end

        %------------------------------------------------------------------
        % Update the scalar internal variables and the current maximum
        % values of the normal and shear jumps.
        function weq = evaluateStateVariable(this,w,pt)
            
            % Compute the maximum shear and normal jumps in the history of
            % this integration point
            wsMax = max(pt.statevarOld(1),abs(w(1)));
            wnMax = max(pt.statevarOld(2),max(w(2),0.0));

            % Shear contribution factor
            beta = this.parameters(5);

            % Equivalent displacement jump
            weq = wnMax + beta * wsMax;

            % Save the new state variables
            pt.statevar(1) = wsMax;
            pt.statevar(2) = wnMax;

        end

        %------------------------------------------------------------------
        function Td = damageConstitutiveMatrix(this,w,DdDweq)

            % Get the material parameters
            ks = this.parameters(1); 
            kn = this.parameters(2);

            % Shear contribution factor
            beta = this.parameters(5);

            % Jump shear and normal components
            ws = w(1); wn = w(2);

            % Constitutive matrix
            Td = DdDweq * [ ks*ws*beta*sign(ws)  ks*ws;
                            kn*wn*beta*sign(ws)  kn*wn ];
        end

    end

end