%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with linear-elastic analysis.
%
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Anl_Linear(result)
            this = this@Anl('Linear',result);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function process(this,mdl)

            % Compute the global stiffness matrix
%             U0 = mdl.U;       % To check the computation of the stresses
            [K,~] = mdl.globalKF(mdl.U); 

            % Add contribution of the nodal forces to the external force
            % vector
            F = mdl.addNodalLoad(mdl.F);
            
            % Solve linear-elastic analysis
            [mdl.U, mdl.F] = this.solveSystem(mdl,K,F,mdl.U);

            % To check the computation of the stresses 
%             [~,~] = mdl.globalKF(mdl.U-U0); 


        end
    end
    
end