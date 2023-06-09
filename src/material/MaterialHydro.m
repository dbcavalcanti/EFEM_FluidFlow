%% Material class
%
% This class defines an abstract stress-strain constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
%% Class definition
classdef MaterialHydro < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        model       = 'saturated';   
        parameters  = [];
        anm         = 'Hydro';
        nStVar      = 0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialHydro(model, parameters, anm)
            this.model      = model;
            this.parameters = parameters;
            this.anm        = anm;
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Get the number of state variables associated with the model 
        % besides the strains and stresses, like plastic deformations or
        % damage
        function nStVar = getNumberStateVar(this)
            nStVar = this.nStVar;
        end

    end

    %% Abstract methods
    methods(Abstract)

        %------------------------------------------------------------------
        % Compute the permeability matrix
        De = permeabilityMtrx(this, dStrain, pt);
        
    end
end