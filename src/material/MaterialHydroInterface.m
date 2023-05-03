%% MaterialInterface class
%
% This class defines an abstract traction-displacement constitutive law
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: January 2023
%%%
%% Class definition
classdef MaterialHydroInterface < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        model       = 'elastic';   
        parameters  = [];
        penal       = false; 
        nStVar      = 0;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = MaterialHydroInterface(model, parameters, penal)
            if (nargin > 0)
                this.model      = model;
                this.parameters = parameters;
                this.penal      = penal;
            end
        end
    end

    %% Public methods
    methods

        %------------------------------------------------------------------
        % Get the number of state variables associated with the model
        function nStVar = getNumberStateVar(this)
            nStVar = this.nStVar;
        end

    end
    
    %% Public methods
    methods(Abstract)
        
        %------------------------------------------------------------------
        % Compute the constitutive matrix
        K = permeabilityMtrx(this, dStrain);

    end
end