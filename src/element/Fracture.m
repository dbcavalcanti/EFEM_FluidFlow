%% Fracture class
%
% This class defines a fracture element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
% Initially prepared for the course CIV 2801 - Fundamentos de Computação
% Gráfica, 2022, second term, Department of Civil Engineering, PUC-Rio.
%
%% Class definition
classdef Fracture < handle    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        shape      = [];              % Object of the Shape class
        idelem     = 0;               % Identify the element that the fracture belongs
        node       = [];              % Nodes of the fem mesh
        connect    = [];              % Nodes connectivity
        m          = [];              % Tangent orientation vector
        n          = [];              % Normal orientation vector
        ld         = 0.0;             % Fracture length
        Xref       = [];              % Reference point
        t          = 1.0;             % Thickness
        matModel   = 'interfaceFlow'; % Material model
        mat        = [];              % Vector with material properties
        nnd_el     = 2;               % Number of nodes per element
        ndof_nd    = 2;               % Number of dof per node
        ndof       = 1;               % Number of dofs
        glw        = [];              % Vector of the degrees of freedom
        nIntPoints = 2;               % Number of integration points
        intPoint   = [];              % Vector with integration point objects
        result     = [];              % Result object
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = Fracture(node, elem, t, matModel, mat, glw)
            if (nargin > 0)

                this.node     = node;
                this.connect  = elem;
                this.t        = t;
                this.matModel = matModel;
                this.mat      = mat;
                this.glw      = glw;
                this.ndof     = length(glw);
                this.shape    = Shape_Bar();

                % Initialize the geometry properties
                this.initializeGeometry();

                % Initialize the integration points
                this.initializeIntPoints();
            end
        end
    end

    %% Abstract methods
    methods(Abstract)

        % Compute the shape function matrix
        N = shapeFncMtrx(this,xn)
        
        % Compute the jump transmission matrix M
        M = jumpTransmissionMtrx(this,X);

        % Compute the rotation matrix
        R = rotationMtrx(this);

    end  
    %% Public methods
    methods

        %------------------------------------------------------------------
        % Initialize the fracture element geometry. Computes the length,
        % the orientation vectors and the reference point.
        function initializeGeometry(this)

            % Fracture length
            dx = this.node(2,1) - this.node(1,1);
            dy = this.node(2,2) - this.node(1,2);
            this.ld = sqrt(dx.^2 + dy.^2);
            
            % Fracture orientation
            sn = dy./this.ld;
            cs = dx./this.ld;
            
            % Tangential vector to the discontinuity
            this.m = [ cs   sn];  
            
            % Normal vector to the discontinuity
            % Defined considering n = ez x m, where ez = [0 0 1]
            this.n = [-sn   cs];

            % Reference point
            this.Xref = 0.5*(this.node(1,:) + this.node(2,:));
            
        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function initializeIntPoints(this)

            % Get integration points coordinates and weights
            [X,w,this.nIntPoints] = this.shape.getIntegrationPoints();

            % Initialize the integration points objects
            intPts(this.nIntPoints,1) = IntPoint();
            for i = 1:this.nIntPoints

                % Initialize the constitutive model
                if strcmp(this.matModel,'interfaceFlow')
                    constModel = MaterialHydroInterface_CubicLaw(this.mat);
                end

                % Create the integration points
                intPts(i) = IntPoint(X(:,i),w(i),'Interface', constModel);

            end

            this.intPoint = intPts;

        end

        %------------------------------------------------------------------
        % Initialize the elements integration points
        function updateStateVar(this)
            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStrainVct();
            end
        end

        %------------------------------------------------------------------
        % This function assembles the element stiffness matrix and internal
        % force vector
        % 
        % Input:
        %   dUe: vector with increment of the nodal displacement vector
        %        associated with the element
        %
        % Output:
        %   ke : element stiffness matrix
        %   fe : element internal force vector
        %
        function [L1, L2, L3, Hf] = elementFluidFlowMtrcs(this,dPe,elem)

            % Initialize the element stiffness matrix and internal force
            % vector
            L1pp  = zeros(elem.nglp,elem.nglp);
            L1pa  = zeros(elem.nglp,elem.nglp);
            L1aa  = zeros(elem.nglp,elem.nglp);
            L2ppf = zeros(elem.nglp,elem.nglpenr/2);
            L2apf = zeros(elem.nglp,elem.nglpenr/2);
            L3    = zeros(elem.nglpenr/2,elem.nglpenr/2);
            Hf    = zeros(elem.nglpenr/2,elem.nglpenr/2);

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Cartesian coordinates of the integration point 
                X = this.cartesianCoordinate(this.intPoint(i).X);

                % Natural coordinates associated with the continuum element
                % of this point
                Xn = elem.shape.coordCartesianToNatural(elem.node,X);

                % Shape function matrix of the continuum
                N = elem.shape.shapeFncMtrx(Xn);
                
                % Enhanced shape function matrix
                Ntop = elem.topEnhancedShapeFncMtrx(N);
                Nbot = elem.bottomEnhancedShapeFncMtrx(N);

                % Shape function matrix of the fracture
                Nd = this.shapeFncMtrx(this.intPoint(i).X);

                % Gradient of the shape function matrix
                dNdds = this.gradShapeFncMtrx(this.intPoint(i).X);

                % Evaluate the jump at the integration point in the local
                % coordinate system
                dp = Nd * dPe(1:2);
           
                % Compute the stress vector and the constitutive matrix
                Kf = this.intPoint(i).getPermeabilityMtrx(dp);
        
                % Numerical integration term. The determinant is ld/2.
                c = this.intPoint(i).w * this.ld/2 * this.t;

                % Get transverse flow leak-off parameters
                cT = Kf(2); cB = Kf(3);

                % Numerical integration of the coupling matrix associated
                % with the porous-media dofs
                L1pp = L1pp + (N' * N) * (cT + cB) * c;
                L1pa = L1pa + ((N' * Ntop) * cT + (N' * Nbot) * cB) * c;
                L1aa = L1aa + ((Ntop' * Ntop) * cT + (Nbot' * Nbot) * cB) * c;

                % Numerical integration of the coupling matrix associated
                % with the discontinuity dofs
                L2ppf = L2ppf + (N' * Nd) * (cT + cB)  * c;
                L2apf = L2apf + ((Ntop' * Nd) * cT + (Nbot' * Nd) * cB) * c;
        
                % Numerical integration of the fracture fluid-flow matrices
                Hf = Hf  + (dNdds' * dNdds) * Kf(1) * c;
                L3 = L3  + (Nd' * Nd) * (cT + cB)  * c;

            end

            L1 = [ L1pp , L1pa;
                   L1pa', L1aa];

            L2 = [ L2ppf; L2apf ];
            
        end

        %------------------------------------------------------------------
        % This function computes the tangential local coordinate s for a
        % given parametric coordinate in the natural coordinate sistem.
        function [s,X] = tangentialLocCoordinate(this,xn)

            % Cartesian coordinates of the given point
            X = this.shape.coordNaturalToCartesian(this.node,xn);

            % Relative position vector
            DX = X - this.Xref;

            % Tangential coordinate
            s = this.m*DX';

        end

        %------------------------------------------------------------------
        % Compute the cartesian coordinates of a point given the tangential
        % coordinate s.
        function X = cartesianCoordinate(this,s)

            % Cartesian coordinates of the given point
            X = this.Xref + s*this.m;

        end

        %------------------------------------------------------------------
        % This function compute the stress interpolation vector
        function S = stressIntVctFnc(this, shape, node, intpOrder)

            % Initialize the Gram matrix
            S = zeros(shape.getSizeStressIntVct(), intpOrder+1);

            % Get the centroid of the element
            X0 = this.Xref;
%             X0 = shape.computeCentroid(node);
 
            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
           
                % Tangential coordinate and point in the global coordinate
                % system
                [s,X] = this.tangentialLocCoordinate(this.intPoint(i).X);

                % Relative position coordinate
                Xrel = X - X0;

                % Compute the integrand of the stress interpolation vector
                dS = shape.integrandStressIntVct(s,Xrel,intpOrder);
        
                % Numerical integration term. The determinant is ld/2.
                c = this.intPoint(i).w * this.ld/2;
        
                % Numerical integration of the stiffness matrix and the
                % internal force vector
                S = S + dS * c;

            end
            
        end

    end

end