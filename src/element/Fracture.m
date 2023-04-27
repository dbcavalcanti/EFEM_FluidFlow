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
% Initially prepared for the course CIV 2801 - Fundamentos de Computa��o
% Gr�fica, 2022, second term, Department of Civil Engineering, PUC-Rio.
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
        N = interpJumpShapeMtrx(this,xn, enrVar)
        
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
                    constModel = MaterialInterface_CubicLaw(this.mat);
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
        function [Hfl,Hft,Hftt,Hftb] = elementFluidFlowMtrcs(this,dUe,enrVar)

            % Initialize the element stiffness matrix and internal force
            % vector
            ke = zeros(this.ndof,this.ndof);
            fe = zeros(this.ndof,1);

            % Compute the rotation matrix
            R = this.rotationMtrx();

            % Transform the enrichment dofs to the local coordinate system.
            % [x y] => [shear normal]
            dUe = R*dUe;

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints

                % Shape function matrix
                N = this.interpJumpShapeMtrx(this.intPoint(i).X,enrVar);

                % Evaluate the jump at the integration point in the local
                % coordinate system
                dw = N * dUe;
           
                % Compute the stress vector and the constitutive matrix
                [td,T] = this.intPoint(i).constitutiveModel(dw);
        
                % Numerical integration term. The determinant is ld/2.
                c = this.intPoint(i).w * this.ld/2 * this.t;
        
                % Numerical integration of the stiffness matrix and the
                % internal force vector
                Hft  = Hft  + N' * (ct + cb)  * N * c;
                Hftt = Hftt + N' * ct  * N * c;
                Hftb = Hftb + N' * cb  * N * c;

            end

            % Rotate to the global coordinate system
            ke = R' * ke * R;
            fe = R' * fe;
            
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