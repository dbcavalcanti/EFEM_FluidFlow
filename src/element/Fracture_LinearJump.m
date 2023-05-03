%% Fracture_ConstantJump class
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: February, 2023
%
%% Class definition
classdef Fracture_LinearJump < Fracture
    %% Constructor method
    methods
        function this = Fracture_LinearJump(node, elem, t, matModel, mat, glw)
            this = this@Fracture(node, elem, t, matModel, mat, glw);
        end
    end
    %% Public methods
    % Implementation of the abstract methods declared in super-class
    methods

        %------------------------------------------------------------------
        % This function computes the matrix of the shape function
        % to evaluate the pressure field along the discontinuity
        function N = shapeFncMtrx(this,xn)

            % Cartesian coordinates of the given point
            X = this.shape.coordNaturalToCartesian(this.node,xn);

            % Relative position vector
            DX = X - this.Xref;

            % Normalized tangential coordinate
            s = this.m*DX';

            % Length of the discontinuity
            ld = this.ld;

            % Shape function matrix
            N = [0.5-s/ld, 0.5+s/ld];

        end

        %------------------------------------------------------------------
        % This function computes the matrix of the derivatives of the shape
        % functions with respect to the tangential coordinate s
        function dNds = gradShapeFncMtrx(this,~)

            % Gradient of the shape function matrix
            dNds = [-1, 1]/this.ld;

        end

        % -----------------------------------------------------------------
        % Compute the jump transmission matrix M. This matrix relates the
        % enrichment degrees of freedom alpha with the enhanced
        % displacement field.
        function M = jumpTransmissionMtrx(this,X)

            % Fracture geometric properties
            m    = this.m;
            ld   = this.ld;
            Xref = this.Xref;

            % Relative position vector
            DX = X - Xref;

            % Mapping matrix (2D)
            M = [0.5  0.5];
            M = M + [-m*DX'/ld , m*DX'/ld];

        end

        %------------------------------------------------------------------
        % This function computes the element's rotation matrix. Change from
        % the local coordinate system mn to the global system xy
        function R = rotationMtrx(this)

            % Rotation of a point
            r = this.rotationPointMtrx();

            % Rotation matrix of the element (2 points)
            R = blkdiag(r,r);

        end

        %------------------------------------------------------------------
        % This function computes the element's rotation matrix. Change from
        % the local coordinate system mn to the global system xy
        function R = rotationPointMtrx(this)

            % Rotation of a point
            R = [ this.m(1)   this.m(2);
                  this.n(1)   this.n(2) ];

        end

        %------------------------------------------------------------------
        % This function compute the stress interpolation vector
        function S = stressIntVct(this, shape, node)
            S = this.stressIntVctFnc(shape, node, 1);
        end
 
    end

    %% Public methods
    methods

        % -----------------------------------------------------------------
        % Mapping matrix of the translation part of the jump displacement
        % field 
        function M = addTranslationMtrx(~,M)

            % Translation matrix
            Mt = [1.0   0.0   0.0   0.0;
                  0.0   1.0   0.0   0.0]; 

            % Add to the jump matrix
            M = M + Mt;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the relative rotation part of the jump 
        % displacement field. Evaluated at one point X.
        function M = addRelativeRotationMtrx(this,X,M)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;
            
            % Relative position vector
            DX = X - Xref;

            % Relative rotation transmission matrix matrix
            Mrot = [ 0.0  0.0   sn*DX(2)  -cs*DX(2);
                     0.0  0.0  -sn*DX(1)   cs*DX(1)];

            % Add to the jump matrix
            M = M + Mrot;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the stretching part of the jump displacement 
        % field. Evaluated at one point X.
        function M = addStretchingMtrx(this,X,M)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;

            % Tangential relative coordinate
            s = m*(X' - Xref');

            % Mapping matrix
            Ms = [ 0.0  0.0  s*cs*cs  s*cs*sn;
                   0.0  0.0  s*cs*sn  s*sn*sn];

            % Add to the jump matrix
            M = M + Ms;

        end

        % -----------------------------------------------------------------
        % Mapping matrix of the stretching part in the normal direction due
        % to the Poisson effect. Evaluated at one point X.
        function M = addPoissonEffectMtrx(this,X,M,nu)
            
            % Fracture geometric properties
            m    = this.m;
            cs   = m(1);
            sn   = m(2);
            Xref = this.Xref;

            % Relative position vector
            DX = X - Xref;

            % Auxiliary coefficients
            a1 = sn*sn*DX(1) - cs*sn*DX(2);
            a2 = cs*cs*DX(2) - cs*sn*DX(1);

            % Mapping matrix
            Mnu = [ 0.0  0.0  -nu*cs*a1  -nu*sn*a1;
                    0.0  0.0  -nu*cs*a2  -nu*sn*a2];

            % Add to the jump matrix
            M = M + Mnu;

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from alpha
        % to w
        function Se = transformAlphaToW(this)

            % Find the tangential coordinates along the discontinuity of
            % the initial and final nodes
            s1 = this.tangentialLocCoordinate(-1.0);
            s2 = this.tangentialLocCoordinate(1.0);

            % Matrix Se
            Se = [-s2    0.0   s1    0.0;
                   0.0   -s2   0.0    s1;
                   1.0   0.0  -1.0   0.0;
                   0.0   1.0   0.0  -1.0]/(s1 - s2);

        end

        % -----------------------------------------------------------------
        % Matrix to transform the enrichment degrees of freedom from w
        % to alpha
        function S = transformWToAlpha(this)

            % Find the tangential coordinates along the discontinuity of
            % the initial and final nodes
            s1 = this.tangentialLocCoordinate(-1.0);
            s2 = this.tangentialLocCoordinate(1.0);

            % Matrix S
            S = [ 1.0   0.0    s1   0.0;
                  0.0   1.0   0.0    s1;
                  1.0   0.0    s2   0.0;
                  0.0   1.0   0.0    s2]/(s1 - s2);

        end 


    end

end