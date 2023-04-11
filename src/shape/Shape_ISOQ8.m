%% Shape_ISOQ8 class
%
% This is class defines the behavior of a quadratic quadrilateral
% isoparametric element.
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version:March 2023
%
classdef Shape_ISOQ8 < Shape
    %% Constructor method
    methods
        function this = Shape_ISOQ8()
            this = this@Shape('ISOQ8');
        end
    end

    %% Public methods: methods defined in the abstract superclass
    methods

        % -----------------------------------------------------------------
        % Evaluate the shape function at a given point X of a linear 
        % quadrilateral isoparametric element.
         function N = shapeFnc(~,Xn)

            % Natural coordinates of the given point
            xi = Xn(1); eta = Xn(2);

            % Shape functions
            N5 = 0.50*(1.0 - xi*xi)*(1.0 - eta);
            N6 = 0.50*(1.0 + xi)*(1.0 - eta*eta);
            N7 = 0.50*(1.0 - xi*xi)*(1.0 + eta);
            N8 = 0.50*(1.0 - xi)*(1.0 - eta*eta);
            N1 = 0.25*(1 - xi)*(1 - eta) - 0.50*(N8 + N5);
            N2 = 0.25*(1 + xi)*(1 - eta) - 0.50*(N5 + N6);
            N3 = 0.25*(1 + xi)*(1 + eta) - 0.50*(N6 + N7);
            N4 = 0.25*(1 - xi)*(1 + eta) - 0.50*(N7 + N8);
            N  = [ N1  N2  N3  N4  N5  N6  N7  N8 ];

         end

         % -----------------------------------------------------------------
         % Get the shape function matrix
         function Nm = shapeFncMtrx(this,Xn)

             % Vector with the shape functions
             N = this.shapeFnc(Xn);
            
             % Shape function matrix
             Nm = [N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)  0.0   N(5)  0.0   N(6)  0.0   N(7)  0.0   N(8)  0.0;
                   0.0   N(1)  0.0   N(2)  0.0   N(3)  0.0   N(4)  0.0   N(5)  0.0   N(6)  0.0   N(7)  0.0   N(8)];

         end

         % -----------------------------------------------------------------
         % Compute the derivatives of the shape functions wrt to the
         % natural coordinates
         function dNdxn = shapeFncDrv(~,Xn)

            % Natural coordinates of the given point
            xi = Xn(1); eta = Xn(2);

            % Derivatives of the shape functions
            dN1_dxi  = eta/4 + xi/2 - (eta*xi)/2 - eta^2/4;
            dN2_dxi  = xi/2 - eta/4 - (eta*xi)/2 + eta^2/4;
            dN3_dxi  = eta/4 + xi/2 + (eta*xi)/2 + eta^2/4;
            dN4_dxi  = xi/2 - eta/4 + (eta*xi)/2 - eta^2/4;
            dN5_dxi  = eta*xi - xi;
            dN6_dxi  = 0.5 - eta*eta/2;
            dN7_dxi  = - xi - eta*xi;
            dN8_dxi  = eta*eta/2 - 0.5;

            dN1_deta = eta/2 + xi/4 - (eta*xi)/2 - xi^2/4;  
            dN2_deta = eta/2 - xi/4 + (eta*xi)/2 - xi^2/4; 
            dN3_deta = eta/2 + xi/4 + (eta*xi)/2 + xi^2/4;
            dN4_deta = eta/2 - xi/4 - (eta*xi)/2 + xi^2/4;
            dN5_deta  = xi*xi/2 - 0.5;
            dN6_deta  = - eta - eta*xi;
            dN7_deta  = 0.5 - xi*xi/2;
            dN8_deta  = eta*xi - eta;

            dNdxn   = [ dN1_dxi   dN2_dxi   dN3_dxi   dN4_dxi   dN5_dxi   dN6_dxi   dN7_dxi   dN8_dxi ;
                        dN1_deta  dN2_deta  dN3_deta  dN4_deta  dN5_deta  dN6_deta  dN7_deta  dN8_deta ];

         end

         % -----------------------------------------------------------------
         % Compute the jacobian matrix
         function J = JacobianMtrx(this,X,Xn)

            % Compute the shape function derivatives wrt to the natural
            % coordinate system
            dNdxn = this.shapeFncDrv(Xn);
              
            % Jacobian matrix
            J = dNdxn * X;

         end

         % -----------------------------------------------------------------
         % Compute the determinant of the jacobian
         function detJ = detJacobian(this,X,Xn)
              
            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);
            detJ = det(J);

         end

         % -----------------------------------------------------------------
         % Compute the strain-displacement matrix
         function [B,detJ] = BMatrix(this,X,Xn)

            % Jacobian matrix
            J = this.JacobianMtrx(X,Xn);

            % Determinant of the Jacobian matrix
            detJ = det(J);

            % Compute the derivatives of the shape functions wrt to the
            % natural coordinate system
            dNdxn = this.shapeFncDrv(Xn);

            % Compute the derivatives of the shape functions wrt to the
            % global cartesian coordinate system
            dNdx = J\dNdxn;

            B = zeros(3,8*2);
            for i = 1:8
                B(1,2*i-1) = dNdx(1,i); 
                B(2,2*i)   = dNdx(2,i);
                B(3,2*i-1) = dNdx(2,i); B(3,2*i) = dNdx(1,i);
            end

         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of an
         %           ISOQ8 element
         %   Xn   : vector with the xi and eta coordinates of a point in 
         %          the natural coordinate system
         %
         % Output:
         %   X : vector with the x and y coordinates of a point in the 
         %       global coordinate system
         function X = coordNaturalToCartesian(this,NODE,Xn)

            % Extract the nodal coordinates
            x = NODE(:,1);
            y = NODE(:,2);
            
            % Vector with the shape functions
            Nv = this.shapeFnc(Xn);
            
            % Initialize output
            X = [0.0, 0.0];
            
            % Interpolation the position
            X(1) = Nv(1)*x(1) +  Nv(2)*x(2) +  Nv(3)*x(3) +  Nv(4)*x(4) + Nv(5)*x(5) +  Nv(6)*x(6) +  Nv(7)*x(7) +  Nv(8)*x(8);
            X(2) = Nv(1)*y(1) +  Nv(2)*y(2) +  Nv(3)*y(3) +  Nv(4)*y(4) + Nv(5)*y(5) +  Nv(6)*y(6) +  Nv(7)*y(7) +  Nv(8)*y(8);

         end

         % -----------------------------------------------------------------
         % Transform a point from the natural coordinate system to the
         % global cartesian coordinate system
         % Input:
         %   NODE : matrix with the x and y coordinates of the nodes of the 
         %          element
         %   X    : vector with the x and y coordinates of a point in the 
         %          global cartesian coordinate system
         %
         % Output:
         %   xi, eta: coordinates xi and eta of the point X in the natural
         %            coordinate system.
         %
         function Xn = coordCartesianToNatural(this,NODE,Xp)
            
            % Initial guess for the natural coordinate point
            Xn = [0.0; 0.0];

            % Cartesian coordinates of the initial guess
            X = this.coordNaturalToCartesian(NODE,Xn)';

            % Parameters for the Newton-Raphson algorithm
            tol      = 1e-8;
            it       = 0;
            nIterMax = 10; 

            % Newton-Raphson iterations
            while norm(X-Xp') > tol

                % Compute the jacobian matrix
                J = this.JacobianMtrx(NODE,Xn);

                % Update the natural coordinate
                Xn = Xn + J\(Xp' - X);

                % Cartesian coordinates of the current point
                X = this.coordNaturalToCartesian(NODE,Xn)';

                % To avoid a infinite loop
                it = it + 1;
                if it > nIterMax, break, end

            end 
            
         end

        %------------------------------------------------------------------
        % Function to get the matrices with the coordinates and weights of
        % integration points in the natural coordinate system. 
        % The definition of the integration points depends on the choice of
        % the order of the quadrature rule and the option for subdividing
        % the domain or not.
        % Output:
        %
        %   X : matrix with the coordinates of the integration points in
        %       the natural coordinate system. Each column of this matrix
        %       gives the coordinates of a point.
        %   W : vector with the weights associated with each point
        %   n : number of integration points
        %
        function [X,W,n] = getIntegrationPoints(this,intOrder,elem)
            if nargin == 2
                subDivInt = false;
            else
                subDivInt = elem.subDivInt;
            end

            if subDivInt == false
                % Compute the integration points
                [w,x] = this.getlineQuadrature(intOrder);

                % Assemble the matrix with the integration points 
                [Xip,Yip] = meshgrid(x,x);
                X = [Xip(:)';Yip(:)'];

                % Assemble the vector with the weights
                [Wxip,Wyip] = meshgrid(w,w);
                W = Wxip(:)' .* Wyip(:)';

                % Number of integration points
                n = size(X,2);

            else
                % Heaviside function evaluated in the nodes of the element
                h = elem.heavisideFnc(elem.node);

                % Get the nodes located at each region of the discontinuity
                Node_minus = [elem.node((h < 1),:); elem.fracture.node];
                Node_plus  = [elem.node((h > 0),:); elem.fracture.node];
                
                % Add the center of each region to the list of nodes 
                Node_minus = [Node_minus; mean(Node_minus(:,1))  mean(Node_minus(:,2))];
                Node_plus  = [Node_plus;  mean(Node_plus(:,1))   mean(Node_plus(:,2))];

                % Apply Delaunay's triangularization of each sub-region
                tri_minus = delaunay(Node_minus);
                tri_plus  = delaunay(Node_plus);

                % Get the integration points of a triangle
                shapeCST = Shape_CST();
                [x,w,nIntPoints] = shapeCST.getIntegrationPoints(intOrder,[]);

                % Total number of integration points
                n = nIntPoints*(size(tri_plus,1) + size(tri_minus,1));

                % Initialize the matrix of coordinates and the vector of
                % weights
                W = zeros(1,n);
                X = zeros(2,n);

                % Fill the matrix X and the vector W
                count = 1;
                for subelem = 1:size(tri_minus,1)
                    for ip = 1:length(w)
                        % Find the global cartesian coordinates of the integration point
                        [Xi] = shapeCST.coordNaturalToCartesian(...
                            Node_minus(tri_minus(subelem,:),:),x(:,ip));

                        % Weight
                        detJ = 1;%shapeCST.detJacobian(Node_minus(tri_minus(subelem,:),:),x(:,ip));
                        W(count) = w(ip)*detJ;
                
                        % Mapping the integration point to the quadrilateral
                        % element natural system
                        X(:,count) = this.coordCartesianToNatural(elem.node,Xi)';

                        count = count + 1;

                    end
                end
                for subelem = 1:size(tri_plus,1)
                    for ip = 1:length(w)
                        % Find the global cartesian coordinates of the integration point
                        [Xi] = shapeCST.coordNaturalToCartesian(...
                            Node_plus(tri_plus(subelem,:),:),x(:,ip));
                        
                        % Weight
                        detJ = 1;%shapeCST.detJacobian( Node_plus(tri_plus(subelem,:),:),x(:,ip));
                        W(count) = w(ip)*detJ;
                
                        % Mapping the integration point to the quadrilateral
                        % element natural system
                        X(:,count) = this.coordCartesianToNatural(elem.node,Xi)';

                        count = count + 1;

                    end
                end

            end

        end

        % -----------------------------------------------------------------
        % Integrand to compute the Gram Matrix
        function dH = integrandGramMtrx(this, node, X)

            X    = this.coordNaturalToCartesian(node,X);
            X0   = this.coordNaturalToCartesian(node,[0.0;0.0]);
            Xrel = X - X0;

            dH = [  1.0          Xrel(1)           Xrel(2);
                   Xrel(1)    Xrel(1)*Xrel(1)   Xrel(2)*Xrel(1);
                   Xrel(2)    Xrel(1)*Xrel(2)   Xrel(2)*Xrel(2)];
        end

        % -----------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function n = getSizeStressIntVct(~)
            n = 3;
        end

        % -----------------------------------------------------------------
        % Integrand to compute the stress interpolation vector
        function dS = integrandStressIntVct(~,s,Xrel,jumpOrder)

            if jumpOrder == 0
                dS = [  1.0;
                      Xrel(1);
                      Xrel(2)];
            elseif jumpOrder == 1
                dS = [  1.0       s;
                      Xrel(1)  s*Xrel(1);
                      Xrel(2)  s*Xrel(2)];
            end
        end

    end

    methods (Static)
        %------------------------------------------------------------------
        % Compute the Gauss integration points for a given order
        function [w,gp] = getlineQuadrature(order)
            % Number of Gauss points
            intOrder = order;
        
            % Define the coefficients of the three-term recurrence
            % relationship
            a = @(intOrder)(2*intOrder+1)./(intOrder+1);
            b = @(intOrder)0;
            c = @(intOrder)intOrder./(intOrder+1);
        
            % Constructe the symmetric tridiagonal matrix
            A = -b(0:intOrder-1)./a(0:intOrder-1); B = sqrt(c(1:intOrder-1)./(a(0:intOrder-2).*a(1:intOrder-1)));
            J = diag(B,1) + diag(A) + diag(B,-1);  
        
            % Compute the eigenvalues and eigenvectors
            [V,D] = eig(J,'vector');
        
            % Save (sorted) points and weights
            [gp,I] = sort(D);
            w = (2*V(1,I).^2)'; 
        
            % Note: The next three lines insure zero is zero and the points
            % and weights are perfectly symmetric
            gp(abs(gp)<10*eps) = 0; 
            gp(ceil(end/2)+1:end) = -flipud(gp(1:floor(end/2)));
            w(ceil(end/2)+1:end) = flipud(w(1:floor(end/2)));

        end  
    end
end