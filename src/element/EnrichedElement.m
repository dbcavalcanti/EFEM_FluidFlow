%% EnrichedElement class
%
% This class defines an enriched finite element
%
%% Author
% Danilo Cavalcanti
%
%% History
% @version 1.00
%
% Initial version: December 2022
%%%
%
%% Class definition
classdef EnrichedElement < RegularElement
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        idEnr              = [];            % Vector identifying the intersections
        fracture           = [];            % Object fracture
        glpenr             = [];            % Vector of the enhancement degrees of freedom
        nglpenr            = 0;             % Number of enhancement dof
        subDivInt          = false;         % Flag to apply a subdivision of the domain to define the quadrature points
        jumpOrder          = 1;             % Order of the interpolation of the jump displacement field
        staticCondensation = false          % Flag to apply a static condensation of the additional dofs
    end
    %% Constructor method
    methods
        function this = EnrichedElement(type, node, elem, anm, t, matModel, mat, intOrder, glp, fracture, glpenr, subDivInt, jumpOrder, staticCondensation)
            this = this@RegularElement(type, node, elem, anm, t, matModel, mat, intOrder, glp);
            if nargin > 6
                this.fracture         = fracture;
                this.glpenr           = glpenr;
                if staticCondensation == true
                    this.gle = glp;
                else
                    this.gle = [glp glpenr];
                end
                this.nglp               = length(this.glp);
                this.nglpenr            = length(this.glpenr);
                this.ngle               = length(this.gle);
                this.subDivInt          = subDivInt;
                this.jumpOrder          = jumpOrder;
                this.staticCondensation = staticCondensation;
                this.createResults();
            end
        end
    end
    %% Public methods
    methods

        %------------------------------------------------------------------
        % This function assembles the element fluid-flow matrix and 
        % internal discharge vector
        % 
        % Input:
        %   dPe: vector with increment of the nodal pressure vector
        %        associated with the element. The first nglp components of
        %        this vector are associated with the regular dofs and the
        %        final components are associated to the enrichment dofs.
        %
        % Output:
        %   He : element fluid-flow matrix
        %   qe : element internal discharge vector
        %
        function [He,qe] = elementHeQint(this, dPe)

            % Initialize the element's stiffness matrix
            Hpp = zeros(this.nglp,      this.nglp);
            Hpa = zeros(this.nglp,      this.nglpenr/2);
            Hap = zeros(this.nglpenr/2, this.nglp);
            Haa = zeros(this.nglpenr/2, this.nglpenr/2);

            % Initialize the internal force sub-vectors
            qp = zeros(this.nglp, 1);
            qa = zeros(this.nglpenr, 1);

            % Compute the increment of the enrichment dofs
            dPeEnr = this.computeIncrEnrichmentDofs(dPe);

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
            
                % Compute the B matrix at the int. point and the detJ
                [B, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the enriched B-matrix
                Benr = this.enhancedKinematicCompatibilityMtrx(B,this.intPoint(i).X);
        
                % Compute the increment of the strain vector
                dStrain = B*dPe(1:this.nglp);% + Benr*dPeEnr;
        
                % Compute the stress vector and the constitutive matrix
                [stress,K] = this.intPoint(i).constitutiveModel(dStrain);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness sub-matrices
                Hpp = Hpp + B'   * K * B    * c;
                Hpa = Hpa + B'   * K * Benr * c;
                Hap = Hap + Benr'* K * B    * c;
                Haa = Haa + Benr'* K * Benr * c;

                % Numerical integration of the internal force sub-vectors
                qp = qp + B'    * stress * c;
                qa = qa + Benr' * stress * c;

            end

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [Hfl,Hft,Hftt,Hftb] = this.fracture.elementFluidFlowMtrcs(dPeEnr,this.enrVar);

            % Assemble the element fluid flow matrix
            [He,qe] = this.assembleElemHeQe(Hpp,Hpa,Hap,Haa,Hfl,Hft,Hftt,Hftb,qp,qa);
            
        end

        % -----------------------------------------------------------------
        % Function to assemble the element stiffness matrix and internal
        % force vector.
        % The assembly is based on the flag to apply a static condensation
        % or not.
        function [He,qe] = assembleElemHeQe(this,Hpp,Hpa,Hap,Haa,Hfl,Hft,Hftt,Hftb,fp,fa)

            if this.staticCondensation == true

                He = kaa - kaw*(kww\kwa);
                qe = fa  - kaw*(kww\fw);

            else

                % Auxialiary matrices
                Hfp = - (Hftt + Hftb)*Tp;
                Hfa = -(Hftt*Tat + Hftb*Tab)*T;
                Hff = Hfl + Hft;
                Hpf = zeros();
                Haf = zeros();

                He = [Hpp, Hpa, Hpf;
                      Hap, Haa, Haf;
                      Hfp, Hfa, Hff];

                qe = [fp;fa;ff];

            end
        end

        % -----------------------------------------------------------------
        % Function to update the state variables
        function updateStateVar(this)

            for i = 1:this.nIntPoints
                this.intPoint(i).updateStateVar();
                this.intPoint(i).updateStrainVct();
            end

            this.fracture.updateStateVar();

        end

        % -----------------------------------------------------------------
        % Function to compute the increment of the enrichment dofs.
        % When a static condensation of the additional dofs is applied,
        % this increments must be computed.
        function dPeEnr = computeIncrEnrichmentDofs(this,dUe)

            if this.staticCondensation == true

                dPeEnr = this.staticCondensation_ComputeIncrEnrDofs(dUe);

            else

                dPeEnr = dUe((1+this.nglp):end);

            end

        end

        %------------------------------------------------------------------
        % Compute the increment of the enrichment dof based on the static
        % condensation process.
        function dPeEnr = staticCondensation_ComputeIncrEnrDofs(this,DUe)

            % Newton-Raphson parameters
            maxiter = 10;
            conv    = false;
            tol     = 1.0e-9;
            iter    = 0;

            % Initial increment of the enrichment dofs
            dPeEnr = zeros(this.nglpenr,1);

            % Initialize the iterative process
            while (conv == false) && (iter < maxiter)

                % Compute the stiffness matrix and internal force vector
                [kww,fw] = this.discontinuityKFint(DUe,dPeEnr);

                % Check convergence
                if norm(fw) < tol, conv = true; end

                % Compute the total increment of the enrichment dofs
                dPeEnr = dPeEnr - kww\fw;

                % Update the iteration counter
                iter = iter + 1;

            end
        end

        %------------------------------------------------------------------
        % This function assembles the discontinuity stiffness matrix and 
        % internal force vector
        % 
        % Input:
        %   due: vector with the total increment of the nodal displacement
        %        vector associated with the element.
        %   dPeEnr: vector with the total increment of the nodal displacement 
        %        vector associated with the element.
        %
        % Output:
        %   kww : tangent discontinuity stiffness matrix
        %   fw  : discontinuity internal force vector
        %
        function [kww,fw] = discontinuityKFint(this, due, dPeEnr)

            % Initialize the discontinuity stiffness matrix
            kww = zeros(this.nglpenr, this.nglpenr);

            % Initialize the internal force vector
            fw = zeros(this.nglpenr, 1);

            % Numerical integration of the stiffness matrix components
            for i = 1:this.nIntPoints
            
                % Compute the B matrix at the int. point and the detJ
                [B, detJ] = this.shape.BMatrix(this.node,this.intPoint(i).X);

                % Compute the matrices Gr and Gv
                [Gr, Gv] = this.enhancedStrainCompatibilityMtrcs(B,this.intPoint(i).X);

                % Compute the increment of the strain vector
                dStrain = B*due + Gr*dPeEnr;
        
                % Compute the stress vector and the constitutive matrix
                [stress,D] = this.intPoint(i).constitutiveModel(dStrain);
        
                % Numerical integration coefficient
                c = this.intPoint(i).w * detJ * this.t;
        
                % Numerical integration of the stiffness sub-matrices
                kww = kww + Gv'* D * Gr * c;

                % Numerical integration of the internal force sub-vectors
                fw = fw + Gv' * stress * c;

            end

            % Get the discontinuity stiffness matrix and internal force
            % vector
            [kd, fd] = this.fracture.elementKeFint(dPeEnr,this.enrVar);

            % Add the fracture stiffness contribution
            kww = kww + kd;
            fw  = fw  + fd;
            
        end

        % -----------------------------------------------------------------
        % Compute the matrix that discretizes the real/virtual bounded 
        % enhanced strain field based on the enhancement part of the 
        % displacement field. This matrix is used in the KOS and the KSON 
        % formulation.
        function Benr = enhancedKinematicCompatibilityMtrx(this, B, Xn)

            % Integration point in the cartesian coordinates
            X = this.shape.coordNaturalToCartesian(this.node, Xn);

            % Evaluate the Heaviside function at the point X
            h = this.heavisideFnc(X);

            % Compute the Heaviside matrix
            Hd = this.heavisideMtrx();

            % Compute the mapping matrix
            M = this.elementMappingMtrx();

            % Identity matrix
            I = eye(size(Hd,1));

            % Compute the matrix Gr
            Benr = B*(h*I - Hd)*M;

        end

        % -----------------------------------------------------------------
        % Mapping matrix associated to a element. This matrix is
        % constructed by stacking by rows the mapping matrices evaluated at
        % the element's nodes.
        function Me = elementMappingMtrx(this)
            
            % Initialize the element's mapping matrix
            Me = zeros(this.nnd_el*this.ndof_nd,this.nglpenr/2);

            % Construct the matrix by stacking by-rows the mapping matrix
            % evaluated at each node
            for i = 1:this.nnd_el

                % Compute the mapping matrix at node X
                X  = this.node(i,:);
                Mi = this.fracture.jumpTransmissionMtrx(X);

                % Assemble the element mapping matrix
                Me(i,:) = Mi;
            end

        end

        % -----------------------------------------------------------------
        % Heaviside function evaluated at the point (or set of points) X,
        % wrt to the reference point of the fracture (Xref).
        function h = heavisideFnc(this,X)

            % Reference point at the fracture that crosses the element
            Xref = this.fracture.Xref;

            % Fracture normal vector
            n = this.fracture.n;

            % Relative distance of the nodes of the element wrt to the reference point
            DX = X - repmat(Xref,size(X,1),1);
            
            % Heaviside function evaluated in the nodes of the element
            h = max(sign(DX*n'),0.0);

        end

        % -----------------------------------------------------------------
        % Diagonal matrix with the Heaviside function evaluated at the
        % nodes associeted to each regular dof of the element, wrt to the
        % reference point of the fracture (Xref).
        function Hde = heavisideMtrx(this)
 
            % Heaviside function evaluated in the nodes of the element
            h = this.heavisideFnc(this.node);

            % Create the vector
            hv = repmat(h,[1,this.ndof_nd])';
            
            % Matrix with the Heaviside function evaluated in the node of the dofs
            Hde = diag(reshape(hv,numel(hv),1));

        end

        %------------------------------------------------------------------
        % Function to compute the displacement field inside a given element
        function u = displacementField(this,X)
        %
        % Input:
        %   X   : position vector in the global cartesian coordinate system
        %
        % Output:
        %   u   : displacement vector evaluated in "X"
        
            % Natural coordinate system
            Xn = this.shape.coordCartesianToNatural(this.node,X);
            
            % Vector with the shape functions
            Nm = this.shape.shapeFncMtrx(Xn);

            % Heaviside function evaluated in X
            h = this.heavisideFnc(X);

            % Compute the Heaviside matrix
            Hd = this.heavisideMtrx();

            % Compute the mapping matrix
            M = this.elementMappingMtrx();

            % Get enrichment dofs
            w = this.getEnrichmentDofs();

            % Displacement field
            u = Nm*this.ue(1:this.nglp) + Nm*(h*eye(size(Hd)) - Hd)*M*w;
        
        end

        %------------------------------------------------------------------
        % Function to get the enrichment dofs in a element.
        % Considers only one discontinuity inside an element.
        function we = getEnrichmentDofs(this)

            if this.staticCondensation == true
                we = zeros(this.fracture.ndof,1);
                n  = this.fracture.ndof/this.fracture.ndof_nd;
                for i = 1:n 

                    % Get the displacement jumps at the integration points.
                    % Since we are using a Newton-Cotes quadrature rule,
                    % the integration points are coincident with the nodes
                    % of the discontinuity. 
                    w = this.fracture.intPoint(i).strain;

                    % Transform the enrichment dofs to the global 
                    % coordinate system. [shear normal] => [x y] 
                    R = this.fracture.rotationPointMtrx();
                    w = R'*w;
  
                    % Store the jump of node i
                    we((2*i-1):(2*i),:) = w;
                end

                % Transform w to alpha
                if strcmp(this.enrVar,'alpha')
                    S = this.fracture.transformAlphaToW();
                    we = S*we;
                end
            else
                we = this.ue(this.nglp+1:end);
            end

        end

        %------------------------------------------------------------------
        % Function to get the nodes that will be used to plot the results.
        function [resNodes,fractNodes] = getResultNodes(this)

            % Define the edges that the fracture crosses. The following
            % methodology works for one fracture only, then to take into
            % account multiple fractures, a loop needs to be performed.
            
            % Order the nodes in counterclockwise order
            Nodes = [this.node; this.fracture.node];
            order = this.sortCounterClockWise(Nodes);

            % Get the position of the fracture nodes in the order
            nf1 = find(order == size(this.node,1)+1);
            nf2 = find(order == size(this.node,1)+2);

            % The nodes before and after each node of the fracture define
            % the edge of the element.
            % Edge associated to the first node:
            if nf1 == 1
                edge1 = [order(end) order(nf1+1)];
            elseif nf1 == length(order)
                edge1 = [order(nf1-1) order(1)];
            else
                edge1 = [order(nf1-1) order(nf1+1)];
            end
            % Edge associated to the second node:
            if nf2 == 1
                edge2 = [order(end) order(nf2+1)];
            elseif nf2 == length(order)
                edge2 = [order(nf2-1) order(1)];
            else
                edge2 = [order(nf2-1) order(nf2+1)];
             end

            % Length of each edge:
            len1 = sqrt((Nodes(edge1(2),1)-Nodes(edge1(1),1))^2+(Nodes(edge1(2),2)-Nodes(edge1(1),2))^2);
            len2 = sqrt((Nodes(edge2(2),1)-Nodes(edge2(1),1))^2+(Nodes(edge2(2),2)-Nodes(edge2(1),2))^2);

            % Compute the vectors with the origin at the nodes of the
            % fracture oriented along the edges.
            vf1 = Nodes(edge1(1),:) - this.fracture.node(1,:);
            vf2 = Nodes(edge2(1),:) - this.fracture.node(2,:);
            vf1 = vf1/norm(vf1);
            vf2 = vf2/norm(vf2);

            % Compute the increment vector:
            df1 = (len1/1000)*dot(this.fracture.n,vf1)*vf1;
            df2 = (len2/1000)*dot(this.fracture.n,vf2)*vf2;

            % Find a node along the Omega^plus region along each edge:
            n1_plus = this.fracture.node(1,:) + df1;
            n2_plus = this.fracture.node(2,:) + df2;

            % Find a node along the Omega^minus region along each edge:
            n1_minus = this.fracture.node(1,:) - df1;
            n2_minus = this.fracture.node(2,:) - df2; 
           
            % Get the order of the nodes for the element connectivity be
            % done in a counterclockwise way
            fractNodes = [n1_plus;n1_minus; n2_plus;n2_minus];
            resNodes = [this.node; fractNodes];
            orderResNodes = this.sortCounterClockWise(resNodes);
            resNodes = resNodes(orderResNodes,:);

            % Nodes for the fracture result object
            orderFract = this.sortCounterClockWise(fractNodes);
            fractNodes = fractNodes(orderFract,:);

        end

        %------------------------------------------------------------------
        % Function to create an array of result objects based on the
        % children elements.
        function createResults(this)
            [resNodes,fractNodes] = this.getResultNodes();
            nNodes      = size(resNodes,1);
            nFractNodes = size(fractNodes,1);
            this.result = Result(resNodes ,1:nNodes ,0.0*ones(nNodes,1) ,'');
            this.fracture.result = Result(fractNodes ,1:nFractNodes ,0.0*ones(nNodes,1) ,'');
        end

    end

    methods(Static)

        %------------------------------------------------------------------
        % This function sorts counterclockwise a set of nodes.
        % It uses as a reference point the centroid defined by the nodes.
        function order = sortCounterClockWise(NODE)
            
            % Centroid coordinate
            cx = mean(NODE(:,1));
            cy = mean(NODE(:,2));
            
            % Compute the angle that the relative vector of the vertices 
            % from the centroid has with the horizontal axis
            a = atan2(NODE(:,2) - cy, NODE(:,1) - cx);
            
            % Sort the angles
            [~, order] = sort(a);
            
        end

    end
end