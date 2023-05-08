function [NODE_D, FRACT] = fractureMesh(NODE, ELEM, XD, SEGD)
% 
% This function divides a fracture segment according to the continuum
% finite element mesh.
%
% It does not assume that the continuum mesh is structured.
%

% Initialize the output arrays
NODE_D = [];
FRACT  = [];

nFract = size(SEGD,1);
nElem  = size(ELEM,1);

globalCounter = 1;

% Loop through the fracture segments --------------------------------------
for i = 1:nFract

    % Points that the define the fracture segment
    p1 = XD(SEGD(i,1),:);
    p2 = XD(SEGD(i,2),:);

    % 
    localCounter = 1;

    % Loop through the continuum elements ---------------------------------
    for el = 1:nElem 

        % Get the number of edges of the element
        nEdges = size(ELEM,2);

        % Get the coordinates of the element (repeat the first one to close
        % the polygon)
        cX = [NODE(ELEM(el,:),1); NODE(ELEM(el,1),1)];
        cY = [NODE(ELEM(el,:),2); NODE(ELEM(el,1),2)];

        % Loop through the edges of the element ---------------------------
        for j = 1:nEdges

            % Points that defined the element edge
            p3 = [cX(j)  , cY(j)];
            p4 = [cX(j+1), cY(j+1)];

            % Evaluate if the segments p1-p2 and p3-p4 intersect
            [flagInt,pint,t] = intersectionSegment(p1,p2,p3,p4);

            % Update the intersection vector points
            if flagInt == true
                NODE_Di(localCounter,:) = pint;
                tNd(localCounter)      = t;
                localCounter           = localCounter + 1;
            end
        end
    end

    % Get the unique nodes of the intersection of the fracture "i" with the
    % mesh
    [tNd,I] = sort(tNd);
    NODE_Di = NODE_Di(I,:);
    [~,I]   = uniquetol(tNd,1e-9);
    NODE_Di = NODE_Di(I,:);

    nAddedNd = size(NODE_Di,1);

    if nAddedNd > 1
        FRACT_i = zeros(nAddedNd-1,2);
        for k = 1:nAddedNd-1
            FRACT_i(k,:)  = [globalCounter globalCounter+1];
            globalCounter = globalCounter + 1;
        end
    end

    % Assemble the global arrays
    NODE_D = [NODE_D; NODE_Di];
    FRACT  = [FRACT;  FRACT_i];

end

end