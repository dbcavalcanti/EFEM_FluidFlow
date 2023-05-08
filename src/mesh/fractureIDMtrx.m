function IDenr= fractureIDMtrx(NODE,ELEM,NODE_D,FRACT)

% Initialize the matrix
IDenr = zeros(size(ELEM,1),size(FRACT,1)); 

nFract  = size(FRACT,1);
nElem   = size(ELEM,1);
counter = 0;

PINT = [];

% Loop through the fracture segments --------------------------------------
for i = 1:nFract

    % Points that the define the fracture segment
    p1 = NODE_D(FRACT(i,1),:);
    p2 = NODE_D(FRACT(i,2),:);

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
            [flagInt,pint] = intersectionSegment(p1,p2,p3,p4);

            % Update the intersection vector points
            if flagInt == true
                PINT = [PINT;pint];
                counter = counter + 1;
            end

        end
        
        % If the fracture is crossing the element domain, it will have two
        % intersection points
        if size(unique(PINT,'rows')) == 2
            IDenr(el,i) = 1;
        end
        PINT = [];
        counter = 0;
    end
end

end