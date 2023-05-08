function [SUPP, LOAD] = boundaryConditions(NODE,CoordSupp,CoordLoad,Lx, Ly, Nx, Ny)

SUPP = zeros(size(NODE,1),1);
LOAD = zeros(size(NODE,1),1);

% --- Generate loads
epsx= Lx/Nx/10; epsy= Ly/Ny/10;
for i=1:size(CoordLoad,1)
    Px = CoordLoad(i,1);
    cx = CoordLoad(i,2);
    cy = CoordLoad(i,3);
    if cx<0         % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE(:,2)-cy)<=epsy);
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE(:,1)-cx)<=epsx);
    else
        k= find(abs(NODE(:,1)-cx)<=epsx & abs(NODE(:,2)-cy)<=epsy);
    end
    for j= 1:length(k)
        LOAD(k(j),:)= Px;
    end
end

% --- Generate supports
for i=1:size(CoordSupp,1)
    rtx = CoordSupp(i,1); %if rtx==0,rtx= NaN; end
    cx  = CoordSupp(i,2);
    cy  = CoordSupp(i,3);
    if cx<0          % If cx<0, generate load along x at nodes with y==cy
        k= find(abs(NODE(:,2)-cy)<=epsy);
    elseif cy<0      % If cy<0, generate load along y at nodes with x==cx
        k= find(abs(NODE(:,1)-cx)<=epsx);
    else
        k= find(abs(NODE(:,1)-cx)<=epsx & abs(NODE(:,2)-cy)<=epsy);
    end
    for j= 1:length(k)
        SUPP(k(j),:)= [rtx];
    end
end