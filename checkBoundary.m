%% Code for Johnston and Painter to check if an individual crosses the boundary
%% in the map-based simulation. Called by Homing_Script_WorldMap.m

% Check if individual is located with coarse polygon approximation of the
% coastline.
crossedBoundary = inpolygon(position(:,1),position(:,2),coarseLandX,coarseLandY);

% If so, abort the jump.
position(crossedBoundary,:) = oldPosition(crossedBoundary,:);