%% Code for Johnston and Painter to save relevant data at specific time points.
%% Called by Homing_Script_WorldMap.m

xPosition(tSaveCount,iRepeat) = mean(position(:,1));                        % Mean position (x) of the population.
yPosition(tSaveCount,iRepeat) = mean(position(:,2));                        % Mean position (y) of the population.

distanceToGoal(tSaveCount,iRepeat) = sqrt((xPosition(tSaveCount,iRepeat)-goalLocationX)^2+(yPosition(tSaveCount,iRepeat)-goalLocationY)^2);         % Mean distance of the population to the goal.
clusterMeasure(tSaveCount,iRepeat) = sum(pairDistanceVec)/max((nIndividuals*(nIndividuals-1)),1);                                                   % Measure of population clustering.                                                                      
nNeighbours = zeros(nIndividuals,1);                                        % Current number of observed neighbours.
diffDirection = zeros(nIndividuals,1);                                      % Difference in direction between heading and target.

% Loop over individuals
for i = 1:nIndividuals
    
    neighbours = find(pairDistances(i,:)>0&pairDistances(i,:)<sensingRangeField(lonPosition(i),latPosition(i)));
    nNeighbours(i) = numel(neighbours);                                     % Number of observed neighbours.
    diffDirection(i) = abs(heading(i)+pi - mod(navigationField(position(i,1),position(i,2))+pi,2*pi));
    % Direction between heading and target.
    if diffDirection(i) > pi
        diffDirection(i) = pi - mod(diffDirection(i),pi);
    end
end

meanNeighbours(tSaveCount,iRepeat) = mean(nNeighbours);                     % Average number of observed neighbours.
meanDifferenceDirection(tSaveCount,iRepeat) = mean(diffDirection);          % Average difference between heading and target.
nIndividualsRemaining(tSaveCount,iRepeat) = nIndividuals;                   % Number of individuals yet to arrive at target.
savePositionLat(tSaveCount,individualList) = latPosition';                  % Store latitude.         
savePositionLon(tSaveCount,individualList) = lonPosition';                  % Store longtitude.
savePositionLat(tSaveCount,removalStore) = savePositionLat(max(1,tSaveCount-1),removalStore);   % Remove latitude of arrived individuals.
savePositionLon(tSaveCount,removalStore) = savePositionLon(max(1,tSaveCount-1),removalStore);   % Remove longtitude of arrived individuals.     
tSaveCount = tSaveCount + 1;                                                % Increase counter of number of saved ata points.
directionHist = directionHist + histcounts(mod(heading,2*pi),linspace(0,2*pi,nHistDirection))'; % Generate histogram of headings.

% Check if 90% of individals have arrived at target and if so, store time.
if nIndividuals/nIndividualsStart <= 0.1 && majorityCheck == 0              
    majorityGone(iRepeat) = t;
    majorityCheck = 1;
end