%% Code to perform the collective navigation described in "Modelling collective
%% navigation via nonlocal communication" by Johnston and Painter. This is
%% highest level script for the "real world" information fields. Note this file
%% requires the look-up table for the concentration parameter, which can found
%% at https://melbourne.figshare.com/articles/dataset/kappaCDFLookupTable_mat/14551614.

clear all
close all

% Choose appropriate terrain map.
terrainMap = 'EastAtlantic'; % 'NorthSea' or 'EastAtlantic'

load(sprintf('%sBackground.mat',terrainMap)); % Load appropriate background field.
% Choices are "NorthSeaBackground" "EastAtlanticBackground" "EastAtlanticBackground20pcNoise" "EastAtlanticBackground40pcNoise" "EastAtlanticBackground100pcNoise" 

nRepeats = 10;                                                  % Number of realisations of the model.
nSavePoints = 501;                                              % Number of time points to save model output.

load('kappaCDFLookupTable.mat');                                % Load the lookup table for estimating the vM concentration parameter.

finalTime = zeros(nRepeats,1);                                  % Time for all individuals to arrive at the goal.
xPosition = zeros(nSavePoints,nRepeats);                        % Mean position (x) of the population.
yPosition = zeros(nSavePoints,nRepeats);                        % Mean position (y) of the population.
clusterMeasure = zeros(nSavePoints,nRepeats);                   % Measure of clustering of the population.
meanNeighbours = zeros(nSavePoints,nRepeats);                   % Mean number of neighbours within perceptual range.
distanceToGoal = zeros(nSavePoints,nRepeats);                   % Mean distance to goal of the population.
meanDifferenceDirection = zeros(nSavePoints,nRepeats);          % Mean error in heading relative to target.
nIndividualsRemaining = zeros(nSavePoints,nRepeats);            % Number of individuals remaining in the simulation (i.e. yet to arrive at goal).
majorityGone = zeros(nRepeats,1);                               % Time for 90% of the individuals to arrive at the goal.

backgroundFieldType = 'Drilling';   %Type of background field - only drilling here.
noiseInfluence = 'Information'; % Choose type of noise influence either 'Information' or 'Range'. All results generated with 'Information' except for F9.
flowField = 0;                  % Flow field (unused).
flowDirection = 0;              % Flow direction (unused).
flowVelocity = 0;               % Flow velocity (unused).

totalStepCountLoop = 0;         % Number of reorientation events.
nHistDirection = 60;            % Number of points in histograms.
directionHist = zeros(nHistDirection-1,1); %Predefine direction histogram.
cbar = [linspace(40,115,20)',linspace(36,213,20)',linspace(108,236,20)']/255; %Define colormap.
projection = projcrs(54030,'Authority','ESRI');                             % Appropriate lat-lon x-y projection.

% Calculate the relevant latitude and longtitude ranges based on terrain
% map.
relevantLand = (coarseLat>(min(latLimit)-5)&coarseLat<(max(latLimit)+5)&coarseLon>(min(lonLimit)-5)&coarseLon<(max(lonLimit)+5))|isnan(coarseLat);

% Calculate x,y co-ordinates based on projection of latitude and
% longtitude.
[coarseLandX,coarseLandY] = projfwd(projection,coarseLat(relevantLand),coarseLon(relevantLand));

nIndividualsStart = 100;                            % Number of individuals in the simulation at start.
domainWidth = lonLimit(2)-lonLimit(1);              % Width of the domain.
domainHeight = latLimit(2)-latLimit(1);             % Height of the domain.    
velocity = 8000;                                    % Speed of individuals (m/h).
runTime = 0.5;                                      % Mean reorientation time (h).
tEnd = 5000;                                        % End of simulation (h).
alpha = 10/20;                                      % Weighting of observations for heading calculation.
beta = 10/20;                                       % Weighting of observations for concentration calculation.

sensingRange = 114000;                              % Perceptual range of individuals (m).
backgroundStrength = 1;                             % Background information level.
repulsionDistance = 0;                              % Repulsion mechanism (unused).
alignDistance = sensingRange;                       % Alignment distance (always = sensing range).
attractDistance = sensingRange;                     % Attraction mechanism (unused).

goalDistance = 50000;                               % Distance from goal (m) to be counted as "arrived".

if strcmpi(terrainMap,'NorthSea')
    [goalLocationX,goalLocationY] = projfwd(projection,62,-1);  % Calculate x,y position of goal location.
    startLat = 55;                                              % Starting latitude for north sea map.
    startLon = 3;                                               % Starting longtitide for north sea map.
elseif strcmpi(terrainMap,'EastAtlantic')
    [goalLocationX,goalLocationY] = projfwd(projection,65,5); 	% Calculate x,y position of goal location.
    startLat = 50;                                              % Starting latitude for east atlantic map.
    startLon = -19;                                             % Starting longtitide for east atlantic map.
end

navigationField = @(x,y) atan2(goalLocationY-y,goalLocationX-x) ;           % Direction of target.

for iRepeat = 1:nRepeats
    majorityCheck = 0;                      % Check if 90% of population has arrived at the target.
    iRepeat                                 % Print the realisation.

    nIndividuals = nIndividualsStart;       % Number of indivuduals in the simulation.

    t = 0;                                  % Initialise time counter.

    defineBackgroundFields;                 % Define noise and background fields.

    savePositionLon = zeros(nSavePoints,nIndividuals);          % Store longtitude position.
    savePositionLat = zeros(nSavePoints,nIndividuals);          % Store latitude position.
    removalStore = [];                                          % Keep track of individuals to remove.
    
    initialPosition = zeros(nIndividuals,2);
    initialPosition(:,2) = startLat+rand(nIndividuals,1);       % Initial latitude of individuals.
    initialPosition(:,1) = startLon+rand(nIndividuals,1);       % Initial longtitude of individuals.
    individualList = 1:nIndividuals;                            % List of individuals.
    lonPosition = initialPosition(:,1);                         % Longtitude position of individuals.
    latPosition = initialPosition(:,2);                         % Latitude position of individuals.
    pairDistances = zeros(nIndividuals);
    [xPositionStart,yPositionStart] = projfwd(projection,initialPosition(:,2),initialPosition(:,1)); %x,y locations of individuals.
    position = [xPositionStart,yPositionStart];                         % Position of individuals.
    pairDistanceVec = pdist(position);                                  % Calculate distances between all pairs of individuals.
    pairDistances(triu(ones(nIndividuals)==1,1)) = pairDistanceVec;     % Set pair distances for i =/= j.
    pairDistances(tril(ones(nIndividuals)==1,-1)) = pairDistanceVec;    % Set pair distances for i =/= j.
    
    turningTime = exprnd(runTime,nIndividuals,1);                       % Calculate durations of run events.
    timeToUpdate = turningTime;                                         % Calculate time until reorientation events.
    
    heading = zeros(nIndividuals,1);                                    % Headings of individuals.
    
    % Sample individual headings based on inherent information.
    for i = 1:nIndividuals
        heading(i) = circ_vmrnd(navigationField(position(i,1),position(i,2)), ...
            navigationStrengthField(lonPosition(i),latPosition(i)),1); 
    end
    
    meanPosition = zeros(tEnd,2);                                       % Calculate mean position of the population
    tSave = linspace(0,tEnd,nSavePoints);                               % Time points where the data will be saved.   
    tSaveCount = 1;                                                     % Count of time points saved.
    totalStepCount = 0;                                                 % Number of steps taken.
    
    % Main loop of the individual simulations, run until end of simulation
    % or all individuals have arrived at the target.
    while t < tEnd && nIndividuals > 0
        
        totalStepCount = totalStepCount + 1;                            % Keep track of total steps taken.

        [nextUpdate,nextAgent] = min(timeToUpdate);                     % Calculate next reorientation event time and agent.
        timeToUpdate = timeToUpdate - nextUpdate;                       % Update time to update for all individuals.
        timeElapsed = nextUpdate;                                       % Calculate time step length.
        t = t+nextUpdate;                                               % Update time.
        oldPosition = position;                                         % Retain previous positions.
        
        % Update the position of all individuals. Flow field is not used.
        position = position + velocity*timeElapsed*[cos(heading),sin(heading)] + flowField*flowVelocity*[cos(flowDirection),sin(flowDirection)];
        
        % Check that no individuals crossed onto land.
        checkBoundary;
        [latPosition,lonPosition] = projinv(projection,position(:,1),position(:,2)); %Calculate latitude and longtitude from x,y positions
        pairDistanceUpdate;                                             % Update pair distances for all pairs of individuals.
        pairDistances(1:nIndividuals+1:end) = 1e10;                     % Avoid influence of pairs of identical individuals.
        
        % Find individuals within the perceptual range of the individual
        % undergoing reorientation.
        neighbours = find(pairDistances(nextAgent,:)>0&pairDistances(nextAgent,:)<sensingRangeField(lonPosition(nextAgent),latPosition(nextAgent)));
        nNeighbours = numel(neighbours);                                % Number of individuals within perceptual range.
        [minDistance,closestAgent] = min(pairDistances(nextAgent,:));   % Find closest agent.
        
        % Calculate sample heading based on inherent information only.
        potentialHeading = circ_vmrnd(navigationField(position(nextAgent,1),position(nextAgent,2)),...
            navigationStrengthField(lonPosition(nextAgent),latPosition(nextAgent)),1);
        
        % Update heading based on other observed individuals if number of
        % neighbours exceeds zero.
        if nNeighbours > 0 && minDistance < sensingRangeField(lonPosition(nextAgent),latPosition(nextAgent))
            % Repulsion mechanism unused.
            if minDistance < repulsionDistance
                heading(nextAgent) = atan2(-position(closestAgent,2)+position(nextAgent,2), ...
                    -position(closestAgent,1)+position(nextAgent,1));
            % Alignment mechanism.    
            elseif minDistance < alignDistance
                bestGuessHeading = circ_mean([circ_mean(heading(neighbours));potentialHeading],[1-alpha;alpha]);    % MLE of heading.
                alphaLookup = [heading(neighbours);potentialHeading];                                               % Set of observed headings.
                w = [(1-beta)*ones(size(neighbours'));beta*nNeighbours];                                            % Weighting of observed headings.
                circ_kappa_script;                                                                                  % Calculate estimate of concentration parameter.
                bestGuessStrength = kappa;                                                                          % Estimate of concentration parameter.
                heading(nextAgent) = circ_vmrnd(bestGuessHeading,bestGuessStrength,1);                              % Set new heading.
            % Attraction mechanism unused.
            elseif minDistance < attractDistance
                heading(nextAgent) = atan2(mean(position(neighbours,2))-position(nextAgent,2),...
                    mean(position(neighbours,1))-position(nextAgent,1));
            end
        else
            heading(nextAgent) = potentialHeading;
        end
        
        timeToUpdate(nextAgent) = exprnd(runTime,1);                        % New duration of run.
        pairDistances(1:nIndividuals+1:end) = 0;                            % Set pair distances to zeros for identical individuals.
        
        % Storage of data at specific time points
        if t > tSave(tSaveCount)
            saveDataWorldMap;
        end
        
        % Determine which individuals have arrived at the target and remove
        % from simulation.
        removal = [];
        for i = 1:nIndividuals
            if sqrt((position(i,1)-goalLocationX)^2+(position(i,2)-goalLocationY)^2) < goalDistance
                removal = [removal;i];
            end
        end
        removalStore = [removalStore;individualList(removal)'];             % Store list of removed individuals.
        individualList(removal) = [];                                       % Remove individuals from list.
        position(removal,:) = [];                                           % Remove individuals from position.
        heading(removal) = [];                                              % Remove individuals from heading.
        timeToUpdate(removal) = [];                                         % Remove individuals from reorientation.
        nIndividuals = nIndividuals - numel(removal);                       % Number of individuals remaining.
    end
    
    finalTime(iRepeat) = t;                                                 % Final time in the simulation.
    
end

xPositionMean = mean(xPosition,2);                                          % Mean of average x position across realisation loop.
clusterMeasure = mean(clusterMeasure,2);                                    % Mean of clustering across realisation loop.
distanceToGoal = mean(distanceToGoal,2);                                    % Mean of average distance to goal across realisation loop.
meanNeighbours = mean(meanNeighbours,2);                                    % Mean of average number of neighbours across realisation loop.
meanDifferenceDirection = mean(meanDifferenceDirection,2);                  % Mean of difference between heading and target across realisation loop.
nIndividualsRemaining = mean(nIndividualsRemaining,2);                      % Mean of number individuals remaining across realisation loop.

clear kappaCDF                                                              % Clear CDF to avoid saving over and over.

% Plot relevant results
plotResults;