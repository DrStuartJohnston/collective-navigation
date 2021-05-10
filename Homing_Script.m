%% Code to perform the collective navigation described in "Modelling collective
%% navigation via nonlocal communication" by Johnston and Painter. This is
%% highest level script for the idealised information fields. Note this file
%% requires the look-up table for the concentration parameter, which can found
%% at https://melbourne.figshare.com/articles/dataset/kappaCDFLookupTable_mat/14551614.

clear all
close all

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


backgroundFieldType = 'Fixed';   % Choose type of background field, choice of 'Void', 'Fixed','Random','Void', 'Increasing', 'Decreasing', 'Brownian'.
noiseInfluence = 'Information'; % Choose type of noise influence either 'Information' or 'Range'. All results generated with 'Information' except for F9.
flowField = 0;                  % Flow field (unused).
flowDirection = 0;              % Flow direction (unused).
flowVelocity = 0;               % Flow velocity (unused).

totalStepCountLoop = 0;         % Number of reorientation events.
nHistDirection = 60;            % Number of points in histograms.
directionHist = zeros(nHistDirection-1,1); %Predefine direction histogram.
cbar = [linspace(40,115,20)',linspace(36,213,20)',linspace(108,236,20)']/255; %Define colormap.

nIndividualsStart = 100;        % Number of individuals in the simulation at start.
domainWidth = 400;              % Width of the domain.
domainHeight = 300;             % Height of the domain.    
velocity = 1;                   % Speed of individuals.
runTime = 1;                    % Mean reorientation time.
tEnd = 1000;                    % End of simulation.
alpha = 10/20;                  % Weighting of observations for heading calculation.
beta = 10/20;                   % Weighting of observations for concentration calculation.

sensingRange = 20;              % Perceptual range of individuals.
backgroundStrength = 1;         % Background information level.
repulsionDistance = 0;          % Repulsion mechanism (unused).
alignDistance = sensingRange;   % Alignment distance (always = sensing range).
attractDistance = sensingRange; % Attraction mechanism (unused).

goalDistance = 10;              % Distance from goal to be counted as "arrived".
noiseWavelength = 6;            % Frequency of noise structure in the Brownian noise field only.
    
goalLocation = [0,0];           % Location of target.
holeLocation = [125,175];        % Location of information void.
    
navigationField = @(x,y) atan2(goalLocation(2)-y,goalLocation(1)-x) ;       % Direction of target.

%% Main body of simulation, loops over number of realisations.
for iRepeat = 1:nRepeats
    majorityCheck = 0;                      % Check if 90% of population has arrived at the target.
    iRepeat                                 % Print the realisation.

    nIndividuals = nIndividualsStart;       % Number of indivuduals in the simulation.

    t = 0;                                  % Initialise time counter.

    defineBackgroundFields;                 % Define noise and background fields.
    
    initialPosition = zeros(nIndividuals,2);                            % Initial location of individuals.
    initialPosition(:,2) = -20+40*rand(nIndividuals,1);                 % Initial position (y) of individuals.
    initialPosition(:,1) = domainWidth-120+40*rand(nIndividuals,1);     % Initial position (x) of individuals.
    position = initialPosition;                                         % Position of individuals.
    pairDistances = zeros(nIndividuals);                                
    pairDistanceVec = pdist(position);                                  % Calculate distances between all pairs of individuals.
    pairDistances(triu(ones(nIndividuals)==1,1)) = pairDistanceVec;     % Set pair distances for i =/= j.
    pairDistances(tril(ones(nIndividuals)==1,-1)) = pairDistanceVec;    % Set pair distances for i =/= j.
    
    turningTime = exprnd(runTime,nIndividuals,1);                       % Calculate durations of run events.
    timeToUpdate = turningTime;                                         % Calculate time until reorientation events.
    
    heading = zeros(nIndividuals,1);                                    % Headings of individuals.
    
    % Sample individual headings based on inherent information.
    for i = 1:nIndividuals
        heading(i) = circ_vmrnd(navigationField(position(i,1),position(i,2)), ...
            navigationStrengthField(position(i,1),position(i,2)),1);
    end
    
    meanPosition = zeros(tEnd,2);                                       % Calculate mean position of the population
    tSave = linspace(0,tEnd,nSavePoints);                               % Time points where the data will be saved.   
    tSaveCount = 1;                                                     % Count of time points saved.
    totalStepCount = 0;                                                 % Number of steps taken.
    
    % Main loop of the individual simulations, run until end of simulation
    % or all individuals have arrived at the target.
    while t < tEnd && nIndividuals > 0
        
        totalStepCount = totalStepCount + 1;                            % Keep track of total steps taken.
        totalStepCountLoop = totalStepCountLoop + 1;                    % Keep track of overall total steps (i.e. over repeats)
        
        % If sufficiently many steps taken, add extra preallocated vectors.
        if mod(totalStepCountLoop,1e6) == 0
            reorientation = [reorientation;zeros(1e6,1)];
            navError = [navError;zeros(1e6,1)];
            distToGoal = [distToGoal;zeros(1e6,1)];
            neighboursOut = [neighboursOut;zeros(1e6,1)];
        end

        [nextUpdate,nextAgent] = min(timeToUpdate);                     % Calculate next reorientation event time and agent.
        timeToUpdate = timeToUpdate - nextUpdate;                       % Update time to update for all individuals.
        timeElapsed = nextUpdate;                                       % Calculate time step length.
        t = t+nextUpdate;                                               % Update time.
        % Update the position of all individuals. Flow field is not used.
        position = position + velocity*timeElapsed*[cos(heading),sin(heading)] + flowField*flowVelocity*[cos(flowDirection),sin(flowDirection)];
        pairDistanceUpdate;                                             % Update pair distances for all pairs of individuals.
        pairDistances(1:nIndividuals+1:end) = 1e10;                     % Avoid influence of pairs of identical individuals.
        
        % Find individuals within the perceptual range of the individual
        % undergoing reorientation.
        neighbours = find(pairDistances(nextAgent,:)>0&pairDistances(nextAgent,:)<sensingRangeField(position(nextAgent,1),position(nextAgent,2)));
        nNeighbours = numel(neighbours);                                % Number of individuals within perceptual range.
        [minDistance,closestAgent] = min(pairDistances(nextAgent,:));   % Find closest agent.
        oldHeading = heading;                                           % Retain previous heading.
        
        % Calculate sample heading based on inherent information only.
        potentialHeading = circ_vmrnd(navigationField(position(nextAgent,1),position(nextAgent,2)),...
            navigationStrengthField(position(nextAgent,1),position(nextAgent,2)),1);
        
        % Update heading based on other observed individuals if number of
        % neighbours exceeds zero.
        if nNeighbours > 0 && minDistance < sensingRangeField(position(nextAgent,1),position(nextAgent,2))
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
         
        timeToUpdate(nextAgent) = exprnd(runTime,1);                    % New duration of run.
        pairDistances(1:nIndividuals+1:end) = 0;                        % Set pair distances to zeros for identical individuals.
        
        % Storage of data at specific time points
        if t > tSave(tSaveCount)
            saveData;
        end
        
        % Determine which individuals have arrived at the target and remove
        % from simulation.
        removal = [];
        for i = 1:nIndividuals
            if sqrt((position(i,1)-goalLocation(1))^2+(position(i,2)-goalLocation(2))^2) < goalDistance
                removal = [removal;i];
            end
        end
        
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