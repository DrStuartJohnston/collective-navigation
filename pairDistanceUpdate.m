%% Code for Johnston and Painter to update distances between pairs of
%% individuals. Called by Homing_Script.m.

pairDistances = zeros(nIndividuals);
pairDistanceVec = pdist(position);                                  % Calculate distances between all pairs of individuals.
pairDistances(triu(ones(nIndividuals)==1,1)) = pairDistanceVec;     % Set pair distances for i =/= j.
pairDistances(tril(ones(nIndividuals)==1,-1)) = pairDistanceVec;    % Set pair distances for i =/= j.     
