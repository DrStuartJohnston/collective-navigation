%% Code for Johnston and Painter to generate a sample concentration parameter
%% using the look-up table. Modified version of circ_kappa from CircStat.

alphaLookup = alphaLookup(:);

N = length(alphaLookup);

% Calculate R value from observed headings
if N>1
    R = circ_r(alphaLookup,w);
else
    R = alphaLookup;
end

% Calculate kappa value from R based on original circ_kappa.m
if R < 0.53
    kappa = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
    kappa = -.4 + 1.39*R + 0.43/(1-R);
else
    kappa = 1/(R^3 - 4*R^2 + 3*R);
end

% Sample from look-up table if concentration or observation parameters
% sufficiently small
if kappa < 25 && N < 25
    [~,kappaLookUpIndex] = min(abs(kappa-kappaInput));          % Find location of kappa in grid
    cdfSample = rand;
    tmp = find(cdfSample<kappaCDF(:,N-1,kappaLookUpIndex),1);   % Kapppa sample from CDF
    kappa = kappaInput(tmp)+rand*(kappaInput(2)-kappaInput(1)); % Apply noise based on grid spacing
end
