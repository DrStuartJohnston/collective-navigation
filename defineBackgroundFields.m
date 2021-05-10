%% Code for Johnston and Painter to define the background fields. Called by
%% Homing_Script.m

% Check if noise influences information or perceptual range.
if strcmpi(noiseInfluence,'Information')
    sensingRangeField = @(x,y) sensingRange;
    % Check each type of field.
    if strcmpi(backgroundFieldType,'Brownian')
        backgroundField = brownianNoise(domainWidth,domainHeight,noiseWavelength);
        backgroundField = backgroundField';
        newbackgroundField = zeros(size(backgroundField));
        newbackgroundField(backgroundField<=0.5) = backgroundField(backgroundField<=0.5).^2;
        newbackgroundField(backgroundField>0.5) = 1 - (1-backgroundField(backgroundField>0.5)).^2;
        backgroundField = newbackgroundField;
    end
    if strcmpi(backgroundFieldType,'Fixed')
        navigationStrengthField = @(x,y) backgroundStrength;
    elseif strcmpi(backgroundFieldType,'Random')
        navigationStrengthField = @(x,y) 2*backgroundStrength*rand;
    elseif strcmpi(backgroundFieldType,'Void')
        navigationStrengthField = @(x,y) backgroundStrength*(1-(sqrt(x^2+y^2)>holeLocation(1))*(sqrt(x^2+y^2)<holeLocation(2)));
    elseif strcmpi(backgroundFieldType,'Increasing')
        navigationStrengthField = @(x,y) backgroundStrength*(0.5 + 1.5*(tanh(-0.02*(sqrt((x-goalLocation(1))^2+(y-goalLocation(2))^2)-50))+1)/2);
    elseif strcmpi(backgroundFieldType,'Decreasing')
        navigationStrengthField = @(x,y) backgroundStrength*(0.5 + 1.5*(tanh(0.02*(sqrt((x-goalLocation(1))^2+(y-goalLocation(2))^2)-50))+1)/2);
    elseif strcmpi(backgroundFieldType,'Brownian')
        navigationStrengthField = @(x,y) 2*backgroundStrength*backgroundField(round(x)+1+domainWidth,round(y)+1+domainHeight);
    elseif strcmpi(backgroundFieldType,'Drilling')
        navigationStrengthField = @(x,y) backgroundStrength*(1-sum(noiseMap(x,y))/max(max(noiseGrid)));
    end
elseif strcmpi(noiseInfluence,'Range')
    navigationStrengthField = @(x,y) backgroundStrength;
    % Check each type of field.
    if strcmpi(backgroundFieldType,'Brownian')
        backgroundField = brownianNoise(domainWidth,domainHeight,noiseWavelength);
        backgroundField = backgroundField';
        newbackgroundField = zeros(size(backgroundField));
        newbackgroundField(backgroundField<=0.5) = backgroundField(backgroundField<=0.5).^2;
        newbackgroundField(backgroundField>0.5) = 1 - (1-backgroundField(backgroundField>0.5)).^2;
        backgroundField = newbackgroundField;
    end
    if strcmpi(backgroundFieldType,'Fixed')
        sensingRangeField = @(x,y) sensingRange;
    elseif strcmpi(backgroundFieldType,'Random')
        sensingRangeField = @(x,y) 2*sensingRange*rand;
    elseif strcmpi(backgroundFieldType,'Void')
        sensingRangeField = @(x,y) sensingRange*(1-(sqrt(x^2+y^2)>holeLocation(1))*(sqrt(x^2+y^2)<holeLocation(2)));
    elseif strcmpi(backgroundFieldType,'Increasing')
        sensingRangeField = @(x,y) sensingRange*(0.5 + 1.5*(tanh(-0.02*(sqrt((x-goalLocation(1))^2+(y-goalLocation(2))^2)-50))+1)/2);
    elseif strcmpi(backgroundFieldType,'Decreasing')
        sensingRangeField = @(x,y) sensingRange*(0.5 + 1.5*(tanh(0.02*(sqrt((x-goalLocation(1))^2+(y-goalLocation(2))^2)-50))+1)/2);
    elseif strcmpi(backgroundFieldType,'Brownian')
        sensingRangeField = @(x,y) 2*sensingRange*backgroundField(round(x)+1+domainWidth,round(y)+1+domainHeight);
    elseif strcmpi(backgroundFieldType,'Drilling')
        sensingRangeField = @(x,y) sensingRange*(1-sum(noiseMap(x,y))/max(max(noiseGrid)));
    end
end