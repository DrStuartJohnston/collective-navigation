%% Code for Johnston and Painter to define the Brownian noise field. Called by
%% defineBackgroundFields.m

function interpNoise = brownianNoise(domainWidth,domainHeight,noiseWavelength)

%% Parameters

m = floor(4*domainHeight/(2^noiseWavelength))+1; %Number of node points in y
n = floor(4*domainWidth/(2^noiseWavelength))+1; %Number of node points in x
smoothing_Iterations = 1; %Number of iterations to smooth noise

min_X = -domainWidth; %Minimum x value
max_X = 3*domainWidth; %Maximum x value
min_Y = -domainHeight; %Minimum y value
max_Y = domainHeight; %Maximum y values

%% Initial Calculations

x_Centre = (max_X-min_X)/2; %Centre of x domain
y_Centre = (max_Y-min_Y)/2; %Centre of y domain
y_Values = linspace(min_Y,max_Y,m); %Domain of y
x_Values = linspace(min_X,max_X,n); %Domain of x

noise = zeros([m,n]);

%% Generate Fractal Noise

w = sqrt(m*n);
i = 0;
while w > 3
    i = i + 1;
    d = interp2(randn(ceil((m-1)/(2^(i-1))+1),ceil((n-1)/(2^(i-1))+1)), i-1, 'spline');
    noise = noise + i * d(1:m, 1:n);
    w = w - ceil(w/2 - 1);
end
noise = (noise - min(min(noise))) ./ (max(max(noise)) - min(min(noise)));

%% Smoothing

z1 = 2:m-1;
z2 = 2:n-1;

for ii = 1:smoothing_Iterations %Uses five point hat filter
    smoothed_Noise = noise;
    smoothed_Noise(z1,z2) =(noise(z1,z2)+noise(z1-1,z2)+noise(z1+1,z2)+noise(z1,z2-1)+noise(z1,z2+1))/5;
    smoothed_Noise(1,z2) = (noise(1,z2)+noise(2,z2)+noise(1,z2-1)+noise(1,z2+1))/4;
    smoothed_Noise(m,z2) = (noise(m,z2)+noise(m-1,z2)+noise(m,z2-1)+noise(m,z2+1))/4;
    smoothed_Noise(z1,1) = (noise(z1,1)+noise(z1-1,1)+noise(z1+1,1)+noise(z1,2))/4;
    smoothed_Noise(z1,n) = (noise(z1,n)+noise(z1-1,n)+noise(z1-1,n)+noise(z1,n-1))/4;
    smoothed_Noise(1,1) = (noise(1,1)+noise(2,1)+noise(1,2))/3;
    smoothed_Noise(m,1) = (noise(m,1)+noise(m-1,1)+noise(m,2))/3;
    smoothed_Noise(1,n) = (noise(1,n)+noise(2,n)+noise(1,n-1))/3;
    smoothed_Noise(m,n) = (noise(m,n)+noise(m-1,n)+noise(m,n-1))/3;
    noise = smoothed_Noise;
end

interpNoise = interp2(repmat(x_Values,numel(y_Values),1),repmat(y_Values',1,numel(x_Values)),noise,min_X:max_X,[min_Y:max_Y]');
