

% Define the input parameter values and saves them into a structure array.
% The user can change the values.

% Unit ball diameter for neighborhood smoothing (number of voxels)
% Can have multiple values but only one value usually works well
BallDiam = 3;

% Neighborhood filtering threshold values are selected automatically based
% on derivative of pass(threshold) function.
% 1st and 2nd values are the maximum and minimum derivative values, 3rd
% value is the number of values defined for filtering threshold based on
% the derivative values
FilterPara = [-1 -0.15 4]; 

% Minimum relative data value for acceptable neighbor in the vessel expansion
MinNei0 = [0.6 0.5 0.4 0.3];

% Minimum data value for a starting point of vessel expansion.
% 1st and 2nd values give the upper and lower proportions, 3rd number of 
% values, if more than 2 then interpolated between the upper and lower
% values.
StartPara = [90 30 4]; 

% Maximum distance between vessel tips to be joined by adding new points
MaxDist = [5 10 15 20]; % voxels

% Maximum angle between directions of vessel tips to be joined by adding new points
MaxAngle = [20 30 40 50]; % degrees

clear input
input.BallDiam = BallDiam;
input.FilterPara = FilterPara;
input.MinNei = MinNei0;
input.StartPara = StartPara;
input.MaxDist = MaxDist;
input.MaxAngle = MaxAngle;