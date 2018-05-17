function VesselPoints = vessel_segmentation(Data,input)

% ---------------------------------------------------------------------
% VESSEL_SEGMENTATION.M     Determines which voxels of the input 
%                           photoacoustic image are part of vessels
%
% Version 1.0.0
% Updated       17 May 2018
% Copyright (C) 2017-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% The photoacoustic image is segmented so that the output image describes 
% level of certainty/reliability of a voxel being part of a vessel. 
%
% Input: 
% Data      Photoacoustic image, 3d-array of intensity values
%
% Output:
% VesselPoints  3d-array similar to input "Data", intensity describes level
%               of reliability of a voxel being part of a vessel

if nargin == 1
    %% Define the parameter values
    % !!!!!! THESE PARAMETERS CAN BE CHANGED
    % Unit ball diameter for neighborhood smoothing (number of voxels)
    BallDiam = 3;
    % Neighborhood filtering threshold values are selected automatically based 
    % on derivative of pass(threshold) function
    FilterPara = [-1 -0.15 4];
    % Minimum relative data value for acceptable neighbor in the vessel expansion
    MinNei0 = [0.6 0.5 0.4 0.3];
    % Minimum data value for a starting point of vessel expansion
    StartPara = [90 30 4]; % upper and lower proportions, number of values
    % Maximum distance between vessel tips to be joined by adding new points
    MaxDist = [5 10 15 20];
    % Maximum angle between directions of vessel tips to be joined by adding new points
    MaxAngle = [20 30 40 50]; % degrees
    
    clear input
    input.BallDiam = BallDiam;
    input.FilterPara = FilterPara;
    input.MinNei = MinNei0;
    input.StartPara = StartPara;
    input.MaxDist = MaxDist;
    input.MaxAngle = MaxAngle;
else
    BallDiam = input.BallDiam;
    FilterPara = input.FilterPara;
    MinNei0 = input.MinNei;
    StartPara = input.StartPara;
    MaxDist = input.MaxDist;
    MaxAngle = input.MaxAngle;
end

tic
%% Prepare data and output
% Normalize the data
Data = double(Data);
if max(max(max(Data))) > 1
    Data = Data/max(max(max(Data)));
end
n = size(Data); % Dimensions of the image data
% Add zero values to the boundaries of the data:
data = zeros(n(1)+2,n(2)+2,n(3)+2);
data(2:n(1)+1,2:n(2)+1,2:n(3)+1) = Data;
% Define the 3d-array for counting classification in different
% segmentations
VesselPoints = zeros(n(1),n(2),n(3),'single'); 

%% Segmentation iterations
for h = 1:length(BallDiam)
    %% Neighbourhood Smoothing
    % Make a ball-like unit convolution kernel
    F = zeros(BallDiam(h),BallDiam(h),BallDiam(h));
    c = ceil(BallDiam(h)/2);
    r = BallDiam(h)/2;
    for i = 1:BallDiam(h)
        for j = 1:BallDiam(h)
            for k = 1:BallDiam(h)
                if norm([i-c j-c k-c]) <= r
                    F(i,j,k) = 1;
                end
            end
        end
    end
    
    % Convolute the data with the kernel
    DataConv = convn(data,F,'same');
    DataConv = DataConv/max(max(max(DataConv))); % normalize
    
    %% Determine the neighbourhood filtering parameters
    % Compute the intensity thresholds for the filtering based on given
    % uniform distribution of proportions passing the filtering 
    % Use only non-zero intensity data
    D = DataConv(:);
    D = D(D > 0);
    ndata = length(D); % number of non-zero voxels in the data
    
    % Compute the proportions of voxels passing the filtering for different
    % threshold values
    PropPassFilt = zeros(401,1);
    for i = 1:401
        I = DataConv >= i*0.001;
        PropPassFilt(i) = nnz(I)/ndata;
    end
    
    % Estimate derivative first with 0.01 steps of threshold
    Der = zeros(400,1);
    for i = 1:400
        if i-5 > 0 && i+5 < 402 
            Der(i) = (PropPassFilt(i+5)-PropPassFilt(i-5))/0.001;
        elseif i-5 <= 0
            Der(i) = (PropPassFilt(i+5)-PropPassFilt(1))/(0.001*(i+4));
        else
            Der(i) = (PropPassFilt(401)-PropPassFilt(i-5))/(0.001*(402-i+4));
        end
    end
    
    % Select the threshold limits based on the derivative
    lower = 1;
    while Der(lower) < FilterPara(1)
        lower = lower+1;
    end
    upper = lower;
    while upper < 400 && Der(upper) < FilterPara(2)
        upper = upper+1;
    end
    
    % Determine the threshold values corresponding to uniformly distributed
    % proportion values
    Proportions = linspace(PropPassFilt(lower),PropPassFilt(upper),FilterPara(3));
    
    FilterThreshold = zeros(1,FilterPara(3));
    FilterThreshold(1) = lower*0.001;
    FilterThreshold(FilterPara(3)) = upper*0.001;
    if FilterPara(3) > 2
        f = lower;
        for i = 2:FilterPara(3)-1
            while PropPassFilt(f) > Proportions(i)
                f = f+1;
            end
            FilterThreshold(i) = f*0.001;
        end
    end
    input.FilterThreshold = FilterThreshold;
            
    for f = 1:length(FilterThreshold)
        %% Neighborhood filtering
        I = DataConv < FilterThreshold(f);
        Data1 = DataConv;
        Data1(I) = 0;
        a = round(nnz(Data1 > 0)/ndata*1000)/10;
        disp(['  ',num2str(a),' % of the points pass the filtering'])
        
        %% Define MinStart parameter values
        MinNei = MinNei0;
        % Use only non-zero intensity data
        D = Data1(:);
        D = D(D > 0);
        
        % The intensity bin count in percents of the total (0.001-width bins)
        IntensityCount = histcounts(D,0:0.001:1)/length(D)*100;
        CumIntensityCount = cumsum(IntensityCount); % Cumulative count
        
        % Define the uniform distibution of proportions
        Proportions = linspace(StartPara(1),StartPara(2),StartPara(3));
        
        % Determine which intensity values correspond to the proportions
        MinStart = Proportions;
        a = length(CumIntensityCount);
        b = 1;
        while b <= StartPara(3)
            while CumIntensityCount(a) > MinStart(b)
                a = a-1;
            end
            MinStart(b) = a/1000;
            b = b+1;
        end
        input.MinStart = MinStart;
        
        for i = 1:length(MinStart)
            
            for j = 1:length(MinNei)
                
                %% Segment the data based on MinNei and MinStart values
                [Q,Vessels,Tips,TipPoints] = segment_data(Data1,MinNei(j),MinStart(i));
                
                % Update the vessel frequency for voxels
                Q = [Q(:,1)-1 Q(:,2)-1 Q(:,3)-1]; % Back to original "coordinates"
                I = [Q(:,1) Q(:,2)-1 Q(:,3)-1]*[1 n(1) n(1)*n(2)]';
                VesselPoints(I) = VesselPoints(I)+length(MaxDist)*length(MaxAngle);
                
                if ~isempty(Tips)
                    %% Add points between tips of vessels
                    NewPoints = combine_vessels(Data1,Q,Vessels,Tips,...
                        TipPoints,MaxDist(end),MaxAngle(end));
                    
                    % Update the vessel frequence for voxels
                    nc = size(NewPoints,1);
                    for k = 1:nc
                        Q = NewPoints{k,1};
                        D = NewPoints{k,2};
                        Q = [Q(:,1)-1 Q(:,2)-1 Q(:,3)-1]; % Back to original "coordinates"
                        I = [Q(:,1) Q(:,2)-1 Q(:,3)-1]*[1 n(1) n(1)*n(2)]'; % linear indexes
                        J = I > 0;
                        I = I(J);
                        n1 = nnz(D(1) <= MaxDist);
                        n2 = nnz(D(2) <= MaxAngle);
                        VesselPoints(I) = VesselPoints(I)+n1*n2;
                    end
                end
            end
        end
        toc
    end
end
VesselPoints = VesselPoints/max(max(max(VesselPoints)));


%% Visualise the vessel segmentation
% Generate point cloud with the computed reliability as its intensity
n = size(VesselPoints);
np = nnz(VesselPoints > 0);
P = zeros(np,3);
Values = zeros(np,1);
np = 0;
for i = 1:n(1)
    for j = 1:n(2)
        for k = 1:n(3)
            if VesselPoints(i,j,k) > 0
                np = np+1;
                P(np,1:3) = [i j k];
                Values(np) = VesselPoints(i,j,k);
            end
        end
    end
end
P = P(1:np,:);

figure(1)
scatter3(P(:,1),P(:,2),P(:,3),20*ones(np,1),Values,'Filled')
axis equal
colorbar
colormap('default')
caxis([0 1])

end % End of main function
