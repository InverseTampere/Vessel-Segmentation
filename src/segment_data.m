function [P,Vessels,Tips,TipPoints] = segment_data(Data,MinNei,MinStart)

% ---------------------------------------------------------------------
% SEGMENT_DATA.M     Segments the data into vessels and defines the tips of
%                       the vessels
%
% Version 1.0.0
% Updated       17 May 2018
% Copyright (C) 2017-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Segments a voxel image into vessels based on Minimum starting point
% intensity and region growing based on minimum neighbor intensity.
%
% Input: 
% Data          Photoacoustic image, 3d-array of normalized intensity values
% MinNei        Parameter defining the minimum relative neighbor intensity
%                   for neighbourhood growing in the clustering process
% MinStart      Parameter defining the minimum acceptable intensity for
%                   starting cluster
%
% Output:
% P             3D point cloud of the non-zero data points
% Vessels       Vector giving the vessel/cluster of each point in the point cloud
% Tips          Matrix whose rows contain information of each determined
%                   vessel tip: Location, direction, vessel/cluster it
%                   belongs
% TipPoints     Cell-array of points of P forming the tips
% ---------------------------------------------------------------------

%% Transform the data into point cloud data
Data = double(Data);
[P,Nei,Data] = transform_to_point_cloud(Data);

%% Determine the initial estimate of the vessels
[P,Vessels] = initial_vessels(Data,Nei,MinNei,MinStart,P);

if ~isempty(Vessels)
    
    % Define input parameters for the two covers
    clear inputs
    inputs.PatchDiam = 2;
    inputs.BallRad = 3;
    
    clear inputs2
    inputs2.PatchDiam = 1.5;
    inputs2.BallRad = 2.5;
    
    %% Segment the point cloud
    % Cover the point cloud
    cover = cover_vessels(P,inputs,Vessels);
    nb = size(cover.ball,1);
    
    % Determine a base for each component
    nv = max(Vessels); % number of vessel-components
    Base = cell(nv,1);
    ind = (1:1:nb)';
    vessels = Vessels(cover.center);
    for i = 1:nv
        Vessel = ind(vessels == i);
        L = cellfun('length',cover.neighbor(Vessel));
        [~,J] = min(L);
        Base{i} = Vessel(J);
    end
    % Segment the point cloud
    Forb = false(nb,1);
    segment = vessel_segments(cover,Base,Forb,2);
        
    %% Segment the point cloud again
    % Select the bases again
    Base = select_bases(segment,cover,nv);
    
    % Cover the point cloud again
    cover = cover_vessels(P,inputs2,Vessels);
    
    % Define the bases in terms of the new cover
    for i = 1:nv
        B = unique(cover.BallOfPoint(Base{i}));
        Base{i} = B;
    end
    
    % Segment the point cloud
    Forb = false(size(cover.ball,1),1);
    segment = vessel_segments(cover,Base,Forb,1);
    
    %% Determine the vessel tips
    [Tips,TipPoints] = vessel_tips(P,cover,segment,nv);
    
else
    Tips = zeros(0,7);
    TipPoints = cell(0,1);
end

end % End of main function


function [P,Nei,Data] = transform_to_point_cloud(Data)

% Transform the grid data into point cloud
n = size(Data);
PointInCell = zeros(n(1),n(2),n(3),'uint32');
np = ceil(0.2*prod(n));
P = zeros(np,3);
Data1 = zeros(np,1);
p = 0;
for i = 2:n(1)-1
    for j = 2:n(2)-1
        for k = 2:n(3)-1
            if Data(i,j,k) > 0.001
                p = p+1;
                P(p,:) = [i j k];
                PointInCell(i,j,k) = p;
                Data1(p) = Data(i,j,k);
            end
        end
    end
end
np = p;
P = P(1:p,:);
Data = Data1(1:p); % transform the data into vector

% Generate the neighbors
Nei = cell(np,1);
for i = 1:np
    M = PointInCell(P(i,1)-1:P(i,1)+1,P(i,2)-1:P(i,2)+1,P(i,3)-1:P(i,3)+1);
    Nei{i} = M(M > 0);
end

end % End of subfunction


function [P,Vessels,I] = initial_vessels(Data,Nei,MinNei,MinStart,P)

% Fix the values for minimum acceptable initial vessel (component) size and
% for the maximum number of initial components (large number of components 
% indicates that the segmentation parameters MinNei, MinStart do nor have 
% proper values).
MinCompSize = 20;
MaxNumOfComps = 1000;

% Search the element with the maximum value 
np = size(Data,1);
Vessel = zeros(np,1);
[MaxDataValue,Point] = max(Data);
n = 1; % number of points for the vessel expansion
nv = 0; % vessel index / number of vessels found
% Define vessels as long as there exist suitable starting points
Points = zeros(np,1);
while MaxDataValue > MinStart
    nv = nv+1;
    % Expands as long as suitable neighbors
    t = 1;
    while n > 0
        t0 = t;
        for i = 1:n
            N = Nei{Point(i)};
            I = Data(N) > MinNei*Data(Point(i));
            J = Vessel(N) == 0;
            N = N(I&J);
            if ~isempty(N)
                Vessel(N) = nv;
                m = length(N);
                Points(t:t+m-1) = N;
                t = t+m;
            end
        end
        if t > 1
            Point = Points(t0:t-1);
        end
        n = t-t0;
    end
    if t > 1
        Data(Points(1:t-1)) = 0; % Nullify the segmented vessel points in the data
    else
        Data(Point) = 0;
    end
    [MaxDataValue,Point] = max(Data);
    n = 1;
end

% Change the neighbors to only vessel points
for i = 1:np
    N = Nei{i};
    I = Vessel(N) > 0;
    Nei{i} = N(I);
end

% Determine the vessels as connected components
VesselComps = connected_components(Nei,0,MinCompSize);

% Determine vessels only if there is not too many components
nc = size(VesselComps,1);
if nc < MaxNumOfComps
    Vessels = zeros(np,1,'uint16');
    for i = 1:nc
        Vessels(VesselComps{i}) = i;
    end
else
    Vessels = zeros(1,0);
end
I = Vessels > 0;
P = P(I,:);
Vessels = Vessels(I);

end % End of subfunction


function Base = select_bases(segment,cover,nv)

% Select the base for each vessel-cluster based on the vessel-segmentation 
% of the cluster

% The vessel tips are the tips of the segments
% Determine the vessel-component of each segment
ns = max(size(segment.segments));
VesOfSeg = zeros(ns,1);
for i = 1:nv
    C = segment.ChildSegment{i};
    VesOfSeg(i) = i;
    while ~isempty(C)
        VesOfSeg(C) = i;
        C = vertcat(segment.ChildSegment{C});
    end
end

D = zeros(2*ns,4);
T = cell(2*ns,1);
t = 0;
for i = 1:ns
    Seg = segment.segments{i};
    N = size(Seg,1);
    C = segment.ChildSegment{i};
    if ~isempty(C) && segment.ParentSegment(i,1) == 0
        % First-order segment with children
        L = segment.ParentSegment(C,2);
        if all(L < N)
            % Potential base at the top of the segment as there are no
            % branches there
            t = t+1;
            B = vertcat(cover.ball(vertcat(Seg{max(1,N-2):N})));
            D(t,1) = length(B);
            D(t,2) = 2;
            D(t,3) = i;     D(t,4) = VesOfSeg(i);
            T{t} = Seg{N};
        end
        if all(L > 1)
            % Potential base at the bottom if there are no branches there
            t = t+1;
            B = vertcat(cover.ball(vertcat(Seg{1:min(3,N)})));
            D(t,1) = length(B);
            D(t,2) = 1;
            D(t,3) = i;     D(t,4) = VesOfSeg(i);
            T{t} = Seg{1};
        end
    elseif ~isempty(C)
        % higher-order segment
        L = segment.ParentSegment(C,2);
        if all(L < N)
            % Potential base at the top of the segment as there are no
            % branches there
            t = t+1;
            B = vertcat(cover.ball(vertcat(Seg{max(1,N-2):N})));
            D(t,1) = length(B);
            D(t,2) = 2;
            D(t,3) = i;     D(t,4) = VesOfSeg(i);
            T{t} = Seg{N};
        end
    elseif segment.ParentSegment(i,1) == 0
        % No children and first-order segment (no parent segment)
        % The top
        t = t+1;
        B = vertcat(cover.ball(vertcat(Seg{max(1,N-2):N})));
        D(t,1) = length(B);
        D(t,2) = 2;
        D(t,3) = i;     D(t,4) = VesOfSeg(i);
        T{t} = Seg{N};
        
        % The bottom
        t = t+1;
        B = vertcat(cover.ball(vertcat(Seg{1:min(3,N)})));
        D(t,1) = length(B);
        D(t,2) = 1;
        D(t,3) = i;     D(t,4) = VesOfSeg(i);
        T{t} = Seg{1};
    else
        % No children and higher-order segment
        % The top
        t = t+1;
        B = vertcat(cover.ball(vertcat(Seg{max(1,N-2):N})));
        D(t,1) = length(B);
        D(t,2) = 2;
        D(t,3) = i;     D(t,4) = VesOfSeg(i);
        T{t} = Seg{N};
    end
end
D = D(1:t,:);

% Define the output "Base"
ind = (1:1:t)';
Base = cell(nv,1);
for i = 1:nv
    Ind = ind(D(:,4) == i);
    [~,K] = max(D(Ind,1));
    Ind = Ind(K);
    S = segment.segments{D(Ind,3)};
    if D(Ind,2) == 2
        B = S{end};
    else
        B = S{1};
    end
    Base{i} = vertcat(cover.ball{B});
end

end


function [Tips,TipPoints] = vessel_tips(P,cover,segment,nv)

% The vessel tips are the tips of the segments

% Determine the vessel-component of each segment
ns = max(size(segment.segments));
VesOfSeg = zeros(ns,1);
for i = 1:nv
    C = segment.ChildSegment{i};
    VesOfSeg(i) = i;
    while ~isempty(C)
        VesOfSeg(C) = i;
        C = vertcat(segment.ChildSegment{C});
    end
end

Tips = zeros(2*ns,7);
TipPoints = cell(2*ns,1);
t = 0;
for i = 1:ns
    Seg = segment.segments{i};
    N = size(Seg,1);
    C = segment.ChildSegment{i};
    
    if ~isempty(C) && segment.ParentSegment(i,1) == 0
        % First-order segment with children
        L = segment.ParentSegment(C,2);
        if all(L < N)
            % define tip at the top if there are no branches
            Tip = vertcat(cover.ball{Seg{N}});
            Top = mean(P(Tip,:),1);
            if N > 2
                Bot = mean(P(vertcat(cover.ball{Seg{N-2}}),:),1);
            elseif N > 1
                Bot = mean(P(vertcat(cover.ball{Seg{N-1}}),:),1);
            else
                Nei = vertcat(cover.neighbor{Seg{N}});
                Nei = setdiff(Nei,Seg{N});
                Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
            end
            V = Top-Bot;
            V = V/norm(V);
            t = t+1;
            Tips(t,:) = [Top V VesOfSeg(i)];
            TipPoints{t} = Tip;
        end
        if all(L > 1)
            % define tip at the bottom if there are no branches
            Tip = vertcat(cover.ball{Seg{1}});
            Top = mean(P(Tip,:),1);
            if N > 2
                Bot = mean(P(vertcat(cover.ball{Seg{3}}),:),1);
            elseif N > 1
                Bot = mean(P(vertcat(cover.ball{Seg{2}}),:),1);
            else
                Nei = vertcat(cover.neighbor{Seg{1}});
                Nei = setdiff(Nei,Seg{1});
                Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
            end
            V = Top-Bot;
            V = V/norm(V);
            t = t+1;
            Tips(t,:) = [Top V VesOfSeg(i)];
            TipPoints{t} = Tip;
        end
    elseif ~isempty(C)
        % higher-order segment
        L = segment.ParentSegment(C,2);
        if all(L < N)
            % define tip at the top if there are no branches
            Tip = vertcat(cover.ball{Seg{N}});
            Top = mean(P(Tip,:),1);
            if N > 2
                Bot = mean(P(vertcat(cover.ball{Seg{N-2}}),:),1);
            elseif N > 1
                Bot = mean(P(vertcat(cover.ball{Seg{N-1}}),:),1);
            else
                Nei = vertcat(cover.neighbor{Seg{N}});
                Nei = setdiff(Nei,Seg{N});
                Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
            end
            V = Top-Bot;
            V = V/norm(V);
            t = t+1;
            Tips(t,:) = [Top V VesOfSeg(i)];
            TipPoints{t} = Tip;
        end
    elseif segment.ParentSegment(i,1) == 0
        % No children and first-order segment (no parent segment)
        % Define tip for the top
        Tip = vertcat(cover.ball{Seg{N}});
        Top = mean(P(Tip,:),1);
        if N > 2
            Bot = mean(P(vertcat(cover.ball{Seg{N-2}}),:),1);
        elseif N > 1
            Bot = mean(P(vertcat(cover.ball{Seg{N-1}}),:),1);
        else
            Nei = vertcat(cover.neighbor{Seg{N}});
            Nei = setdiff(Nei,Seg{N});
            Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
        end
        V = Top-Bot;
        V = V/norm(V);
        t = t+1;
        Tips(t,:) = [Top V VesOfSeg(i)];
        TipPoints{t} = Tip;
        
        % Define tip for the bottom
        Tip = vertcat(cover.ball{Seg{1}});
        Top = mean(P(Tip,:),1);
        if N > 2
            Bot = mean(P(vertcat(cover.ball{Seg{3}}),:),1);
        elseif N > 1
            Bot = mean(P(vertcat(cover.ball{Seg{2}}),:),1);
        else
            Nei = vertcat(cover.neighbor{Seg{1}});
            Nei = setdiff(Nei,Seg{1});
            Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
        end
        V = Top-Bot;
        V = V/norm(V);
        t = t+1;
        Tips(t,:) = [Top V VesOfSeg(i)];
        TipPoints{t} = Tip;
    else
        % No childer and higher-order segment
        % Define tip for the top
        Tip = vertcat(cover.ball{Seg{N}});
        Top = mean(P(Tip,:),1);
        if N > 2
            Bot = mean(P(vertcat(cover.ball{Seg{N-2}}),:),1);
        elseif N > 1
            Bot = mean(P(vertcat(cover.ball{Seg{N-1}}),:),1);
        else
            Nei = vertcat(cover.neighbor{Seg{N}});
            Nei = setdiff(Nei,Seg{N});
            Bot = mean(P(vertcat(cover.ball{Nei}),:),1);
        end
        V = Top-Bot;
        V = V/norm(V);
        t = t+1;
        Tips(t,:) = [Top V VesOfSeg(i)];
        TipPoints{t} = Tip;
    end
end
TipPoints = TipPoints(1:t);
Tips = Tips(1:t,:);
end
