function cover = cover_vessels(P,inputs,Vessels)

% ---------------------------------------------------------------------
% COVER_VESSELS.M       Creates cover sets and their neighbor-relation for 
%                           a point cloud
%
% Version 1.0.0
% Latest update     17 May 2018
% Copyright (C) 2017-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Partitions the point cloud with small sets such that each set contains 
% only points from a single vessel network. 
%
% Inputs: 
% P         Point cloud
% inputs    Input stucture, the following fields are needed:
%   PatchDiam   Minimum distance between centers of cover sets; i.e. the
%                   minimum diameter of a cover set. If "RelSize" is given
%                   as input, then there is minimum and maximum PatchDiam
% 	BallRad     Radius of the balls used to generate the cover sets, these 
%                   balls are also used to determine the neighbors and the 
%                   cover set characteristics
% Vessels   Vector giving the vessel network the point belongs
%
% Outputs:
% cover     Structure array containing the followin fields:
%   ball        Cover sets, (n_sets x 1)-cell
%   center      Center points of the cover sets, (n_sets x 1)-vector
%   neighbor    Neighboring cover sets of each cover set, (n_sets x 1)-cell
%   BallOfPoint The cover set the point belons, (n_points x 1)-vector

if ~isa(P,'double')
    P = double(P);
end

%% Large balls and centers
np = size(P,1);
Ball = cell(np,1);  % the large balls, used to generate the cover sets and their neighbors
Cen = zeros(np,1,'uint32');  % the center points of the balls/cover sets
NotExa = true(np,1); % the points not yet examined
Dist = 1e8*ones(np,1,'single');  % distance of point to the closest center 
BoP = zeros(np,1,'uint32');  % the balls/cover sets the points belong
nb = 0;             % number of sets generated

%% Same size cover sets everywhere
BallRad = inputs.BallRad;
PatchDiamMax = inputs.PatchDiam;
% Partition the point cloud into cubes for quick neighbor search
[partition,CC] = cubical_partition(P,BallRad);

% Generate the balls
Radius = BallRad^2;
MaxDist = PatchDiamMax^2;
RandPerm = uint32(randperm(np)); % random permutation of points,
% results in different covers for same input
for i = 1:np
    if NotExa(RandPerm(i))
        Q = RandPerm(i);
        points = partition(CC(Q,1)-1:CC(Q,1)+1,CC(Q,2)-1:CC(Q,2)+1,CC(Q,3)-1:CC(Q,3)+1);
        points = vertcat(points{:});
        V = [P(points,1)-P(Q,1) P(points,2)-P(Q,2) P(points,3)-P(Q,3)];
        dist = sum(V.*V,2);
        J = dist < Radius & Vessels(points) == Vessels(Q);
        I = points(J);
        d = dist(J);
        J = (dist < MaxDist);
        NotExa(points(J)) = false;
        nb = nb+1;
        Ball{nb} = I;
        Cen(nb) = Q;
        D = Dist(I);
        L = d < D;
        I = I(L);
        Dist(I) = d(L);
        BoP(I) = nb;
    end
end
Ball = Ball(1:nb,:);
Cen = Cen(1:nb);


%% Cover sets
% Number of points in each ball and index of each point in its ball
Num = zeros(nb,1,'uint32');
Ind = zeros(np,1,'uint32');
for i = 1:np
    if BoP(i) > 0
        Num(BoP(i)) = Num(BoP(i))+1;
        Ind(i) = Num(BoP(i));
    end
end

% Initialization of the "Bal"
Bal = cell(nb,1);
for i = 1:nb
    Bal{i} = zeros(Num(i),1,'uint32');
end

% Define the "Bal"
for i = 1:np
    if BoP(i) > 0
        Bal{BoP(i),1}(Ind(i)) = i;
    end
end

%% Neighbors
% Define neighbors. Sets A and B are neighbors if the large ball of A 
% contains points of B. Notice that this is not a symmetric relation.
Nei = cell(nb,1);
Fal = false(nb,1);
for i = 1:nb
    B = Ball{i};        % the points in the big ball of cover set "i" 
    I = (BoP(B) ~= i);  
    N = B(I);           % the points of B not in the cover set "i"
    N = BoP(N);
    
    % select the unique elements of N:
    n = length(N);
    if n > 2
        Include = true(n,1);
        for j = 1:n
            if ~Fal(N(j))
                Fal(N(j)) = true;
            else
                Include(j) = false;
            end
        end
        Fal(N) = false;
        N = N(Include);
    elseif n == 2
        if N(1) == N(2)
            N = N(1);
        end
    end
    
    Nei{i} = uint32(N);
end

% Make the relation symmetric by adding, if needed, A as B's neighbor 
% in the case B is A's neighbor
for i = 1:nb
    N = Nei{i};
    for j = 1:length(N)
        K = (Nei{N(j)} == i);
        if ~any(K)
            Nei{N(j)} = uint32([Nei{N(j)}; i]);
        end
    end
end

% Define output
clear cover
cover.ball = Bal;
cover.center = Cen;
cover.neighbor = Nei;
cover.BallOfPoint = BoP;
