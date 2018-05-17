function NewPoints = combine_vessels(Data,P,VesselOfPoint,Tips,...
    TipPoints,MaxDist,MaxAngle)

% ---------------------------------------------------------------------
% COMBINE_VESSELS.M     Generate new points by combining vessels that are 
%                           close to each other and have almost parallel 
%                           directions at their tips.
%
% Version 1.0.0
% Updated       17 May 2018
% Copyright (C) 2017-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Input: 
% Data          Photoacoustic image, 3d-array of normalized intensity values
% P             3D point cloud of the non-zero data points
% VesselOfPoint Vector giving the vessel of each point in the point cloud P
% Tips          Matrix containing information about the tips
% TipPoints     Cell-array containing the points forming the tips
% MaxDist       Parameter defining the maximum distance between tips that
%                   can be joined
% MaxAngle      Parameter defining the maximum angle between tips that
%                   can be joined
%
% Output:
% NewPoints     Cell-array, contains the new points that close gabs between
%                   vessel tips and information (distance, angle)
% ---------------------------------------------------------------------

% Generate new points by combining vessels that are close to each
% other and have almost parallel directions at their tips.

% First modify data so that only those points segmented so far are data
np = size(P,1);
for i = 1:np
    if VesselOfPoint(i) == 0
        Data(P(i,1),P(i,2),P(i,3)) = 0;
    end
end

% Partition the tips for neighbor search
[partition,CC] = cubical_partition(Tips(:,1:3),MaxDist);

% Search the neighboring tips
nt = size(TipPoints,1);
Nei = cell(max(Tips(:,7)),1);
NewPoints = cell(nt,2);
nadded = 0;
for i = 1:nt
    tips = partition(CC(i,1)-1:CC(i,1)+1,CC(i,2)-1:CC(i,2)+1,CC(i,3)-1:CC(i,3)+1);
    tips = vertcat(tips{:}); % potential tips
    I = tips ~= i; % other tips than the current one
    J = Tips(tips,7) ~= Tips(i,7); % tips from other vessels
    I = I&J;
    tips = tips(I); % the potential tips
    if ~isempty(tips)
        % Compute angles and select suitable tips
        V = [Tips(tips,1)-Tips(i,1) Tips(tips,2)-Tips(i,2) Tips(tips,3)-Tips(i,3)];
        L = sqrt(sum(V.*V,2));
        W = [V(:,1)./L V(:,2)./L V(:,3)./L];
        A1 = 180/pi*acos(W*Tips(i,4:6)');
        A2 = 180/pi*acos(sum(-W.*Tips(tips,4:6),2));
        I = A1 <= MaxAngle & A2 <= MaxAngle & L <= MaxDist;

        if any(I)
            % Select the closest tip
            n = length(A1);
            ind = (1:1:n)';
            L = L(I);
            A1 = A1(I);
            A2 = A2(I);
            ind = ind(I);
            [L,I] = min(L);
            A1 = A1(I);
            A2 = A2(I);
            ind = ind(I);
            t = tips(ind); % the selected tip to be joined to the current tip
            comp1 = Tips(i,7); % the component of the current tip
            comp2 = Tips(t,7); % the component of the selected tip
            Join = ~any(Nei{comp1} == comp2); % Join if the components are not already joined
            if Join && L <= MaxDist
                Bot = TipPoints{i};
                Top = TipPoints{t};
                Q = generate_points(Data,P,Bot,Top);
                Nei{comp1} = [Nei{comp1}; comp2]; % Join the components
                Nei{comp2} = [Nei{comp2}; comp1]; % Join the components
                nadded = nadded+1;
                NewPoints{nadded,1} = Q;
                NewPoints{nadded,2} = [L max(A1,A2)];
            end
        end
    end
end
NewPoints = NewPoints(1:nadded,:);

end % End of main function


function Q = generate_points(Data,P,Bot,Top)

% Generate points on the lines connecting every point of the Top to every
% point of the Bot
n = length(Bot);
m = length(Top);
Q = zeros(m*n*20,3);
t = 0;
for i = 1:n
    for j = 1:m
        V = P(Top(j),:)-P(Bot(i),:);
        L = ceil(norm(V));
        l = L-1;
        for k = 1:l
            t = t+1;
            Q(t,:) = P(Bot(i),:)+k/L*V;
        end
    end
end
Q = Q(1:t,:);
Q = round(Q); % Only whole number values as in grid points
Q = unique(Q,'rows'); % Only unique points
n = size(Q,1);

% Keep only new points not in the original data
Keep = true(n,1);
for i = 1:n
    if Data(Q(i,1),Q(i,2),Q(i,3)) > 0
        Keep(i) = false;
    elseif Q(i,1) == 0 || Q(i,2) == 0 || Q(i,3) == 0
        Keep(i) = false;
    end
end
Q = Q(Keep,:);

end % End of subfunction

