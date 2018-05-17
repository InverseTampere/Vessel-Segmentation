function [Partition,CubeCoord,Info,Cubes] = cubical_partition(P,EL,NE)

% Partitions the input point cloud into cubes.
%
% Inputs:
% P             Point cloud, (n_points x 3)-matrix
% EL            Length of the cube edges
% NE            Number of empty edge layers
%
% Outputs:              
% Partition     Point cloud partitioned into cubical cells,
%                   (nx x ny x nz)-cell, where nx,ny,nz are the number
%                   of cubes in x,y,z-directions, respectively
% CC            (n_points x 3)-matrix whose rows are the cube coordinates 
%                   of each point: x,y,z-coordinates
% Info          The minimum coordinate values and number of cubes in each
%                   coordinate direction

if nargin == 2
    NE = 3;
end

% The vertices of the big cube containing P
Min = double(min(P));
Max = double(max(P));

% Number of cubes with edge length "EdgeLength" in the sides 
% of the big cube
N = double(ceil((Max-Min)/EL)+2*NE+1);
while 8*N(1)*N(2)*N(3) > 4e9
    EL = 1.1*EL;
    N = double(ceil((Max-Min)/EL)+2*NE+1);
end
Info = [Min N EL NE];

% Calculates the cube-coordinates of the points
CubeCoord = floor([P(:,1)-Min(1) P(:,2)-Min(2) P(:,3)-Min(3)]/EL)+NE+1;

% Sorts the points according a lexicographical order
LexOrd = [CubeCoord(:,1) CubeCoord(:,2)-1 CubeCoord(:,3)-1]*[1 N(1) N(1)*N(2)]';
CubeCoord = uint16(CubeCoord);
[LexOrd,SortOrd] = sort(LexOrd);
SortOrd = uint32(SortOrd);
LexOrd = uint32(LexOrd);

if nargout <= 3
    % Define "Partition"
    Partition = cell(N(1),N(2),N(3));
    np = size(P,1);     % number of points
    p = 1;              % The index of the point under comparison
    while p <= np
        t = 1;
        while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        q = SortOrd(p);
        Partition{CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)} = SortOrd(p:p+t-1);
        p = p+t;
    end
    
else
    nc = size(unique(LexOrd),1);
    
    % Define "Partition"
    Cubes = zeros(N(1),N(2),N(3),'uint32');
    Partition = cell(nc,1);
    np = size(P,1);     % number of points
    p = 1;              % The index of the point under comparison
    c = 0;
    while p <= np
        t = 1;
        while (p+t <= np) && (LexOrd(p) == LexOrd(p+t))
            t = t+1;
        end
        q = SortOrd(p);
        c = c+1;
        Partition{c,1} = SortOrd(p:p+t-1);
        Cubes(CubeCoord(q,1),CubeCoord(q,2),CubeCoord(q,3)) = c;
        p = p+t;
    end
end