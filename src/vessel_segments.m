function segment = vessel_segments(cover,Base,Forb,classification)

% ---------------------------------------------------------------------
% VESSEL_SEGMENTS.M     Segments the covered point cloud into branches/vessels.
%
% Version 1.0.0
% Latest update     17 May 2018
% Copyright (C) 2017-2018 Pasi Raumonen
% ---------------------------------------------------------------------
%
% Segments the vessel clusters into vessels and records their parent-child-relations. 
% Bifurcations are recognized by studying connectivity of a "study" region 
% moving along the vessel network. In case of multiple connected components 
% in "study", the components are classified as the continuation and branches.
%
% Inputs:
% cover         Cover sets
% Base          Base of the tree
% Forb          Cover sets not part of the tree
%
% Outputs:
% segment       Structure array containing the followin fields:
%   segments          Segments found, (n_seg x 1)-cell, each cell contains 
%                       a cell array the cover sets
%   ParentSegment     Parent segment of each segment, (n_seg x 1)-vector,
%                       equals to zero if no parent segment
%   ChildSegment      Children segments of each segment, (n_seg x 1)-cell

Nei = cover.neighbor;
nb = size(Nei,1);           % The number of cover sets
a = max([200000 nb/100]);   % Estimate for maximum number of segments
SBas = cell(a,1);           % The segment bases found
Segs = cell(a,1);           % The segments found
SPar = zeros(a,2,'uint32'); % The parent segment of each segment
SChi = cell(a,1);           % The children segments of each segment

% Initialize SChi
SChi{1} = zeros(5000,1,'uint32');
C = zeros(200,1);
for i = 2:a
    SChi{i} = C;
end
NChi = zeros(a,1);      % Number of child segments found for each segment

Fal = false(nb,1);      % Logical false-vector for cover sets
s = 1;                  % The index of the segment under expansion
b = size(Base,1);                  % The index of the latest found base

SBas(1:b) = Base;
Seg = cell(1000,1);    % The cover set layers in the current segment
Seg{1} = Base{1};

ForbAll = Fal;       % The forbidden sets
ForbAll(Forb) = true;
ForbAll(vertcat(Base{:})) = true;
Forb = ForbAll;      % The forbidden sets for the segment under expansion

Continue = true; % True as long as the component can be segmented further 
NewSeg = true;   % True if the first Cut for the current segment
nl = 1;          % The number of cover set layers currently in the segment

% Segmenting stops when there are no more segments to be found
while Continue && (b < nb)
    
    % Update the forbidden sets
    Forb(Seg{nl}) = true;
    
    % Define the study
    Cut = define_cut(Nei,Seg{nl},Forb,Fal);
    CutSize = length(Cut);
    
    if NewSeg
        NewSeg = false;
        ns = min(CutSize,6);
        ns = max(ns,3);
    end
    
    % Define the components of cut and study regions
    if CutSize > 0
        CutComps = cut_components(Nei,Cut,CutSize,Fal,Fal);
        nc = size(CutComps,1);
        if nc > 1
            [StudyComps,Bases,CompSize] = ...
                study_components(Nei,ns,Cut,CutComps,Forb,Fal,Fal);
            nc = length(CompSize);
        end
    else
        nc = 0;
    end
    
    % Classify study region components
    if nc == 1
        % One component, continue expansion of the current segment
        nl = nl+1;
        if size(Cut,2) > 1
            Seg{nl} = Cut';
        else
            Seg{nl} = Cut;
        end
    elseif nc > 1
        % Classify the components of the Study region
        Class = component_classification(CompSize,classification);
        
        for i = 1:nc
            if Class(i) == 1
                Base = Bases{i};
                ForbAll(Base) = true;
                Forb(StudyComps{i}) = true;
                J = Forb(Cut);
                Cut = Cut(~J);
                b = b+1;
                SBas{b} = Base;
                SPar(b,:) = [s nl];
                NChi(s) = NChi(s)+1;
                SChi{s}(NChi(s)) = b;
            end
        end
        
        % Define the new cut.
        % If the cut is empty, determine the new base
        if isempty(Cut)
            Segs{s} = Seg(1:nl);
            S = vertcat(Seg{1:nl});
            ForbAll(S) = true;

            if s < b
                s = s+1;
                Seg{1} = SBas{s};
                Forb = ForbAll;
                NewSeg = true;
                nl = 1;
            else
                Continue = false;
            end
        else
            if size(Cut,2) > 1
                Cut = Cut';
            end
            nl = nl+1;
            Seg{nl} = Cut;
        end
    
    else
        % If the study region has zero size, then the current segment is
        % complete and determine the base of the next segment
        Segs{s} = Seg(1:nl);
        S = vertcat(Seg{1:nl});
        ForbAll(S) = true;
        
        if s < b
            s = s+1;
            Seg{1} = SBas{s};
            Forb = ForbAll;
            NewSeg = true;
            nl = 1;
        else
            Continue = false;
        end
    end
end
Segs = Segs(1:b);
SPar = SPar(1:b,:);
schi = SChi(1:b);

% Define output
SChi = cell(b,1);
for i = 1:b
    if NChi(i) > 0
        SChi{i} = uint32(schi{i}(1:NChi(i)));
    else
        SChi{i} = zeros(0,1,'uint32');
    end
    S = Segs{i};
    for j = 1:size(S,1)
        S{j} = uint32(S{j});
    end
    Segs{i} = S;
end
clear Segment
segment.segments = Segs;
segment.ParentSegment = SPar;
segment.ChildSegment = SChi;

end % End of the main function


function Cut = define_cut(Nei,CutPre,Forb,Fal)

% Defines the "Cut" region
Cut = vertcat(Nei{CutPre});
Cut = unique_elements(Cut,Fal);
I = Forb(Cut);
Cut = Cut(~I);
end % End of function 


function [Components,CompSize] = cut_components(Nei,Cut,CutSize,Fal,False)

% Define the connected components of the Cut
if CutSize == 1
    % Cut is connected and therefore Study is also
    CompSize = 1;
    Components = cell(1,1);
    Components{1} = Cut;
elseif CutSize == 2
    I = Nei{Cut(1)} == Cut(2);
    if any(I)
        Components = cell(1,1);
        Components{1} = Cut;
        CompSize = 1;
    else
        Components = cell(2,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        CompSize = [1 1];
    end
elseif CutSize == 3
    I = Nei{Cut(1)} == Cut(2);
    J = Nei{Cut(1)} == Cut(3);
    K = Nei{Cut(2)} == Cut(3);
    if any(I)+any(J)+any(K) >= 2
        CompSize = 1;
        Components = cell(1,1);
        Components{1} = Cut;
    elseif any(I)
        Components = cell(2,1);
        Components{1} = Cut(1:2);
        Components{2} = Cut(3);
        CompSize = [2 1];
    elseif any(J)
        Components = cell(2,1);
        Components{1} = Cut([1 3]');
        Components{2} = Cut(2);
        CompSize = [2 1];
    elseif any(K)
        Components = cell(2,1);
        Components{1} = Cut(2:3);
        Components{2} = Cut(1);
        CompSize = [2 1];
    else
        CompSize = [1 1 1];
        Components = cell(3,1);
        Components{1} = Cut(1);
        Components{2} = Cut(2);
        Components{3} = Cut(3);
    end
else
    Components = cell(CutSize,1);
    CompSize = zeros(CutSize,1);
    Comp = zeros(CutSize,1);
    Fal(Cut) = true;
    nc = 0;      % number of components found
    m = Cut(1);
    i = 0;
    while i < CutSize
        Added = Nei{m};
        I = Fal(Added);
        Added = Added(I);
        a = length(Added);
        Comp(1) = m;
        Fal(m) = false;
        t = 1;
        while a > 0
            Comp(t+1:t+a) = Added;
            Fal(Added) = false;
            t = t+a;
            Ext = vertcat(Nei{Added});
            Ext = unique_elements(Ext,False);
            I = Fal(Ext);
            Added = Ext(I);
            a = length(Added);
        end
        i = i+t;
        nc = nc+1;
        Components{nc} = Comp(1:t);
        CompSize(nc) = t;
        if i < CutSize
            J = Fal(Cut);
            m = Cut(J);
            m = m(1);
        end
    end
    Components = Components(1:nc);
    CompSize = CompSize(1:nc);
end

end % End of function


function [Components,Bases,CompSize] = ...
    study_components(Nei,ns,Cut,CutComps,Forb,Fal,False)

% Define Study as a cell-array
Study = cell(ns,1);
StudySize = zeros(ns,1);
Study{1} = Cut;
StudySize(1) = length(Cut);
if ns >= 2
    N = Cut;
    i = 1;
    while i < ns
        Forb(N) = true;
        N = vertcat(Nei{N});
        N = unique_elements(N,Fal);
        I = Forb(N);
        N = N(~I);
        if ~isempty(N)
            i = i+1;
            Study{i} = N;
            StudySize(i) = length(N);
        else
            Study = Study(1:i);
            StudySize = StudySize(1:i);
            i = ns+1;
        end
    end
end

% Define study as a vector
studysize = sum(StudySize);
study = vertcat(Study{:});

% Determine the components of study
nc = size(CutComps,1);
i = 1; % index of cut component
j = 0; % number of elements attributed to components
k = 0; % number of study components
Fal(study) = true;
Components = cell(nc,1);
CompSize = zeros(nc,1);
Comp = zeros(studysize,1);
while i <= nc
    C = CutComps{i};
    while j < studysize
        a = length(C);
        Comp(1:a) = C;
        Fal(C) = false;
        if a > 1
            Add = unique_elements(vertcat(Nei{C}),False);
        else
            Add = Nei{C};
        end
        t = a;
        I = Fal(Add);
        Add = Add(I);
        a = length(Add);
        while a > 0
            Comp(t+1:t+a) = Add;
            Fal(Add) = false;
            t = t+a;
            Add = vertcat(Nei{Add});
            Add = unique_elements(Add,False);
            I = Fal(Add);
            Add = Add(I);
            a = length(Add);
        end
        j = j+t;
        k = k+1;
        Components{k} = Comp(1:t);
        CompSize(k) = t;
        if j < studysize
            C = zeros(0,1);
            while i < nc && isempty(C)
                i = i+1;
                C = CutComps{i};
                J = Fal(C);
                C = C(J);
            end
            if i == nc && isempty(C)
                j = studysize;
                i = nc+1;
            end
        else
            i = nc+1;
        end
    end
    Components = Components(1:k);
    CompSize = CompSize(1:k);
end

% Determine Bases
Bases = cell(k,1);
if k > 1
    Fal(study) = false;
    Fal(Study{1}) = true;
    for i = 1:k
        % Determine the size of the base of the components
        Set = unique_elements([Components{i}; Study{1}],False);
        False(Components{i}) = true;
        I = False(Set)&Fal(Set);
        False(Components{i}) = false;
        Set = Set(I);
        Bases{i} = Set;
    end
end

end % End of function


function Class = component_classification(CompSize,classification)

% Classifies study region components:
% Class(i) == 0 continuation
% Class(i) == 1 branch

nc = size(CompSize,1);
Class = ones(nc,1);
if classification == 1
    [s,I] = sort(CompSize,'descend');
    if s(1) > 1.5*s(2)
        Class(I(1)) = 0;
    end
else
    [~,I] = max(CompSize);
    Class(I) = 0;
end
end % End of function
