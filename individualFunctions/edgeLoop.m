%% Identify edge / boundary loops

% Input:    pelvisNum: Numeric identifier used only for logging
%           importData: Struct containing precomputed boundary edges

% Output:   edge: Struct with boundary-loop information

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [edge] = edgeLoop(pelvisNum,importData)

edges = importData.verticesNrCombi;

% Initializing storage cells
verticesLoop = {};
comPoint = {};
i = 1; % loop counter
while ~isempty(edges)
    startPoint = edges(1,1);
    endPoint = edges(1,2);
    loopVertices = edges(1,:);
    edges(1,:) = []; % remove the starting edge

    d = endPoint;
    comPoint{i} = startPoint;
    while d ~= startPoint
        [row,~] = find(edges == d);
        if isempty(row)
            break; % if not found, break from loop
        end

        loopVertices = [loopVertices; edges(row(1), :)];
        comPoint{i} = [comPoint{i}; d];

        if edges(row(1), 1) == d    % new d
            d = edges(row(1), 2);
        else
            d = edges(row(1), 1);
        end

        edges(row(1),:) = []; % remove the used edge
    end

    verticesLoop{i} = loopVertices;
    i = i + 1;
end

% Collect unique vertices from the loops
edgeLoopsUni = cellfun(@(x) unique(x), verticesLoop, 'UniformOutput', false);

% Acetabulum edgeLoop
% Assumption: acetabulum hole with the most vertices (max length)
cellLength = cellfun(@length, edgeLoopsUni);
[~,loopNrAcentre] = max(cellLength);
% EdgeLoop Nr without acetabulum hole
loopNrHoles(:,1) = 1:size(edgeLoopsUni,2);
loopNrHoles(loopNrHoles == loopNrAcentre) = [];
if isempty(loopNrHoles)
    % Only loopNrAcentre exists
    reorderedLoops = edgeLoopsUni{loopNrAcentre};
    reorderedVerticesLoop = verticesLoop{loopNrAcentre};
    reorderedComPoint = comPoint{loopNrAcentre};
else
    % Sort loopNrHoles by their size in descending order
    [~,sortedIndices] = sort(cellLength(loopNrHoles), 'descend');
    sortedLoopNrHoles = loopNrHoles(sortedIndices);
    % Reorder edge loops based on the sorted loopNrAcentre and loopNrHoles
    reorderedLoops = [edgeLoopsUni{loopNrAcentre}, edgeLoopsUni(sortedLoopNrHoles)];
    reorderedVerticesLoop = [verticesLoop{loopNrAcentre}, verticesLoop(sortedLoopNrHoles)];
    reorderedComPoint = [comPoint{loopNrAcentre}, comPoint(sortedLoopNrHoles)];
end

% Store results in the object
edge.edgeLoops = reorderedLoops';
edge.loopsVerticesNrCombi = reorderedVerticesLoop';
edge.comPoint = reorderedComPoint';

disp(['edge loops calculated: pelvis defect ',num2str(pelvisNum)])

end