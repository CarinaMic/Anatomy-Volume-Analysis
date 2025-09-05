%% Vertex-to-Nearest-Neighbour Distance

% Input:    pelvisNum: Numeric identifier used only for logging
%           shrink: Struct to which results will be added/updated
%           verticesRef: Reference point set for nearest-neighbour queries
%           vertices: Query vertices (typically the mesh vertices to evaluate)
%           faces: Triangle connectivity indexing into `vertices`

% Output:   shrink: struct with additional fields under `shrink.comp`

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shapeVertexNeighbour(shrink, pelvisNum, verticesRef, vertices, faces)

% Vertex-to-Nearest-Neighbour calculates the distance from each vertex to its
% nearest neighbor in the vertices list.
% Using a KD-tree searcher for improved performance
% Minimum distances from each vertex to the nearest vertex in reference

% Create a KD-tree searcher for verticesB
MdlB = KDTreeSearcher(verticesRef);

% Find the nearest neighbor in verticesB for each vertex in verticesA
idx = knnsearch(MdlB,vertices);

% Calculate the distances based on the nearest neighbor indices
minDistances = sqrt(sum((vertices - verticesRef(idx,:)).^2, 2));

% Pelvis defect with vertex-to-nearest-neighbour color
% Retrieve nearVertex for each vertex of each face
nearVertexMatrix = minDistances(faces);
% Calculate mean for each face
shrink.comp.nearVertexFace = mean(nearVertexMatrix,2);
% Normalisation
lowerBound = min(shrink.comp.nearVertexFace);
upperBound = max(shrink.comp.nearVertexFace);
% Normalize data within this range
normalData = (shrink.comp.nearVertexFace - lowerBound) / (upperBound - lowerBound);
normalData(normalData < 0) = 0;   % Clamp values below the range to 0
normalData(normalData > 1) = 1;   % Clamp values above the range to 1
% Apply colormap (recommended: viridis)
% 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
colourNum = 32768;
colourMap = viridis(colourNum);
colourIdx = round(normalData * (colourNum-1)) + 1;
rgbColour = colourMap(colourIdx,:);

shrink.comp.nearVertexFaceNorm = normalData;
shrink.comp.nearVertexFaceColour = rgbColour;
shrink.comp.nearVertex = minDistances;
shrink.comp.nearVertexMax = max(minDistances);
shrink.comp.nearVertexMean = mean(minDistances);
shrink.comp.nearVertexStd = std(minDistances);

disp(['Vertex-to-Nearest-Neighbour: pelvis ', num2str(pelvisNum)]);

end