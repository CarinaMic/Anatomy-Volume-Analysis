%% Identify edge / boundary vertices

% Input:    pelvisNum: Numeric identifier used only for logging
%           importData  Struct with triangulated surface (vertices,faces)

% Output:   edge: Struct with boundary information

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [edge] = edging(pelvisNum,importData)

% Identify boundary vertices for calculating new meshes.
% An edge is defined as a triangular side that appears only once in the list of faces.

faces = importData.faces;
vertices = importData.vertices;

% Extract all edges from the faces
edges = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
sortedEdges = sort(edges, 2);

% Identify unique edges and their counts
[uniqueEdges, ~, edgeIndices] = unique(sortedEdges, 'rows');
edgeCounts = accumarray(edgeIndices, 1);

% Boundary edges appear only once
isBoundaryEdge = edgeCounts == 1;
boundaryEdges = uniqueEdges(isBoundaryEdge, :);

% Store results in the object
edge.verticesNrCombi = boundaryEdges;
edge.verticesNr = unique(boundaryEdges(:)); % Flatten the array and then find unique vertices
edge.vertices = vertices(edge.verticesNr,:);

disp(['edge calculated: pelvis defect ',num2str(pelvisNum)])

end