%% Shrink Wrap with alphaShape (without user-input)

% Input:    pelvisNum: Numeric identifier used only for logging
%           type: 'all','acentre','acetabulum'
%           aRadius: Initial alpha radius %
%           numData: Number of rows in the original defect vertex array
%           importData:  Input geometry (vertices,comVertices,faces,acentre)
%           addPoints: Extra points appended to the point set

% Output:   shrink: Results struct with dynamic field named by 'type'

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shrinkAlpha(pelvisNum, type, aRadius, numData, importData, addPoints)

if isfield(importData, 'comVertices') && isfield(importData, 'comFaces')
    verticesImport = importData.comVertices;
else
    verticesImport = importData.vertices;
end
if nargin < 6
    acentre = importData.acentre;
    allVerticesPoints = [verticesImport; acentre];
else
    allVerticesPoints = [verticesImport; addPoints];
end

% Start with an initial small alpha radius and increase iteratively
isClosed = false;
maxIterations = 150;
iterationCount = 0;
while ~isClosed && iterationCount < maxIterations
    try
        % Create alphaShape and get the boundary facets
        shrinkWrap = alphaShape(allVerticesPoints, aRadius);
        shrinkFaces = boundaryFacets(shrinkWrap);

        % Check if the mesh is closed
        edgesList = [shrinkFaces(:,1) shrinkFaces(:,2);
            shrinkFaces(:,2) shrinkFaces(:,3);
            shrinkFaces(:,3) shrinkFaces(:,1)];
        edgesList = sort(edgesList, 2);  % Sort each edge
        [~, ~, edgeOccurrences] = unique(edgesList, 'rows', 'stable');
        edgeCounts = accumarray(edgeOccurrences, 1);

        % Check if any edge is not shared by exactly two triangles
        isClosed = all(edgeCounts == 2);
        % If the mesh is closed, break out of the loop (no increase of alpha radius)
        if isClosed
            break;
        end
    catch
        % If an error occurs, increase aRadius and try again
        isClosed = false;  % Ensure that the loop continues
    end
    % Increment the alpha radius
    % radius 1 may give no alphaShape
    aRadius = aRadius + 1;
    iterationCount = iterationCount + 1;
end

% Store the resulting alpha radius and alpha shape
shrink.(type).radius = aRadius;
shrink.(type).alpha = shrinkWrap;
shrink.(type).faces = shrinkFaces;
shrink.(type).allVerticesPoints = allVerticesPoints; % all vertices

% Logical index for all used vertices/vertices in the mesh
% original data (vertices): without patches, only original mesh
uniVertices = unique(shrinkFaces(:));
usedDefectVerticesImport = uniVertices(uniVertices <= numData);
shrink.(type).usedDefectVerticesMask = false(numData,1); % numData: size of original vertices data
shrink.(type).usedDefectVerticesMask(usedDefectVerticesImport) = true;
shrink.(type).usedDefectVerticesCount = length(usedDefectVerticesImport) / numData;
shrink.(type).usedDefectVertices = verticesImport(usedDefectVerticesImport, :);
% All vertices/points
shrink.(type).usedVerticesPointsMask = false(size(allVerticesPoints, 1), 1);
shrink.(type).usedVerticesPointsMask(uniVertices) = true;
shrink.(type).usedVerticesPoints = allVerticesPoints(shrink.(type).usedVerticesPointsMask, :);
% New faces structure: renumbering faces (faces structure with new numbering of usedVerticesPoints)
vertexIndexMapping = zeros(size(shrink.(type).usedVerticesPointsMask));
vertexIndexMapping(shrink.(type).usedVerticesPointsMask) = 1:nnz(shrink.(type).usedVerticesPointsMask);
remappedFaces = vertexIndexMapping(shrink.(type).faces);
shrink.(type).usedFaces = remappedFaces;

% Volume (divergence theorem)
volumeAlpha = 0; % Initialize volume
% Loop through each face
for i = 1:size(shrinkFaces, 1)
    % Get the vertices of the face
    v1 = allVerticesPoints(shrinkFaces(i, 1), :);
    v2 = allVerticesPoints(shrinkFaces(i, 2), :);
    v3 = allVerticesPoints(shrinkFaces(i, 3), :);
    % Calculate the signed volume of the tetrahedron
    tetraVolume = dot(v1, cross(v2, v3)) / 6;
    % Add to total volume
    volumeAlpha = volumeAlpha + tetraVolume;
end
% The volume should be positive
shrink.(type).volume = abs(volumeAlpha);

disp(['shrink-wrapped mesh (alphaShape) calculated: pelvis defect ', num2str(pelvisNum)])

end