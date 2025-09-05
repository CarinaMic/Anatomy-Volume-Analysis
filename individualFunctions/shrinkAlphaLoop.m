%% Shrink Wrap with alphaShape, more iterations for isclosed (without user-input)

% Input:    pelvisNum: Numeric identifier used only for logging
%           type: 'all','acentre','acetabulum'
%           numLoops: Number of closed-mesh attempts to perform. For each
%                     attempt, the alpha radius is increased until closed
%                     or until maxIterations is reached.
%           aRadius: Initial alpha radius %
%           numData: Number of rows in the original defect vertex array
%           importData:  Input geometry (vertices,comVertices,faces,acentre)
%           addPoints: Extra points appended to the point set

% Output:   shrink: Results struct with dynamic field named by 'type'

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shrinkAlphaLoop(pelvisNum, type, numLoops, aRadius, numData, importData, addPoints)

if isfield(importData, 'comVertices') && isfield(importData, 'comFaces')
    verticesImport = importData.comVertices;
else
    verticesImport = importData.vertices;
end
if nargin < 7
    acentre = importData.acentre;
    allVerticesPoints = [verticesImport; acentre];
else
    allVerticesPoints = [verticesImport; addPoints];
end

% Variables to store results from the three attempts
shrink.(type).successfulRadii = zeros(1,numLoops);
shrink.(type).alphaShapes = cell(1,numLoops);
shrink.(type).shrinkFacesList = cell(1,numLoops);

maxIterations = 150;
for attempt = 1:numLoops
    % Start with an initial small alpha radius and increase iteratively
    isClosed = false;
    iterationCount = 0;
    currentRadius = aRadius;
    while ~isClosed && iterationCount < maxIterations
        try
            % Create alphaShape and get the boundary facets
            shrinkWrap = alphaShape(allVerticesPoints, currentRadius);
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
                % Store successful aRadius, alphaShape, and shrinkFaces
                shrink.(type).successfulRadii(attempt) = currentRadius;
                shrink.(type).alphaShapes{attempt} = shrinkWrap;
                shrink.(type).shrinkFacesList{attempt} = shrinkFaces;
                break;
            end
        catch
            % If an error occurs, increase aRadius and try again
            isClosed = false;  % Ensure that the loop continues
        end
        % Increment the alpha radius
        % radius 1 may give no alphaShape
        currentRadius = currentRadius + 1;
        iterationCount = iterationCount + 1;
    end
    % Increment the starting radius for the next attempt
    aRadius = currentRadius + 1;
end

% Store the resulting alpha radius and alpha shape
shrink.(type).radius = shrink.(type).successfulRadii(1);
shrink.(type).alpha = shrink.(type).alphaShapes{1};
shrink.(type).faces = shrink.(type).shrinkFacesList{1};
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
for i = 1:size(shrink.(type).shrinkFacesList{1}, 1)
    % Get the vertices of the face
    v1 = allVerticesPoints(shrink.(type).shrinkFacesList{1}(i, 1), :);
    v2 = allVerticesPoints(shrink.(type).shrinkFacesList{1}(i, 2), :);
    v3 = allVerticesPoints(shrink.(type).shrinkFacesList{1}(i, 3), :);
    % Calculate the signed volume of the tetrahedron
    tetraVolume = dot(v1, cross(v2, v3)) / 6;
    % Add to total volume
    volumeAlpha = volumeAlpha + tetraVolume;
end
% The volume should be positive
shrink.(type).volume = abs(volumeAlpha);

disp(['shrink-wrapped mesh (alphaShape) iterations calculated: pelvis defect ', num2str(pelvisNum)])

end