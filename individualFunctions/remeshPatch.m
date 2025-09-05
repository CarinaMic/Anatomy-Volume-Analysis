%% Combine mesh with patches

% Input:    pelvisNum: Numeric identifier used only for logging
%           patches: Struct with fields used under patches.(type):
%                    - numComponents : number of connected components in the base mesh
%                    - flip{i,1}     : user flag per loop (0=keep, 1=flip patch face order)
%           importData: Struct with triangulated mesh (vertices,faces)
%           slavesFaces: cell array, each cell i contains a Kx3 face list of the i-th patch
%           slavesVertices: cell array, each cell i contains a Px3 vertex array of the i-th patch
%           type: 'all','acentre','acetabulum'

% Output:   patches: updated struct under patches.(type) with:
%                    - comVertices: combined vertex list [base + patches] 
%                    - comFaces: combined face list 
%                    - comNormals: per-face normals (unit vectors) of comFaces
%                    - comCentreFaces: per-face centroids of comFaces

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function patches = remeshPatch(patches,pelvisNum,importData,slavesFaces,slavesVertices,type)

% Mesh master: imported defect data
currentVertices = importData.vertices;
importFaces = importData.faces;

% Mesh: normal orientation to the center -> defines the outside of the mesh
% -> Flip in- and outside of the mesh -> flip the normal orientation
% The orientation of the normal is defined by the ordering of the vertices in a triangle.
% The ordering follows the right-hand rule.
currentFaces = importFaces(:, [3, 2, 1]);

% Refined patch: new vertices added at the end
% number of edge loops; without acetabulum patch
for i = (patches.(type).numComponents+1):size(slavesFaces,1)
    % Combine faces
    mapSlaveFaces = slavesFaces{i};
    if patches.(type).flip{i,1} == 1
        mapSlaveFaces = slavesFaces{i}(:, [3, 2, 1]);
    end
    % Determine the offset for face indices due to appending new vertices
    offset = size(currentVertices, 1);
    % Adjust the face indices of refined patches
    adjustedFaces = mapSlaveFaces + offset;
    % Add the mapped slaveFaces to the combined mesh
    currentFaces = [currentFaces; adjustedFaces];

    % Combine vertices
    currentVertices = [currentVertices; slavesVertices{i}];
end
% Identify unique vertices in refinedVertices that are not in importData.vertices
[currentVertices, ~, idx] = unique(currentVertices, 'rows', 'stable');
% Adjust refinedFaces based on the removal of duplicates
currentFaces = idx(currentFaces);

% Store results in the object
patches.(type).comVertices = currentVertices;
patches.(type).comFaces = currentFaces;
% Compute the face normals
V1 = patches.(type).comVertices(patches.(type).comFaces(:, 1), :);
V2 = patches.(type).comVertices(patches.(type).comFaces(:, 2), :);
V3 = patches.(type).comVertices(patches.(type).comFaces(:, 3), :);
normals = cross(V2 - V1, V3 - V1, 2);
% Normalize the normals
patches.(type).comNormals = normals ./ vecnorm(normals, 2, 2);
% Calculate the centroid of the triangular faces in the mesh
faceVertices = cat(3, V1, V2, V3);
patches.(type).comCentreFaces = mean(faceVertices, 3);

disp(['combined mesh with patches: pelvis defect ',num2str(pelvisNum)])

end