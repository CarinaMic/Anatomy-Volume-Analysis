%% Close holes and triangulate (patch) + refine patches

% Input:    edge: Struct with boundary-loop info (from edgeLoop)
%           pelvisNum: Numeric identifier used only for logging
%           importData: Struct with triangulated mesh (vertices,faces)
%           varargin: Optional arguments (order-sensitive)
%                 - type: 'all' (default) | 'acentre' | 'acetabulum'
%                 - extra_data: 'acentre': coordinate to be appended before triangulation
%                               'acetabulum': coordinates (loop vertices) appended before triangulation

% Output:   patches: Struct with results under field patches.(type), where 'type'
%               is one of {'all','acentre','acetabulum'}:

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [patches] = meshPatch(edge,pelvisNum,importData,varargin)

vertices = importData.vertices;
faces = importData.faces;

% Default type is 'all'
type = 'all';
extra_data = [];
% Parse varargin for optional extra_data and type
if ~isempty(varargin)
    if ischar(varargin{end})
        type = varargin{end};
        if length(varargin) > 1
            extra_data = varargin{1};
        end
    else
        % extra_data - acentre: coordinates of landmark acentre;
        % acetabulum: coordinates of acetabulum loop (vertices)
        extra_data = varargin{1};
        if isempty(extra_data)
            type = 'all';
        elseif size(extra_data,2) == 3
            type = 'acetabulum';
        else
            type = 'acentre';
        end
    end
end

% Detect seperate components of mesh
% Create adjacence matrix
edges = [faces(:, [1,2]); faces(:, [2,3]); faces(:, [3,1])];
edges = unique(sort(edges, 2), 'rows');
% Create graph
G = graph(edges(:,1), edges(:,2));
% Find connected components using graph theory
componentLabels = conncomp(G, 'Type', 'weak');
numComponents = max(componentLabels);
patches.(type).numComponents = numComponents;

% Determine the number of loops to process based on the type
if any(strcmp(type, {'acentre', 'acetabulum'}))
    numLoops = 1; % Only process the first edge loop (acetabulum loop)
    numComponents = 0; % reset
else
    numLoops = size(edge.comPoint, 1);
end
data_sorted = cell(numLoops, 1);
Tri = cell(numLoops, 1);

% Mesh edge loops
% Matlab File Exchange: discrete_contour_mesh_patch
% https://de.mathworks.com/matlabcentral/fileexchange/78901-discrete-contour-mesh-patch-2d-3d?s_tid=srchtitle
% Mostly efficient for quasi flat contours, cause doesn't yet prevent from self intersecting triangles.
for i = 1:numLoops
    data_sorted{i,1} = vertices(edge.comPoint{i,1}(:),:); % comPoint: sorted edge loop points
    if strcmp(type, 'acentre') || strcmp(type, 'acetabulum') % only one loop
        data_sorted{i,1} = [data_sorted{i,1}; extra_data]; % extra_data = coordinates
    end
    [~,Tri{i,1},~] = discrete_contour_mesh_patch(data_sorted{i,1},'sorted');
end

% Store data to obj
patches.(type).faces = Tri;
patches.(type).vertices = data_sorted;
% Compute the face normals
for i = 1:numLoops
    V1 = patches.(type).vertices{i,1}(patches.(type).faces{i,1}(:,1),:);
    V2 = patches.(type).vertices{i,1}(patches.(type).faces{i,1}(:,2),:);
    V3 = patches.(type).vertices{i,1}(patches.(type).faces{i,1}(:,3),:);
    normals = cross(V2 - V1, V3 - V1, 2);
    % Normalize the normals
    patches.(type).normals{i,1} = normals ./ vecnorm(normals, 2, 2);
    % Calculate the centroid of the triangular faces in the mesh
    faceVertices = cat(3, V1, V2, V3);
    patches.(type).centreFaces{i,1} = mean(faceVertices, 3);
end

% Refine patch using loop subdivision
if ~any(strcmp(type, {'acentre', 'acetabulum'})) % not acetabulum patch
    desired_length = 0.5; % Desired edge length in mm
    tolerance = 0.05;
    maxIterations = 3;
    for i = (numComponents+1):numLoops
        patch.vertices = data_sorted{i};
        patch.faces = Tri{i};
        counter = 0;
        while counter < maxIterations
            % Loop subdivision
            % Map of edges
            patch.edges = unique(sort([patch.faces(:, [1,2]); patch.faces(:, [2,3]); patch.faces(:, [3,1])], 2), 'rows');
            edgeMap = zeros(size(patch.vertices,1), size(patch.vertices,1));
            edgeMidpointIndex = size(patch.vertices,1) + (1:size(patch.edges,1));
            newVertices = (patch.vertices(patch.edges(:,1),:) + patch.vertices(patch.edges(:,2),:)) / 2;
            patch.vertices = [patch.vertices; newVertices];
            % Update edge map
            for j = 1:size(patch.edges,1)
                edgeMap(patch.edges(j,1), patch.edges(j,2)) = edgeMidpointIndex(j);
                edgeMap(patch.edges(j,2), patch.edges(j,1)) = edgeMidpointIndex(j);
            end
            % Create new faces
            oldFaces = patch.faces;
            patch.faces = zeros(size(oldFaces,1) * 4, 3);
            for j = 1:size(oldFaces,1)
                v1 = oldFaces(j,1);
                v2 = oldFaces(j,2);
                v3 = oldFaces(j,3);
                v12 = edgeMap(v1, v2);
                v23 = edgeMap(v2, v3);
                v31 = edgeMap(v3, v1);
                patch.faces(4*j-3,:) = [v1, v12, v31];
                patch.faces(4*j-2,:) = [v2, v23, v12];
                patch.faces(4*j-1,:) = [v3, v31, v23];
                patch.faces(4*j,:) = [v12, v23, v31];
            end
            % Extract edges from the faces
            patch.edges = [patch.faces(:,1), patch.faces(:,2);
                patch.faces(:,2), patch.faces(:,3);
                patch.faces(:,3), patch.faces(:,1)];
            patch.edges = unique(sort(patch.edges, 2), 'rows', 'stable');  % Sort and remove duplicate edges
            % Compute edge lengths from refined mesh
            v1 = patch.vertices(patch.edges(:,1), :);  % Start points of edge
            v2 = patch.vertices(patch.edges(:,2), :);  % End points of edge
            edge_lengths = sqrt(sum((v2 - v1).^2, 2));  % Euclidean distances
            % Test if the majority of edges are close to the desired length
            closeDesiredEdge = abs(edge_lengths - desired_length) < tolerance;
            if mean(closeDesiredEdge) > 0.95  % If more than 95% of edges are close
                break;  % Exit while loop
            end
            counter = counter + 1;  % Increment the counter at the end of the loop
        end

        % Store refined patch data to obj
        patches.(type).refinedVertices{i,1} = patch.vertices;
        patches.(type).refinedFaces{i,1} = patch.faces;
        % Compute the face normals
        V1 = patches.(type).refinedVertices{i,1}(patches.(type).refinedFaces{i,1}(:,1),:);
        V2 = patches.(type).refinedVertices{i,1}(patches.(type).refinedFaces{i,1}(:,2),:);
        V3 = patches.(type).refinedVertices{i,1}(patches.(type).refinedFaces{i,1}(:,3),:);
        normals = cross(V2 - V1, V3 - V1, 2);
        % Normalize the normals
        patches.(type).refinedNormals{i,1} = normals ./ vecnorm(normals, 2, 2);
        % Calculate the centroid of the triangular faces in the mesh
        faceVertices = cat(3, V1, V2, V3);
        patches.(type).refinedCentreFaces{i,1} = mean(faceVertices, 3);
    end
end

disp(['meshed edge loops (patch - ', type, '): pelvis defect ', num2str(pelvisNum)]);

end