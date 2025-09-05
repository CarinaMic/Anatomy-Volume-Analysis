%% Identify clusters and outliers of the point cloud

% Input:    pelvisNum: Numeric identifier used only for logging
%           shrink: Structure contain fields from prior steps (density-filtered inliers)
%           refData: Reference pelvis mesh (vertices,faces)
%           epsilon: DBSCAN neighborhood radius
%           minPts: DBSCAN minimum points per core point
%           importData: For context/visualization only 

% Output:   shrink: struct with additional fields under shrink.inside

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = clusterPoints(shrink,pelvisNum,epsilon,minPts,importData,refData)

verticesMesh = importData.comVertices;
facesMesh = importData.comFaces;
refFaces = refData.faces;
% Vertices / Points
filteredInVerticesPoints = [shrink.inside.filteredInVertices; shrink.inside.filteredInPoints];
% Mark the sources of points: 1 for vertices, 2 for points
sourceLabels = [ones(size(shrink.inside.filteredInVertices, 1), 1); ...
    ones(size(shrink.inside.filteredInPoints, 1), 1) * 2];

% Cluster analysis with DBSCAN
tic
labels = dbscan(filteredInVerticesPoints, epsilon, minPts);
shrink.inside.executionTimeCluster = toc;
% Identify cluster with the most points
unique_labels = unique(labels);
unique_labels(unique_labels == -1) = []; % Remove the label for outliers
% Find the largest cluster
max_cluster_size = 0;
max_cluster_label = NaN; % assumption: at least one cluster
for i = 1:length(unique_labels)
    current_label = unique_labels(i);
    current_cluster_size = sum(labels == current_label);
    if current_cluster_size > max_cluster_size
        max_cluster_size = current_cluster_size;
        max_cluster_label = current_label;
    end
end

% Plot and allow the user to keep or discard clusters
if length(unique_labels) > 1
    beep % sound for user-input
    for i = 1:length(unique_labels)
        if unique_labels(i) ~= max_cluster_label

            % For reference
            if pelvisNum == 1
                labels(labels == unique_labels(i)) = max_cluster_label;
            else

                % Plot the main cluster and other clusters
                figure;
                hold on;
                % Pelvis defect remeshed (scaled and transformed)
                patch('Faces',facesMesh,...
                    'Vertices',verticesMesh,...
                    'FaceColor',[0.9 0.75 0.68], ...    % Face color
                    'FaceAlpha',1,...                   % Transparency of the faces
                    'EdgeColor','none',...    % Edge color TUMcolors.grey50
                    'EdgeAlpha',0.25);                  % Transparency of the edges
                light('Position', [1 1 5], 'Style', 'infinite');
                % Reference pelvis mesh inside
                filteredInVertices = NaN(max(shrink.inside.inVerticesIdx), 3);
                filteredInVertices(shrink.inside.inVerticesIdx, :) = shrink.inside.inVertices;
                patch('Faces',shrink.inside.filteredInFaces,...
                    'Vertices',filteredInVertices,... %%%
                    'FaceColor',[0.8 0.8 0.8], ...      % Face color
                    'FaceAlpha',1,...                   % Transparency of the faces
                    'EdgeColor',[128/255 128/255 128/255],...    % Edge color
                    'EdgeAlpha',0.25);                  % Transparency of the edges
                % Main cluster
                plot3(filteredInVerticesPoints(labels == max_cluster_label, 1), ...
                    filteredInVerticesPoints(labels == max_cluster_label, 2), ...
                    filteredInVerticesPoints(labels == max_cluster_label, 3), ...
                    '.', 'Color', [0 101/255 189/255], 'MarkerSize', 2);
                % Other cluster
                plot3(filteredInVerticesPoints(labels == unique_labels(i), 1), ...
                    filteredInVerticesPoints(labels == unique_labels(i), 2), ...
                    filteredInVerticesPoints(labels == unique_labels(i), 3), ...
                    'x', 'Color', [227/255 114/255 34/255], 'MarkerSize', 5);
                % Format and display properties
                title(['Cluster ', num2str(unique_labels(i))]);
                xlabel('X');
                ylabel('Y');
                zlabel('Z');
                legend('ReferencePelvis', 'Main Cluster', 'Other Cluster');
                daspect([1, 1, 1]); % Equal aspect ratio for the axes
                view(3);
                rotate3d on;        % Enable 3D rotation in the figure

                % Add toggle button and confirm button
                toggle = uicontrol('Style', 'togglebutton', 'String', 'Keep Cluster', ...
                    'Position', [20 60 100 40]);
                uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
                    'Position', [20 20 80 40], ...
                    'Callback', @(src,event)confirmCluster(unique_labels(i), toggle.Value));
                uiwait(gcf); % Wait for the user to respond before proceeding
            end
        end
    end
end

% Updated labels
shrink.inside.clusterLabels = labels;
shrink.inside.clusterMainLabel = max_cluster_label;
% Extract the cluster with the most points
main_cluster = filteredInVerticesPoints(labels == max_cluster_label,:);
main_cluster_sources = sourceLabels(labels == max_cluster_label);
% Extract all other points (other clusters and outliers)
other_points = filteredInVerticesPoints(labels ~= max_cluster_label & labels ~= -1,:);
other_points_sources = sourceLabels(labels ~= max_cluster_label);
% Extract outliers
outliers = filteredInVerticesPoints(labels == -1, :);
outliers_sources = sourceLabels(labels == -1);

% Vertices / Points
% Separate main_cluster into vertices and points (keep vertices/points)
shrink.inside.mainClusterVertices = main_cluster(main_cluster_sources == 1, :);
shrink.inside.mainClusterPoints = main_cluster(main_cluster_sources == 2, :);
% Idx of reference vertices/points inside
numInsideVertices = size(shrink.inside.filteredInVertices,1);
indicesVertices = 1:numInsideVertices;
indicesPoints = numInsideVertices+1 : numInsideVertices+size(shrink.inside.filteredInPoints,1);
isMaxCluster = false(length(labels), 1);
isMaxCluster(labels == max_cluster_label) = true;
% Idx of reference vertices inside
shrink.inside.mainClusterVerticesMask = false(size(shrink.inside.filteredInVerticesMask,1),1);
shrink.inside.mainClusterVerticesIdx = shrink.inside.filteredInVerticesIdx(isMaxCluster(indicesVertices),:);
shrink.inside.mainClusterVerticesMask(shrink.inside.mainClusterVerticesIdx) = true;
% Find faces with vertices (main_cluster)
isVertexInside = false(max(refFaces(:)), 1);
isVertexInside(shrink.inside.mainClusterVerticesIdx) = true;
allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
shrink.inside.mainClusterFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)
% Idx of point base (gridPoints)
shrink.inside.mainClusterPointsMask = false(size(shrink.inside.filteredInPointsMask,1), 1); % Logical mask of point base (gridPoints)
shrink.inside.mainClusterPointsIdx = shrink.inside.filteredInPointsIdx(isMaxCluster(indicesPoints),:);
shrink.inside.mainClusterPointsMask(shrink.inside.mainClusterPointsIdx) = true;

% Extract and separate other clusters into cells
otherClustersVertices = cell(length(unique_labels) - 1, 1); % Minus 1 because of the main cluster
otherClustersPoints = cell(length(unique_labels) - 1, 1); % Minus 1 because of the main cluster
current_index = 1;
for i = 1:length(unique_labels)
    if unique_labels(i) ~= max_cluster_label
        clusterVerticesPoints = filteredInVerticesPoints(labels == unique_labels(i), :);
        clusterSources = sourceLabels(labels == unique_labels(i));
        otherClustersVertices{current_index} = clusterVerticesPoints(clusterSources == 1, :); % Vertices
        otherClustersPoints{current_index} = clusterVerticesPoints(clusterSources == 2, :); % Points
        current_index = current_index + 1;
    end
end
shrink.inside.otherClustersVertices = otherClustersVertices;
shrink.inside.otherClustersPoints = otherClustersPoints;
% Separate outliers into vertices and points
shrink.inside.outliersVertices = outliers(outliers_sources == 1, :);
shrink.inside.outliersPoints = outliers(outliers_sources == 2, :);

    function confirmCluster(clusterLabel, keepCluster)
        if keepCluster
            % Merge the selected cluster with the main cluster by updating labels
            labels(labels == clusterLabel) = max_cluster_label;
        end
        close(gcf); % Close the figure window
    end

disp(['identify boundary points (cluster): pelvis defect ',num2str(pelvisNum)])

end