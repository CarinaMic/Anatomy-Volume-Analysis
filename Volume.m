classdef Volume
    % Volume class generates volume model of a specific surface model (stl)
    
    % Developed by C.Micheler,
    % Department of Orthopaedics and Sportorthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
    % Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich
    
    
    properties 
        edge = struct();    % Identify edges
        patches = struct(); % Repair mesh: fill holes
        shrink = struct();  % Shrink wrap functions
        gridPoints = struct();  % Geometry filled with points 
        acRegion = struct(); % Acetabulum regions for defect loss
    end
    
    methods
        %% Constructer: generate object
        function obj = Volume()
            disp('class volume initialized')
        end

        %% Fill mesh with points (for intersection)
        function obj = fillPoints(obj,pelvisNum,refData,minDist,box)

            % Bounding Box (not aligned to axis)
            % Assuming: cornerPoints is a 8x3 matrix where each row is a corner point
            % Bounding box's axis directions
            dir1 = box.edgeVector(1,:)/norm(box.edgeVector(1,:)); % x % One box edge; right-hand-rule
            dir2 = box.edgeVector(2,:)/norm(box.edgeVector(2,:)); % y
            dir3 = box.edgeVector(3,:)/norm(box.edgeVector(3,:)); % z

            % Create the rotation matrix local -> global (box)
            rotMatrix = [dir1; dir2; dir3];

            % Bounding box's dimensions
            boxSize = [norm(box.edgeVector(1,:)), norm(box.edgeVector(2,:)), norm(box.edgeVector(3,:))];
            % Determine how many points fit in each dimension based on minDist
            % Desired minimum distance between points: minDist
            numCells = ceil(boxSize ./ minDist);
            % Generate points in bounding box's local coordinate system
            numPoints = prod(numCells); % Total number of points
            pointsLocal = zeros(numPoints, 3); % Preallocate the matrix for efficiency
            for idx = 1:numPoints
                [i, j, k] = ind2sub(numCells, idx); % Convert linear index to 3D index
                x = (i-1)*minDist;
                y = (j-1)*minDist;
                z = (k-1)*minDist;
                pointsLocal(idx,:) = [x, y, z];
            end

            % Transform Points to global cosy
            origin = box.cornerpoints(1,:);  % assuming this is the origin of the bounding box in world coordinates
            pointsAll = pointsLocal * rotMatrix;
            pointsAll = pointsAll + origin;
            obj.gridPoints.pointsBox = pointsAll;
            obj.gridPoints.pointsBoxNum = numPoints;

            % Check
            % figure
            % plot3(pointsAll(:,1), pointsAll(:,2), ...
            %     pointsAll(:,3), '.','Color','b', 'MarkerSize', 5);
            % hold on
            % trisurf(box.tri,...
            %     box.cornerpoints(:,1),...
            %     box.cornerpoints(:,2),...
            %     box.cornerpoints(:,3),...
            %     'FaceColor','r','EdgeColor','r','FaceAlpha',0.25);

            % Vertices (V) and faces (F) define the mesh
            V = refData.vertices;
            F = refData.faces;
            % Check which points are inside the mesh 
            inside = inpolyhedron(F, V, pointsAll);
            pointsInside = pointsAll(inside, :);
            obj.gridPoints.inside = pointsInside;
            obj.gridPoints.insideMask = inside; % base: points in box (pointsBox)
            obj.gridPoints.insideMaskIdx = find(inside);
            
            disp(['pelvis filled with points: pelvis ', num2str(pelvisNum)]);

        end
        
        %% Identify edge / boundary vertices 
        function obj = edging(obj,pelvisNum,importData)
            
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
            obj.edge.verticesNrCombi = boundaryEdges;
            obj.edge.verticesNr = unique(boundaryEdges(:)); % Flatten the array and then find unique vertices
            obj.edge.vertices = vertices(obj.edge.verticesNr,:);
                                                            
            disp(['edge calculated: pelvis defect ',num2str(pelvisNum)])
                                   
        end
        
        %% Identify edge / boundary loops 
        function obj = edgeLoop(obj,pelvisNum,importData)
            
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
            obj.edge.edgeLoops = reorderedLoops';
            obj.edge.loopsVerticesNrCombi = reorderedVerticesLoop';
            obj.edge.comPoint = reorderedComPoint';
            
            disp(['edge loops calculated: pelvis defect ',num2str(pelvisNum)])
            
        end
   
        %% Close holes and triangulate (patch) + refine patches
        function obj = meshPatch(obj,pelvisNum,importData,varargin)
            
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
            obj.patches.(type).numComponents = numComponents;

            % Determine the number of loops to process based on the type
            if any(strcmp(type, {'acentre', 'acetabulum'}))
                numLoops = 1; % Only process the first edge loop (acetabulum loop)
                numComponents = 0; % reset
            else
                numLoops = size(obj.edge.comPoint, 1); 
            end
            data_sorted = cell(numLoops, 1);
            Tri = cell(numLoops, 1);

            % Mesh edge loops
            % Matlab File Exchange: discrete_contour_mesh_patch
            % https://de.mathworks.com/matlabcentral/fileexchange/78901-discrete-contour-mesh-patch-2d-3d?s_tid=srchtitle
            % Mostly efficient for quasi flat contours, cause doesn't yet prevent from self intersecting triangles.
            for i = 1:numLoops
                data_sorted{i,1} = vertices(obj.edge.comPoint{i,1}(:),:); % comPoint: sorted edge loop points
                if strcmp(type, 'acentre') || strcmp(type, 'acetabulum') % only one loop
                    data_sorted{i,1} = [data_sorted{i,1}; extra_data]; % extra_data = coordinates 
                end
                [~,Tri{i,1},~] = discrete_contour_mesh_patch(data_sorted{i,1},'sorted');
            end

            % Store data to obj
            obj.patches.(type).faces = Tri;
            obj.patches.(type).vertices = data_sorted;
            % Compute the face normals
            for i = 1:numLoops
                V1 = obj.patches.(type).vertices{i,1}(obj.patches.(type).faces{i,1}(:,1),:);
                V2 = obj.patches.(type).vertices{i,1}(obj.patches.(type).faces{i,1}(:,2),:);
                V3 = obj.patches.(type).vertices{i,1}(obj.patches.(type).faces{i,1}(:,3),:);
                normals = cross(V2 - V1, V3 - V1, 2);
                % Normalize the normals
                obj.patches.(type).normals{i,1} = normals ./ vecnorm(normals, 2, 2);
                % Calculate the centroid of the triangular faces in the mesh
                faceVertices = cat(3, V1, V2, V3);
                obj.patches.(type).centreFaces{i,1} = mean(faceVertices, 3);
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
                    obj.patches.(type).refinedVertices{i,1} = patch.vertices;
                    obj.patches.(type).refinedFaces{i,1} = patch.faces;
                    % Compute the face normals
                    V1 = obj.patches.(type).refinedVertices{i,1}(obj.patches.(type).refinedFaces{i,1}(:,1),:);
                    V2 = obj.patches.(type).refinedVertices{i,1}(obj.patches.(type).refinedFaces{i,1}(:,2),:);
                    V3 = obj.patches.(type).refinedVertices{i,1}(obj.patches.(type).refinedFaces{i,1}(:,3),:);
                    normals = cross(V2 - V1, V3 - V1, 2);
                    % Normalize the normals
                    obj.patches.(type).refinedNormals{i,1} = normals ./ vecnorm(normals, 2, 2);
                    % Calculate the centroid of the triangular faces in the mesh
                    faceVertices = cat(3, V1, V2, V3);
                    obj.patches.(type).refinedCentreFaces{i,1} = mean(faceVertices, 3);
                end
            end

            disp(['meshed edge loops (patch - ', type, '): pelvis defect ', num2str(pelvisNum)]);           
            
        end
        
        %% Check for orientation of the triangulated patch (user-input)
        function obj = flipNormals(obj,pelvisNum,edgeNum,importData,patchVertices,patchFaces,patchNormals,patchCentreFaces,type)

            vertices = importData.vertices;
            faces = importData.faces;

            fig = figure;

            % Defect Mesh
            patch('Faces',faces,...
                'Vertices',vertices,...
                'FaceColor',[0.9 0.75 0.68], ...    % Face color
                'FaceAlpha',1,...                   % Transparency of the faces
                'EdgeColor',[0.502 0.502 0.502],... % Edge color
                'EdgeAlpha',0.25);                  % Transparency of the edges
            hold on
            % Meshed patch
            patch('Faces',patchFaces,...
                'Vertices',patchVertices,...
                'FaceColor',[0.9 0.75 0.68], ...    % Face color
                'FaceAlpha',1,...                   % Transparency of the faces
                'EdgeColor',[0.502 0.502 0.502],... % Edge color
                'EdgeAlpha',0.25);                  % Transparency of the edges
            hold on

            % Normals of the faces
            quiver3(patchCentreFaces(:,1),...
                patchCentreFaces(:,2),...
                patchCentreFaces(:,3),...
                patchNormals(:,1),...
                patchNormals(:,2),...
                patchNormals(:,3),...
                10,'Color',[0.89 0.447 0.133])
            hold on
            % Patch points
            plot3(patchVertices(:,1),patchVertices(:,2), patchVertices(:,3),...
                '.', 'Color', [0 0.396 0.741], 'MarkerSize', 1)

            grid minor
            title('pelvis ',num2str(pelvisNum));
            daspect([1 1 1]); % Axis ratio
            view(3);          % Incl. rotate3d with the mouse
            rotate3d on;      % Enable 3D rotation in the figure

            % Buttons:
            % https://de.mathworks.com/help/matlab/ref/uicontrol.html
            % https://de.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html
            % https://de.mathworks.com/help/matlab/ref/matlab.ui.control.uicontrol-properties.html
            flipButton = uicontrol('Parent',fig, 'Style','togglebutton', 'String','FLIP', ...
                'Units','Pixels', 'Position',[18 70 100 22], 'Callback', @flipCallback);
            confirmButton = uicontrol('Parent',fig,'Style','togglebutton','String','CONFIRM',...
                'FontWeight','bold','Units','Pixels','Position',[18 20 100 22]);

            obj.patches.(type).flip{edgeNum,1} = 0;
            
            while confirmButton.Value == 0
                %drawnow;
                pause(0.05);
            end

            % Plot is closed
            close(fig)

            % Functions of push buttons
            function flipCallback(src,~)
                obj.patches.(type).flip{edgeNum,1} = src.Value;
            end

            disp(['checked orientation of the triangulation (patch): pelvis defect ',num2str(pelvisNum)])

        end
        
        %% Combine mesh with patches
        function obj = remeshPatch(obj,pelvisNum,importData,slavesFaces,slavesVertices,type)

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
            for i = (obj.patches.(type).numComponents+1):size(slavesFaces,1)
                % Combine faces
                mapSlaveFaces = slavesFaces{i};
                if obj.patches.(type).flip{i,1} == 1
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
            obj.patches.(type).comVertices = currentVertices;
            obj.patches.(type).comFaces = currentFaces;
            % Compute the face normals
            V1 = obj.patches.(type).comVertices(obj.patches.(type).comFaces(:, 1), :);
            V2 = obj.patches.(type).comVertices(obj.patches.(type).comFaces(:, 2), :);
            V3 = obj.patches.(type).comVertices(obj.patches.(type).comFaces(:, 3), :);
            normals = cross(V2 - V1, V3 - V1, 2);
            % Normalize the normals
            obj.patches.(type).comNormals = normals ./ vecnorm(normals, 2, 2);
            % Calculate the centroid of the triangular faces in the mesh
            faceVertices = cat(3, V1, V2, V3);
            obj.patches.(type).comCentreFaces = mean(faceVertices, 3);
     
            disp(['combined mesh with patches: pelvis defect ',num2str(pelvisNum)])

        end
        
        %% Shrink Wrap with boundary function over vertices
        function obj = shrinkBound(obj,pelvisNum,acentre,vertices)           

            % Add the additional point acentre
            combinedVertices = vertices;
            combinedVertices = [combinedVertices; acentre];

            % ShrinkWrap with boundary function
            % https://de.mathworks.com/help/matlab/ref/boundary.html#buh3c7k-1-s
            % Shrink factor is between 0 and 1. 0 gives the convex hull, 
            % and 1 gives a compact boundary that envelops the points. (default 0.5)
            obj.shrink.bound.faces = boundary(combinedVertices,0.8); % define shrink factor
            obj.shrink.bound.vertices = combinedVertices;

            disp(['shrinked mesh (boundary): pelvis defect ',num2str(pelvisNum)])
            
        end
        
        %% Shrink Wrap with alphaShape (with user-input)
        function obj = shrinkAlphaUI(obj,pelvisNum,type,aRadius,numData,importData,addPoints)

            if isfield(importData, 'comVertices') && isfield(importData, 'comFaces')
                verticesImport = importData.comVertices;
                facesImport = importData.comFaces;
            else
                verticesImport = importData.vertices;
                facesImport = importData.faces;
            end
            if nargin < 7
                acentre = importData.acentre;
                allVerticesPoints = [verticesImport; acentre];
            else
                allVerticesPoints = [verticesImport; addPoints];
            end

            % First call
            %aRadius = 1; % Start value
            shrinkWrap = alphaShape(allVerticesPoints,aRadius);
            shrinkFaces = boundaryFacets(shrinkWrap);
            
            fig = figure;
            beep;
            % Define view directions
            views = [0, 90; 270, 0; 0, 0; 120, 10]; % [azimuth, elevation]

            % Iterate through the four desired subplots
            for i = 1:4
                % Create subplot
                axHandles(i) = subplot(2, 2, i);  % 2x2 grid of subplots

                % Display defect surface
                p(i) = patch('Faces',shrinkFaces,...
                    'Vertices',shrinkWrap.Points,...
                    'FaceColor',[0.8 0.8 0.8], ...      % Face color
                    'FaceAlpha',1,...                   % Transparency of the faces
                    'EdgeColor','none',... % Edge color [0.502 0.502 0.502]
                    'EdgeAlpha',0.25,...                % Transparency of the edges
                    ... % Ligthing for 3d effect
                    'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
                    'AmbientStrength', 0.5);
                hold on

                % Display vertices
                plot3(verticesImport(:,1), verticesImport(:,2), verticesImport(:,3),...
                    '.','Color',[0.2 0.2 0.2],'MarkerSize', 1)
                hold on
                % Display edge (acetabulum edge)
                for k = 1:obj.patches.all.numComponents
                    plot3(obj.edge.verticesLoops{k}(:,1),obj.edge.verticesLoops{k}(:,2),...
                        obj.edge.verticesLoops{k}(:,3),...
                        '*','Color',[0.89 0.447 0.133],'MarkerSize',3)
                end
                % Display acentre/acetabulum
                if isfield(importData, 'acentre')
                    acentre = importData.acentre;
                    plot3(acentre(:,1), acentre(:,2), acentre(:,3),...
                        '*', 'Color', [0 0.396 0.741], 'MarkerSize', 10);
                end

                % Display defect
                if (~isempty(verticesImport) && ~isempty(facesImport)) == 1
                    % Display defect
                    hold on
                    patch('Faces',facesImport,...
                        'Vertices',verticesImport,...
                        'FaceColor',[0.9 0.75 0.68], ...    % Face color
                        'FaceAlpha',1,...                   % Transparency of the faces
                        'EdgeColor','none',... % Edge color [0.502 0.502 0.502]
                        'EdgeAlpha',0.25,...                % Transparency of the edges
                        ... % Ligthing for 3d effect
                        'FaceLighting', 'gouraud', ...      % Choose a lighting algorithm
                        'AmbientStrength', 0.5);
                    light('Position', [1 1 5], 'Style', 'infinite');
                end

                % Display addPoints
                if nargin > 6 
                   plot3(addPoints(:,1), addPoints(:,2), addPoints(:,3),...
                    '.', 'Color', [0 101/255 189/255], 'MarkerSize', 2); 
                end

                % Axis
                xl = xlim;          % Return x-axis limits
                yl = ylim;          % Return y-axis limits
                zl = zlim;          % Return z-axis limits
                axisAdd = 20;
                axis([xl(1,1)-axisAdd xl(1,2)+axisAdd yl(1,1)-axisAdd yl(1,2)+axisAdd ...
                    zl(1,1)-axisAdd zl(1,2)+axisAdd]);

                % Format and display properties
                grid minor
                title('Shrink Wrap: alpha Radius');
                daspect([1, 1, 1]); % Equal aspect ratio for the axes
                %view(3);
                rotate3d on;        % Enable 3D rotation in the figure
                % Different view for each subplot
                view(views(i,1), views(i,2));
            end

            % Buttons:
            % https://de.mathworks.com/help/matlab/ref/uicontrol.html
            % https://de.mathworks.com/help/matlab/creating_guis/write-callbacks-using-the-programmatic-workflow.html
            % https://de.mathworks.com/help/matlab/ref/matlab.ui.control.uicontrol-properties.html        
            plusButton = uicontrol('Parent',fig,'Style','push','String','Plus alphaRadius','Units','Pixels',...
                'Position',[18 120 100 25],'Callback',@plotPlusView);
            minusButton = uicontrol('Parent',fig,'Style','push','String','Minus alpaRadius','Units','Pixels',...
                'Position',[18 70 100 25],'Callback',@plotMinusView);            
            confirmButton = uicontrol('Parent',fig,'Style','togglebutton','String','CONFIRM',...
                'FontWeight','bold','Units','Pixels','Position',[18 20 100 25]);    
            showButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Show alphaRadius', ...
                'Units', 'Pixels', 'Position', [18 170 100 25], ...
                'Callback', @refreshPlots);
            % Annotation to display the current alphaRadius
            alphaText = annotation('textbox', [0.25, 0, 0.6, 0.05], 'String', ['Set Alpha Radius = ', num2str(aRadius)], ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 10);
                       
            while confirmButton.Value == 0
                % shrinkWrap = alphaShape(vertices,aRadius);
                % shrinkFaces = boundaryFacets(shrinkWrap);
                drawnow  
                pause(0.2);
            end
            
            % Final shrink wrap after user confirmation
            shrinkWrap = alphaShape(allVerticesPoints,aRadius);
            shrinkFaces = boundaryFacets(shrinkWrap);

            % Store results in the object 
            obj.shrink.(type).radius = aRadius;
            obj.shrink.(type).alpha = alphaShape(allVerticesPoints,aRadius);
            % ShrinkWrap: vertices -> not all vertices are used
            shrinkFaces = boundaryFacets(obj.shrink.(type).alpha);
            obj.shrink.(type).faces = shrinkFaces;
            obj.shrink.(type).allVerticesPoints = allVerticesPoints; % all vertices + acentre
            % Logical index for all used vertices/vertices in the mesh
            % original data (vertices): without patches, only original mesh
            uniVertices = unique(shrinkFaces(:)); 
            usedDefectVerticesImport = uniVertices(uniVertices <= numData);
            obj.shrink.(type).usedDefectVerticesMask = false(numData,1); % numData: size of original vertices data
            obj.shrink.(type).usedDefectVerticesMask(usedDefectVerticesImport) = true;
            obj.shrink.(type).usedDefectVerticesCount = length(usedDefectVerticesImport) / numData;
            obj.shrink.(type).usedDefectVertices = verticesImport(usedDefectVerticesImport, :);
            % All vertices/points
            obj.shrink.(type).usedVerticesPointsMask = false(size(allVerticesPoints, 1), 1);
            obj.shrink.(type).usedVerticesPointsMask(uniVertices) = true;
            obj.shrink.(type).usedVerticesPoints = allVerticesPoints(obj.shrink.(type).usedVerticesPointsMask, :);
            % New faces structure: renumbering faces (faces structure with new numbering of usedVerticesPoints)
            vertexIndexMapping = zeros(size(obj.shrink.(type).usedVerticesPointsMask));
            vertexIndexMapping(obj.shrink.(type).usedVerticesPointsMask) = 1:nnz(obj.shrink.(type).usedVerticesPointsMask);
            remappedFaces = vertexIndexMapping(obj.shrink.(type).faces);
            obj.shrink.(type).usedFaces = remappedFaces;
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
            obj.shrink.(type).volume = abs(volumeAlpha);
             
            % Plot is closed
            close(fig)
            
            disp(['Shrink-wrapped mesh (alphaShape): pelvis defect ',num2str(pelvisNum)])
            
            function plotPlusView(plusButton,ax)
                    aRadius = aRadius + 1;
                    updateAlphaText();
            end
            function plotMinusView(minusButton,ax)
                if aRadius > 0
                    aRadius = aRadius - 1;
                    updateAlphaText();
                end
            end
            function updateAlphaText()
                % Update the String property of the alphaText handle
                set(alphaText, 'String', ['Set Alpha Radius = ', num2str(aRadius)]);
            end
            function refreshPlots(src, event)
                % recompute alphaShape
                shrinkWrap = alphaShape(allVerticesPoints,aRadius);
                shrinkFaces = boundaryFacets(shrinkWrap);
                % Refresh
                for i = 1:4
                    set(p(i), 'Faces', shrinkFaces, 'Vertices', shrinkWrap.Points);
                    title(axHandles(i), ['Shrink Wrap: alpha Radius ',num2str(aRadius)]);
                end
                
                % Check if mesh is closed and display the result
                if aRadius > 0 
                    % Create an edge list from the triangle faces.
                    edgesList = [shrinkFaces(:,1) shrinkFaces(:,2);
                        shrinkFaces(:,2) shrinkFaces(:,3);
                        shrinkFaces(:,3) shrinkFaces(:,1)];
                    % Sort each edge to ensure consistent ordering
                    edgesList = sort(edgesList, 2);
                    % Identify unique edges and count occurrences
                    [~, ~, edgeOccurrences] = unique(edgesList, 'rows', 'stable');
                    edgeCounts = accumarray(edgeOccurrences, 1);
                    % Check if any edge is not shared by exactly two triangles
                    isClosed = all(edgeCounts == 2);
                    if isClosed
                        meshStatus = 'Mesh: closed and manifold'; % closed and non-manifold mesh
                    else
                        meshStatus = 'Mesh: open or non-manifold';
                    end
                    % Display mesh status above buttons
                    uicontrol('Style','text', 'Position',[120 165 100 25], ...
                        'String',meshStatus, 'BackgroundColor', get(fig, 'Color'));
                end
                beep;
                % Refresh the figure
                drawnow
            end
                              
        end

        %% Shrink Wrap with alphaShape (without user-input)
        function obj = shrinkAlpha(obj, pelvisNum, type, aRadius, numData, importData, addPoints)

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
            obj.shrink.(type).radius = aRadius; 
            obj.shrink.(type).alpha = shrinkWrap; 
            obj.shrink.(type).faces = shrinkFaces; 
            obj.shrink.(type).allVerticesPoints = allVerticesPoints; % all vertices 

            % Logical index for all used vertices/vertices in the mesh 
            % original data (vertices): without patches, only original mesh
            uniVertices = unique(shrinkFaces(:)); 
            usedDefectVerticesImport = uniVertices(uniVertices <= numData);
            obj.shrink.(type).usedDefectVerticesMask = false(numData,1); % numData: size of original vertices data
            obj.shrink.(type).usedDefectVerticesMask(usedDefectVerticesImport) = true;
            obj.shrink.(type).usedDefectVerticesCount = length(usedDefectVerticesImport) / numData;
            obj.shrink.(type).usedDefectVertices = verticesImport(usedDefectVerticesImport, :);
            % All vertices/points
            obj.shrink.(type).usedVerticesPointsMask = false(size(allVerticesPoints, 1), 1);
            obj.shrink.(type).usedVerticesPointsMask(uniVertices) = true;
            obj.shrink.(type).usedVerticesPoints = allVerticesPoints(obj.shrink.(type).usedVerticesPointsMask, :);
            % New faces structure: renumbering faces (faces structure with new numbering of usedVerticesPoints)
            vertexIndexMapping = zeros(size(obj.shrink.(type).usedVerticesPointsMask));
            vertexIndexMapping(obj.shrink.(type).usedVerticesPointsMask) = 1:nnz(obj.shrink.(type).usedVerticesPointsMask);
            remappedFaces = vertexIndexMapping(obj.shrink.(type).faces);
            obj.shrink.(type).usedFaces = remappedFaces;

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
            obj.shrink.(type).volume = abs(volumeAlpha);

            disp(['shrink-wrapped mesh (alphaShape) calculated: pelvis defect ', num2str(pelvisNum)])

        end

         %% Shrink Wrap with alphaShape, more iterations for isclosed (without user-input) 
         function obj = shrinkAlphaLoop(obj, pelvisNum, type, numLoops, aRadius, numData, importData, addPoints)

            if isfield(importData, 'comVertices') && isfield(importData, 'comFaces')
                verticesImport = importData.comVertices;
            else
                verticesImport = importData.vertices;
            end
            if nargin < 8
                acentre = importData.acentre;
                allVerticesPoints = [verticesImport; acentre];
            else
                allVerticesPoints = [verticesImport; addPoints];
            end

            % Variables to store results from the three attempts
            obj.shrink.(type).successfulRadii = zeros(1,numLoops);
            obj.shrink.(type).alphaShapes = cell(1,numLoops);
            obj.shrink.(type).shrinkFacesList = cell(1,numLoops);

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
                            obj.shrink.(type).successfulRadii(attempt) = currentRadius; 
                            obj.shrink.(type).alphaShapes{attempt} = shrinkWrap; 
                            obj.shrink.(type).shrinkFacesList{attempt} = shrinkFaces;
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
            obj.shrink.(type).radius = obj.shrink.(type).successfulRadii(1); 
            obj.shrink.(type).alpha = obj.shrink.(type).alphaShapes{1}; 
            obj.shrink.(type).faces = obj.shrink.(type).shrinkFacesList{1}; 
            obj.shrink.(type).allVerticesPoints = allVerticesPoints; % all vertices

            % Logical index for all used vertices/vertices in the mesh 
            % original data (vertices): without patches, only original mesh
            uniVertices = unique(shrinkFaces(:)); 
            usedDefectVerticesImport = uniVertices(uniVertices <= numData);
            obj.shrink.(type).usedDefectVerticesMask = false(numData,1); % numData: size of original vertices data
            obj.shrink.(type).usedDefectVerticesMask(usedDefectVerticesImport) = true;
            obj.shrink.(type).usedDefectVerticesCount = length(usedDefectVerticesImport) / numData;
            obj.shrink.(type).usedDefectVertices = verticesImport(usedDefectVerticesImport, :);
            % All vertices/points
            obj.shrink.(type).usedVerticesPointsMask = false(size(allVerticesPoints, 1), 1);
            obj.shrink.(type).usedVerticesPointsMask(uniVertices) = true;
            obj.shrink.(type).usedVerticesPoints = allVerticesPoints(obj.shrink.(type).usedVerticesPointsMask, :);
            % New faces structure: renumbering faces (faces structure with new numbering of usedVerticesPoints)
            vertexIndexMapping = zeros(size(obj.shrink.(type).usedVerticesPointsMask));
            vertexIndexMapping(obj.shrink.(type).usedVerticesPointsMask) = 1:nnz(obj.shrink.(type).usedVerticesPointsMask);
            remappedFaces = vertexIndexMapping(obj.shrink.(type).faces);
            obj.shrink.(type).usedFaces = remappedFaces;

            % Volume (divergence theorem)
            volumeAlpha = 0; % Initialize volume
            % Loop through each face
            for i = 1:size(obj.shrink.(type).shrinkFacesList{1}, 1)
                % Get the vertices of the face
                v1 = allVerticesPoints(obj.shrink.(type).shrinkFacesList{1}(i, 1), :); 
                v2 = allVerticesPoints(obj.shrink.(type).shrinkFacesList{1}(i, 2), :); 
                v3 = allVerticesPoints(obj.shrink.(type).shrinkFacesList{1}(i, 3), :); 
                % Calculate the signed volume of the tetrahedron
                tetraVolume = dot(v1, cross(v2, v3)) / 6;
                % Add to total volume
                volumeAlpha = volumeAlpha + tetraVolume;
            end
            % The volume should be positive
            obj.shrink.(type).volume = abs(volumeAlpha);

            disp(['shrink-wrapped mesh (alphaShape) iterations calculated: pelvis defect ', num2str(pelvisNum)])

        end
        
        %% alphaShape and find points of reference pelvis inside defect 
        function obj = shrinkInside(obj,pelvisNum,shrinkType,refData,addPoints,refGeometryIdx, ...
                pointsBoxNum,addPointsMaskIdx) % point base (gridPoints)

            refVertices = refData.vertices;
            refFaces = refData.faces;
                
            % ShrinkWrap (alphaShape) for rough hull
            % Use of previous function
            %obj = obj.shrinkAlpha(...);

            % Find vertices of the reference pelvis inside hull (alphaShape)
            FV.faces = obj.shrink.(shrinkType).faces; 
            FV.vertices = obj.shrink.(shrinkType).allVerticesPoints; 
            refVerticesInsideLog = inpolyhedron(FV,refVertices);
            % Correction: 
            % Part of the reference geometry part of the alpha shape
            % Not all points on the surface recognised as interior points with the inpolyhedron function
            refVerticesInsideLog(refGeometryIdx) = true;
            refVerticesInside = refVertices(refVerticesInsideLog,:); % Pelvis data (reference)

            % Find faces of the vertices
            allRefFacesInside = all(refVerticesInsideLog(refFaces), 2); % Faces with all vertex indices inside
            refFacesInside = refFaces(allRefFacesInside, :); % Pelvis data (reference)
            
            % Store results in the object
            obj.shrink.inside.refVerticesInside = refVerticesInside;
            obj.shrink.inside.refVerticesInsideMask = refVerticesInsideLog;
            obj.shrink.inside.refVerticesInsideIdx = find(refVerticesInsideLog);
            obj.shrink.inside.refFacesInside = refFacesInside;

            % Check which of the points inside reference pelvis are inside the alphaShape
            addInsideLog = inpolyhedron(FV, addPoints);
            addInside = addPoints(addInsideLog, :);
            obj.shrink.inside.refPointsInside = addInside;
            obj.shrink.inside.refPointsInsideMask = false(pointsBoxNum, 1); 
            obj.shrink.inside.refPointsInsideIdx = addPointsMaskIdx(addInsideLog,:); 
            obj.shrink.inside.refPointsInsideMask(obj.shrink.inside.refPointsInsideIdx) = true;
            
            disp(['shrinked mesh and found points inside: pelvis defect ',num2str(pelvisNum)])

        end

        %% Identify inner and outer points (reference vertices and points inside) - with vertex-to-nearest-neighbour to centreFaces of the Normals
        function obj = inoutPoints(obj,pelvisNum,importData,refData)

            refFaces = refData.faces;
            % Combined mesh: mesh of pelvis defect and refined patches (without acetabulum patch)
            centreFaces = importData.comCentreFaces;
            normals = importData.comNormals;
            verticesMesh = importData.comVertices; %%%
            facesMesh = importData.comFaces; %%%

            % Vector: point to nearest centreFaces -> centreFaces Normal
            verticesRef = centreFaces;
            verticesInside = obj.shrink.inside.refVerticesInside;
            pointsInside = obj.shrink.inside.refPointsInside;
            % Create a KD-tree searcher for verticesB
            MdlB = KDTreeSearcher(verticesRef);
            % Find the nearest neighbor in verticesB for each vertex in verticesA
            K = 1;
            [idxVertices, distVertices] = knnsearch(MdlB, verticesInside, 'K', K);
            [idxPoints, distPoints] = knnsearch(MdlB, pointsInside, 'K', K);
            % Calculate the distances based on the nearest neighbor indices
            minDistancesVertices = sqrt(sum((verticesInside - verticesRef(idxVertices,:)).^2, 2));
            minDistancesPoints = sqrt(sum((pointsInside - verticesRef(idxPoints,:)).^2, 2));

            % Initialize 
            insideOutsideVertices = false(size(verticesInside, 1), 1);
            % Loop through each vertex 
            for i = 1:size(verticesInside, 1)
                % Get the nearest reference point
                nearestRefVertice = verticesRef(idxVertices(i),:);
                % Calculate the vector from the nearest normal origin to the point
                vectorNormalOriginToVertice = verticesInside(i,:) - nearestRefVertice;
                % Get the corresponding normal vector
                normalVectorVertices = normals(idxVertices(i),:);
                % Calculate the dot product
                dotProductVertice = dot(normalVectorVertices, vectorNormalOriginToVertice);
                % Determine if the point is inside or outside
                if dotProductVertice > 0
                    insideOutsideVertices(i) = false;  % Outside the pelvis defect
                else
                    insideOutsideVertices(i) = true; % Inside the pelvis defect
                end
            end
            obj.shrink.inside.inoutVertices = insideOutsideVertices; % Logical mask of reference vertices
            outsideVertices = verticesInside(~insideOutsideVertices,:);
            obj.shrink.inside.outVertices = outsideVertices; 
            insideVertices = verticesInside(insideOutsideVertices,:);
            obj.shrink.inside.inVertices = insideVertices;
            % Logical mask of reference pelvis
            obj.shrink.inside.inVerticesMask = false(size(obj.shrink.inside.refVerticesInsideMask,1),1); 
            obj.shrink.inside.inVerticesIdx = obj.shrink.inside.refVerticesInsideIdx(insideOutsideVertices,:); 
            obj.shrink.inside.inVerticesMask(obj.shrink.inside.inVerticesIdx) = true;
            % Find faces of the vertices
            isVertexInside = false(max(refFaces(:)), 1);
            isVertexInside(obj.shrink.inside.inVerticesIdx) = true;
            allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
            obj.shrink.inside.inFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)

            % Initialize 
            insideOutsidePoints = false(size(pointsInside, 1), 1);
            % Loop through each vertex 
            for i = 1:size(pointsInside, 1)
                % Get the nearest reference point
                nearestRefPoint = verticesRef(idxPoints(i), :);
                % Calculate the vector from the nearest normal origin to the point
                vectorNormalOriginToPoint = pointsInside(i, :) - nearestRefPoint;
                % Get the corresponding normal vector
                normalVectorPoints = normals(idxPoints(i), :);
                % Calculate the dot product
                dotProductPoint = dot(normalVectorPoints, vectorNormalOriginToPoint);
                % Determine if the point is inside or outside
                if dotProductPoint > 0
                    insideOutsidePoints(i) = false;  % Outside the pelvis defect
                else
                    insideOutsidePoints(i) = true; % Inside the pelvis defect
                end
            end
            obj.shrink.inside.inoutPoints = insideOutsidePoints; % Logical mask of reference points
            outsidePoints = pointsInside(~insideOutsidePoints,:);
            obj.shrink.inside.outPoints = outsidePoints; 
            insidePoints = pointsInside(insideOutsidePoints,:);
            obj.shrink.inside.inPoints = insidePoints;
            % Logical mask of point base (gridPoints)
            obj.shrink.inside.inPointsMask = false(size(obj.shrink.inside.refPointsInsideMask,1), 1); % Logical mask of point base (gridPoints)
            obj.shrink.inside.inPointsIdx = obj.shrink.inside.refPointsInsideIdx(insideOutsidePoints,:);
            obj.shrink.inside.inPointsMask(obj.shrink.inside.inPointsIdx) = true; % Idx of points base            

            disp(['identify ínner/outer points (dot product): pelvis defect ',num2str(pelvisNum)])
            
        end

        %% Identify areas with low point density (outliers) - also for clustering
        function obj = densityPoints(obj,pelvisNum,threshold,refData)

            refFaces = refData.faces;
            % Vertices / Points
            insideVerticesPoints = [obj.shrink.inside.inVertices; obj.shrink.inside.inPoints];

            % Density
            bandwidth = 1; % 0.5 / 0.6 / 1
            % Density estimation
            tic
            density = mvksdensity(insideVerticesPoints, insideVerticesPoints, 'Bandwidth', bandwidth);
            obj.shrink.inside.executionTimeDensity = toc;
            obj.shrink.inside.pointCloudDensity = density;

            % Normalisation
            minDensity = min(density);
            maxDensity = max(density);
            normDensity = (density - minDensity) / (maxDensity - minDensity);
            obj.shrink.inside.pointCloudNormDensity = normDensity;
            % Filter out points with low point density
            thresholdLowDensity = prctile(normDensity,threshold); 
            lowDensity = normDensity <= thresholdLowDensity;
            obj.shrink.inside.inoutDensity = ~lowDensity;
            obj.shrink.inside.filteredOutAll = insideVerticesPoints(lowDensity,:);
            obj.shrink.inside.filteredInAll = insideVerticesPoints(~lowDensity,:);
            % Vertices / Points
            numInsideVertices = size(obj.shrink.inside.inVertices,1);
            indicesVertices = 1:numInsideVertices;
            indicesPoints = numInsideVertices+1 : numInsideVertices+size(obj.shrink.inside.inPoints,1);
            obj.shrink.inside.filteredInVertices = insideVerticesPoints(~lowDensity(indicesVertices),:);
            obj.shrink.inside.filteredOutVertices = insideVerticesPoints(lowDensity(indicesVertices),:);      
            obj.shrink.inside.filteredInPoints = obj.shrink.inside.inPoints(~lowDensity(indicesPoints),:);
            obj.shrink.inside.filteredOutPoints = obj.shrink.inside.inPoints(lowDensity(indicesPoints),:); 
            % Idx of reference pelvis (vertices)
            obj.shrink.inside.filteredInVerticesMask = false(size(obj.shrink.inside.inVerticesMask,1),1); 
            obj.shrink.inside.filteredInVerticesIdx = obj.shrink.inside.inVerticesIdx(~lowDensity(indicesVertices),:);
            obj.shrink.inside.filteredInVerticesMask(obj.shrink.inside.filteredInVerticesIdx) = true;
            % Find faces of the vertices
            isVertexInside = false(max(refFaces(:)), 1);
            isVertexInside(obj.shrink.inside.filteredInVerticesIdx) = true;
            allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
            obj.shrink.inside.filteredInFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)
            % Idx of point base (gridPoints)
            obj.shrink.inside.filteredInPointsMask = false(size(obj.shrink.inside.inPointsMask,1), 1); % Logical mask of point base (gridPoints)
            obj.shrink.inside.filteredInPointsIdx = obj.shrink.inside.inPointsIdx(~lowDensity(indicesPoints),:);
            obj.shrink.inside.filteredInPointsMask(obj.shrink.inside.filteredInPointsIdx) = true;         
            
            % Pelvis defect with viridis color
            normDensity(normDensity < 0) = 0;   % Clamp values below the range to 0
            normDensity(normDensity > 1) = 1;   % Clamp values above the range to 1
            % Apply colormap (recommended: viridis)
            % 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
            colourNum = 32768;
            colourMap = viridis(colourNum);
            colourIdx = round(normDensity * (colourNum-1)) + 1;
            rgbColour = colourMap(colourIdx,:);
            obj.shrink.inside.densityRGB = rgbColour;

            disp(['identify inner/outer points (density): pelvis defect ',num2str(pelvisNum)])

        end

        %% Identify clusters and outliers of the point cloud
        function obj = clusterPoints(obj,pelvisNum,epsilon,minPts,importData,refData)

            verticesMesh = importData.comVertices;
            facesMesh = importData.comFaces;
            refFaces = refData.faces;
            % Vertices / Points 
            filteredInVerticesPoints = [obj.shrink.inside.filteredInVertices; obj.shrink.inside.filteredInPoints];
            % Mark the sources of points: 1 for vertices, 2 for points
            sourceLabels = [ones(size(obj.shrink.inside.filteredInVertices, 1), 1); ...
                ones(size(obj.shrink.inside.filteredInPoints, 1), 1) * 2];

            % Cluster analysis with DBSCAN
            tic
            labels = dbscan(filteredInVerticesPoints, epsilon, minPts);
            obj.shrink.inside.executionTimeCluster = toc;
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
                        filteredInVertices = NaN(max(obj.shrink.inside.inVerticesIdx), 3);
                        filteredInVertices(obj.shrink.inside.inVerticesIdx, :) = obj.shrink.inside.inVertices;
                        patch('Faces',obj.shrink.inside.filteredInFaces,...
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
            obj.shrink.inside.clusterLabels = labels;
            obj.shrink.inside.clusterMainLabel = max_cluster_label;
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
            obj.shrink.inside.mainClusterVertices = main_cluster(main_cluster_sources == 1, :);
            obj.shrink.inside.mainClusterPoints = main_cluster(main_cluster_sources == 2, :);
            % Idx of reference vertices/points inside
            numInsideVertices = size(obj.shrink.inside.filteredInVertices,1);
            indicesVertices = 1:numInsideVertices;
            indicesPoints = numInsideVertices+1 : numInsideVertices+size(obj.shrink.inside.filteredInPoints,1);
            isMaxCluster = false(length(labels), 1);
            isMaxCluster(labels == max_cluster_label) = true;
            % Idx of reference vertices inside
            obj.shrink.inside.mainClusterVerticesMask = false(size(obj.shrink.inside.filteredInVerticesMask,1),1);
            obj.shrink.inside.mainClusterVerticesIdx = obj.shrink.inside.filteredInVerticesIdx(isMaxCluster(indicesVertices),:); 
            obj.shrink.inside.mainClusterVerticesMask(obj.shrink.inside.mainClusterVerticesIdx) = true;                    
            % Find faces with vertices (main_cluster) 
            isVertexInside = false(max(refFaces(:)), 1);
            isVertexInside(obj.shrink.inside.mainClusterVerticesIdx) = true;
            allRefFacesInside = all(isVertexInside(refFaces), 2); % Faces with all vertex indices inside
            obj.shrink.inside.mainClusterFaces = refFaces(allRefFacesInside, :); % Pelvis data (reference)
            % Idx of point base (gridPoints)
            obj.shrink.inside.mainClusterPointsMask = false(size(obj.shrink.inside.filteredInPointsMask,1), 1); % Logical mask of point base (gridPoints)
            obj.shrink.inside.mainClusterPointsIdx = obj.shrink.inside.filteredInPointsIdx(isMaxCluster(indicesPoints),:);
            obj.shrink.inside.mainClusterPointsMask(obj.shrink.inside.mainClusterPointsIdx) = true; 

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
            obj.shrink.inside.otherClustersVertices = otherClustersVertices;
            obj.shrink.inside.otherClustersPoints = otherClustersPoints;
            % Separate outliers into vertices and points
            obj.shrink.inside.outliersVertices = outliers(outliers_sources == 1, :);
            obj.shrink.inside.outliersPoints = outliers(outliers_sources == 2, :);

            function confirmCluster(clusterLabel, keepCluster)
                if keepCluster
                    % Merge the selected cluster with the main cluster by updating labels
                    labels(labels == clusterLabel) = max_cluster_label;
                end
                close(gcf); % Close the figure window
            end

            disp(['identify boundary points (cluster): pelvis defect ',num2str(pelvisNum)])

        end

        %% Thicken mesh (point cloud)
        function obj = thickenMesh(obj,pelvisNum, importData, thickness, step_size)

            centreFaces = importData.comCentreFaces;
            normals = importData.comNormals;
            num_steps = thickness / step_size;
            
            % Thicken mesh 
            verticesNew = [];
            for i = 1:num_steps
                % Move vertices inward along the normal direction
                verticesSteps = centreFaces - (i * step_size * normals); % Move points by step_size along normals (inside) 
                verticesNew = [verticesNew; verticesSteps];
            end

            % Filter out points outside mesh
            % Vector: point to nearest centreFaces -> centreFaces Normal
            verticesRef = centreFaces;
            pointsInside = verticesNew;
            % Create a KD-tree searcher for verticesB
            MdlB = KDTreeSearcher(verticesRef);
            % Find the nearest neighbor in verticesB for each vertex in verticesA
            K = 1;
            [idxPoints, distPoints] = knnsearch(MdlB, pointsInside, 'K', K);
            % Calculate the distances based on the nearest neighbor indices
            minDistancesPoints = sqrt(sum((pointsInside - verticesRef(idxPoints,:)).^2, 2));

            % Initialize
            insideOutsidePoints = false(size(pointsInside, 1), 1);
            % Loop through each vertex 
            for i = 1:size(pointsInside, 1)
                % Get the nearest reference point
                nearestRefPoint = verticesRef(idxPoints(i), :);
                % Calculate the vector from the nearest normal origin to the point
                vectorNormalOriginToPoint = pointsInside(i, :) - nearestRefPoint;
                % Get the corresponding normal vector
                normalVectorPoints = normals(idxPoints(i), :);
                % Calculate the dot product
                dotProductPoint = dot(normalVectorPoints, vectorNormalOriginToPoint);
                % Determine if the point is inside or outside
                if dotProductPoint > 0
                    insideOutsidePoints(i) = false;  % Outside the pelvis defect
                else
                    insideOutsidePoints(i) = true; % Inside the pelvis defect
                end
            end
            obj.shrink.inside.inoutVerticesThicken = insideOutsidePoints; % Logical mask
            insidePoints = pointsInside(insideOutsidePoints,:);
            obj.shrink.inside.verticesThicken = insidePoints;
           
            disp(['thicken mesh: pelvis defect ',num2str(pelvisNum)])
        end

        %% Surface reconstruction of point cloud with crust function
        function obj = crust(obj, pelvisNum, importData, addPoints)

            verticesImport = importData.comVertices;
            vertices = [verticesImport; addPoints];

            % Matlab File Exchange: Surface Reconstruction From Scattered Points Cloud 
            % https://de.mathworks.com/matlabcentral/fileexchange/63730-surface-reconstruction-from-scattered-points-cloud
            [faces,norm] = MyRobustCrust(vertices);
            obj.shrink.crust.faces = faces;
            obj.shrink.crust.normals = norm;
            obj.shrink.crust.vertices = vertices;

            figure(1)
            hold on
            title('Output Triangulation','fontsize',14)
            axis equal
            trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
            view(3);
            axis vis3d

            disp(['surface reconstruction: pelvis defect ',num2str(pelvisNum)])
        end

        %% Vertex-to-Nearest-Neighbour Distance
        function obj = shapeVertexNeighbour(obj, pelvisNum, verticesRef, vertices, faces)

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
            obj.shrink.comp.nearVertexFace = mean(nearVertexMatrix,2);
            % Normalisation
            lowerBound = min(obj.shrink.comp.nearVertexFace);
            upperBound = max(obj.shrink.comp.nearVertexFace);
            % Normalize data within this range
            normalData = (obj.shrink.comp.nearVertexFace - lowerBound) / (upperBound - lowerBound);
            normalData(normalData < 0) = 0;   % Clamp values below the range to 0
            normalData(normalData > 1) = 1;   % Clamp values above the range to 1
            % Apply colormap (recommended: viridis)
            % 16-bit colour: alpha: 1bit R: 5bit G: 5bit B: 5bit -> 32768 colors
            colourNum = 32768;
            colourMap = viridis(colourNum);
            colourIdx = round(normalData * (colourNum-1)) + 1;
            rgbColour = colourMap(colourIdx,:);

            obj.shrink.comp.nearVertexFaceNorm = normalData;
            obj.shrink.comp.nearVertexFaceColour = rgbColour;
            obj.shrink.comp.nearVertex = minDistances;
            obj.shrink.comp.nearVertexMax = max(minDistances);
            obj.shrink.comp.nearVertexMean = mean(minDistances);
            obj.shrink.comp.nearVertexStd = std(minDistances);

            disp(['Vertex-to-Nearest-Neighbour: pelvis ', num2str(pelvisNum)]);

        end

        %% Fill pelvis defect with points (Pre-filter with boundary box)
        function obj = fillDefect(obj,pelvisNum,inputPoints,box,meshFaces,meshVertices)

            % Bounding box (not aligned to axis)
            % Bounding box's axis directions
            dir1 = box.edgeVector(1,:)/norm(box.edgeVector(1,:)); % x % One box edge; right-hand-rule
            dir2 = box.edgeVector(2,:)/norm(box.edgeVector(2,:)); % y
            dir3 = box.edgeVector(3,:)/norm(box.edgeVector(3,:)); % z

            % Create the rotation matrix local -> global (box)
            rotMatrix = [dir1; dir2; dir3];

            % Rotate the input points and the box vertices to the aligned space
            rotatedPoints = inputPoints * rotMatrix';
            rotatedBox = box.cornerpoints * rotMatrix';
            % Check if the point is inside the aligned box 
            minX = min(rotatedBox(:, 1));
            maxX = max(rotatedBox(:, 1));
            minY = min(rotatedBox(:, 2));
            maxY = max(rotatedBox(:, 2));
            minZ = min(rotatedBox(:, 3));
            maxZ = max(rotatedBox(:, 3));
            % Check if each coordinate of each point is within the bounds
            isInsideX = rotatedPoints(:, 1) >= minX & rotatedPoints(:, 1) <= maxX;
            isInsideY = rotatedPoints(:, 2) >= minY & rotatedPoints(:, 2) <= maxY;
            isInsideZ = rotatedPoints(:, 3) >= minZ & rotatedPoints(:, 3) <= maxZ;
            % Combine the checks for all coordinates
            isInsideAll = isInsideX & isInsideY & isInsideZ; % logical
            obj.gridPoints.boxInsideMask = isInsideAll;
            % Extract the points that are inside the box
            rotatedPointsInside = rotatedPoints(isInsideAll,:); 
            obj.gridPoints.boxInsideIdx = find(isInsideAll); % idx

            % Transform points to global cosy
            pointsInside = rotatedPointsInside * rotMatrix;
            obj.gridPoints.boxInside = pointsInside;

            % Fill pelvis defects with points
            % Check which points are inside the mesh
            % Vertices/faces define the refined alpha shape mesh
            insideDefect = inpolyhedron(meshFaces, meshVertices, ... % Mesh
                obj.gridPoints.boxInside);
            obj.gridPoints.insideMask = false(size(obj.gridPoints.boxInsideMask,1), 1); % Logical mask of point base (gridPoints)
            obj.gridPoints.insideMaskIdx = obj.gridPoints.boxInsideIdx(insideDefect,:);
            obj.gridPoints.insideMask(obj.gridPoints.insideMaskIdx) = true;
            obj.gridPoints.inside = obj.gridPoints.boxInside(insideDefect,:);

            disp(['pelvis defect filled with points (pre-filter bounding box): pelvis defect ', num2str(pelvisNum)]);

        end
       
        %% Acetabulum region: cuboids for acetabulum region
        function obj = acCuboid(obj, limits) 

            % Large cuboid limits
            xmin = limits.xmin;
            xmax = limits.xmax;
            ymin = limits.ymin;
            ymax = limits.ymax;
            zmin = limits.zmin;
            zmax = limits.zmax;
            xmid = limits.xmid;
            ymid = limits.ymid;
            zmid = limits.zmid;

            % Large cuboid
            obj.acRegion.cuboid.all.limits = struct( ...
                'xmin', xmin, 'xmax', xmax, ...
                'ymin', ymin, 'ymax', ymax, ...
                'zmin', zmin, 'zmax', zmax);
            obj.acRegion.cuboid.all.vertices = [
                xmin ymin zmin;
                xmax ymin zmin;
                xmax ymax zmin;
                xmin ymax zmin;
                xmin ymin zmax;
                xmax ymin zmax;
                xmax ymax zmax;
                xmin ymax zmax];
            obj.acRegion.cuboid.all.faces = [
                1 2 3 4;
                5 6 7 8;
                1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                4 1 5 8];

            % q1: y+ / z+
            obj.acRegion.cuboid.c1.limits = struct( ...
                'xmin', xmin, 'xmax', xmax, ...
                'ymin', ymid, 'ymax', ymax, ...
                'zmin', zmid, 'zmax', zmax);
            obj.acRegion.cuboid.c1.vertices = [
                xmin ymid zmid;
                xmax ymid zmid;
                xmax ymax zmid;
                xmin ymax zmid;
                xmin ymid zmax;
                xmax ymid zmax;
                xmax ymax zmax;
                xmin ymax zmax];
            obj.acRegion.cuboid.c1.faces = [
                1 2 3 4;
                5 6 7 8;
                1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                4 1 5 8];

            % q2: y+ / z-
            obj.acRegion.cuboid.c2.limits = struct( ...
                'xmin', xmin, 'xmax', xmax, ...
                'ymin', ymid, 'ymax', ymax, ...
                'zmin', zmin, 'zmax', zmid);
            obj.acRegion.cuboid.c2.vertices = [
                xmin ymid zmin;
                xmax ymid zmin;
                xmax ymax zmin;
                xmin ymax zmin;
                xmin ymid zmid;
                xmax ymid zmid;
                xmax ymax zmid;
                xmin ymax zmid];
            obj.acRegion.cuboid.c2.faces = [
                1 2 3 4;
                5 6 7 8;
                1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                4 1 5 8];

            % q3: y- / z-
            obj.acRegion.cuboid.c3.limits = struct( ...
                'xmin', xmin, 'xmax', xmax, ...
                'ymin', ymin, 'ymax', ymid, ...
                'zmin', zmin, 'zmax', zmid);
            obj.acRegion.cuboid.c3.vertices = [
                xmin ymin zmin;
                xmax ymin zmin;
                xmax ymid zmin;
                xmin ymid zmin;
                xmin ymin zmid;
                xmax ymin zmid;
                xmax ymid zmid;
                xmin ymid zmid];
            obj.acRegion.cuboid.c3.faces = [
                1 2 3 4;
                5 6 7 8;
                1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                4 1 5 8];

            % q4: y- / z+
            obj.acRegion.cuboid.c4.limits = struct( ...
                'xmin', xmin, 'xmax', xmax, ...
                'ymin', ymin, 'ymax', ymid, ...
                'zmin', zmid, 'zmax', zmax);
            obj.acRegion.cuboid.c4.vertices = [
                xmin ymin zmid;
                xmax ymin zmid;
                xmax ymid zmid;
                xmin ymid zmid;
                xmin ymin zmax;
                xmax ymin zmax;
                xmax ymid zmax;
                xmin ymid zmax];
            obj.acRegion.cuboid.c4.faces = [
                1 2 3 4;
                5 6 7 8;
                1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                4 1 5 8];

            disp('acetabulum cuboid region created')

        end

        %% Acetabulum region: points/vertices inside sub-cuboids
        function obj = acCuboidInside(obj, insidePoints, insidePointsMask, vertices)

            % Cuboid limits
            c1 = obj.acRegion.cuboid.c1.limits;
            c2 = obj.acRegion.cuboid.c2.limits;
            c3 = obj.acRegion.cuboid.c3.limits;
            c4 = obj.acRegion.cuboid.c4.limits;

            % Grid points inside sub-cuboids (clear assignment of limit values)
            % c1: y+ / z+
            isInC1_points = insidePoints(:,1) >= c1.xmin & insidePoints(:,1) <= c1.xmax & ...
                insidePoints(:,2) >= c1.ymin & insidePoints(:,2) <= c1.ymax & ...
                insidePoints(:,3) >= c1.zmin & insidePoints(:,3) <= c1.zmax;
            % c2: y+ / z-
            isInC2_points = insidePoints(:,1) >= c2.xmin & insidePoints(:,1) <= c2.xmax & ...
                insidePoints(:,2) >= c2.ymin & insidePoints(:,2) <= c2.ymax & ...
                insidePoints(:,3) >= c2.zmin & insidePoints(:,3) <  c2.zmax;
            % c3: y- / z-
            isInC3_points = insidePoints(:,1) >= c3.xmin & insidePoints(:,1) <= c3.xmax & ...
                insidePoints(:,2) >= c3.ymin & insidePoints(:,2) <  c3.ymax & ...
                insidePoints(:,3) >= c3.zmin & insidePoints(:,3) <  c3.zmax;
            % c4: y- / z+
            isInC4_points = insidePoints(:,1) >= c4.xmin & insidePoints(:,1) <= c4.xmax & ...
                insidePoints(:,2) >= c4.ymin & insidePoints(:,2) <  c4.ymax & ...
                insidePoints(:,3) >= c4.zmin & insidePoints(:,3) <= c4.zmax;

            % Global indices relative to all grid points / pointsBox
            idxC1_grid = insidePointsMask(isInC1_points);
            idxC2_grid = insidePointsMask(isInC2_points);
            idxC3_grid = insidePointsMask(isInC3_points);
            idxC4_grid = insidePointsMask(isInC4_points);

            % Vertices inside sub-cuboids
            % c1: y+ / z+
            isInC1_vertices = vertices(:,1) >= c1.xmin & vertices(:,1) <= c1.xmax & ...
                vertices(:,2) >= c1.ymin & vertices(:,2) <= c1.ymax & ...
                vertices(:,3) >= c1.zmin & vertices(:,3) <= c1.zmax;
            % c2: y+ / z-
            isInC2_vertices = vertices(:,1) >= c2.xmin & vertices(:,1) <= c2.xmax & ...
                vertices(:,2) >= c2.ymin & vertices(:,2) <= c2.ymax & ...
                vertices(:,3) >= c2.zmin & vertices(:,3) <  c2.zmax;
            % c3: y- / z-
            isInC3_vertices = vertices(:,1) >= c3.xmin & vertices(:,1) <= c3.xmax & ...
                vertices(:,2) >= c3.ymin & vertices(:,2) <  c3.ymax & ...
                vertices(:,3) >= c3.zmin & vertices(:,3) <  c3.zmax;
            % c4: y- / z+
            isInC4_vertices = vertices(:,1) >= c4.xmin & vertices(:,1) <= c4.xmax & ...
                vertices(:,2) >= c4.ymin & vertices(:,2) <  c4.ymax & ...
                vertices(:,3) >= c4.zmin & vertices(:,3) <= c4.zmax;

            % Vertex indices relative to vertices array
            idxC1_vertex = find(isInC1_vertices);
            idxC2_vertex = find(isInC2_vertices);
            idxC3_vertex = find(isInC3_vertices);
            idxC4_vertex = find(isInC4_vertices);

            % Save results
            % Grid points / pointsBox
            obj.acRegion.cuboid.c1.gridIdx = idxC1_grid;
            obj.acRegion.cuboid.c2.gridIdx = idxC2_grid;
            obj.acRegion.cuboid.c3.gridIdx = idxC3_grid;
            obj.acRegion.cuboid.c4.gridIdx = idxC4_grid;
            % Vertices
            obj.acRegion.cuboid.c1.vertexIdx = idxC1_vertex;
            obj.acRegion.cuboid.c2.vertexIdx = idxC2_vertex;
            obj.acRegion.cuboid.c3.vertexIdx = idxC3_vertex;
            obj.acRegion.cuboid.c4.vertexIdx = idxC4_vertex;

            disp('grid points / vertices assigned to acetabulum sub-cuboids')

        end

        %% Acetabulum region: defect points inside sub-cuboids (grid points)
        function obj = acCuboidDefect(obj, isPointInsideDefect, refCuboid, voxelSize, pelvisNum)

            % Sub-cuboids of reference pelvis 
            c1.gridIdx = refCuboid.c1.gridIdx;
            c2.gridIdx = refCuboid.c2.gridIdx;
            c3.gridIdx = refCuboid.c3.gridIdx;
            c4.gridIdx = refCuboid.c4.gridIdx;      
            % Reference grid points per cuboid (total number)
            nRefGridC1 = numel(c1.gridIdx);
            nRefGridC2 = numel(c2.gridIdx);
            nRefGridC3 = numel(c3.gridIdx);
            nRefGridC4 = numel(c4.gridIdx);

            % Voxel size [mm]
            %voxelSize = [0.5 0.5 0.5];
            voxelVolume = prod(voxelSize);   % = 0.125 mm^3

            % global indices of defect grid points
            idxDefectGrid = find(isPointInsideDefect);
            % total number of defect grid points
            nDefectGrid_total = numel(idxDefectGrid);

            % defect grid points assigned to each cuboid 
            idxDefectGridC1 = intersect(idxDefectGrid, c1.gridIdx);
            idxDefectGridC2 = intersect(idxDefectGrid, c2.gridIdx);
            idxDefectGridC3 = intersect(idxDefectGrid, c3.gridIdx);
            idxDefectGridC4 = intersect(idxDefectGrid, c4.gridIdx);
            % defect grid points outside all cuboids
            idxDefectGridInCuboids = [idxDefectGridC1; idxDefectGridC2; idxDefectGridC3; idxDefectGridC4];
            idxDefectGridOutside = setdiff(idxDefectGrid, idxDefectGridInCuboids);
            % number of defect grid points per cuboid
            nDefectGridC1 = numel(idxDefectGridC1);
            nDefectGridC2 = numel(idxDefectGridC2);
            nDefectGridC3 = numel(idxDefectGridC3);
            nDefectGridC4 = numel(idxDefectGridC4);
            nDefectGridOutside = numel(idxDefectGridOutside);

            % Proportion of total defect grid points
            percentDefectGridC1 = 100 * nDefectGridC1 / nDefectGrid_total;
            percentDefectGridC2 = 100 * nDefectGridC2 / nDefectGrid_total;
            percentDefectGridC3 = 100 * nDefectGridC3 / nDefectGrid_total;
            percentDefectGridC4 = 100 * nDefectGridC4 / nDefectGrid_total;
            percentDefectGridOutside = 100 * nDefectGridOutside / nDefectGrid_total;
            % Loss percent relative to reference pelvis grid points inside each cuboid
            percentCuboidGridC1 = 100 * nDefectGridC1 / nRefGridC1;
            percentCuboidGridC2 = 100 * nDefectGridC2 / nRefGridC2;
            percentCuboidGridC3 = 100 * nDefectGridC3 / nRefGridC3;
            percentCuboidGridC4 = 100 * nDefectGridC4 / nRefGridC4;

            % Save (points)
            obj.acRegion.cuboid.totalDefect.defectGridIdx = idxDefectGrid; % total defect, not all sub-cuboids (defect can also be outside sub-cuboids)
            obj.acRegion.cuboid.totalDefect.nGridPoints = nDefectGrid_total; 
            obj.acRegion.cuboid.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume; 

            obj.acRegion.cuboid.c1.defectGridIdx = idxDefectGridC1; % grid points (idx) in sub-cuboid
            obj.acRegion.cuboid.c1.nGridPoints = nDefectGridC1; % number of grid points in sub-cuboid
            obj.acRegion.cuboid.c1.gridLossPercentDefect = percentDefectGridC1; % volume loss (percent) per defect
            obj.acRegion.cuboid.c1.gridLossPercentCuboid = percentCuboidGridC1; % volume loss (percent) per sub-cuboid
            obj.acRegion.cuboid.c1.gridVolumeLoss = nDefectGridC1 * voxelVolume; % volume loss (mm3)

            obj.acRegion.cuboid.c2.defectGridIdx = idxDefectGridC2;
            obj.acRegion.cuboid.c2.nGridPoints = nDefectGridC2;
            obj.acRegion.cuboid.c2.gridLossPercentDefect = percentDefectGridC2;
            obj.acRegion.cuboid.c2.gridLossPercentCuboid = percentCuboidGridC2;
            obj.acRegion.cuboid.c2.gridVolumeLoss = nDefectGridC2 * voxelVolume;

            obj.acRegion.cuboid.c3.defectGridIdx = idxDefectGridC3;
            obj.acRegion.cuboid.c3.nGridPoints = nDefectGridC3;
            obj.acRegion.cuboid.c3.gridLossPercentDefect = percentDefectGridC3;
            obj.acRegion.cuboid.c3.gridLossPercentCuboid = percentCuboidGridC3;
            obj.acRegion.cuboid.c3.gridVolumeLoss = nDefectGridC3 * voxelVolume;

            obj.acRegion.cuboid.c4.defectGridIdx = idxDefectGridC4;
            obj.acRegion.cuboid.c4.nGridPoints = nDefectGridC4;
            obj.acRegion.cuboid.c4.gridLossPercentDefect = percentDefectGridC4;
            obj.acRegion.cuboid.c4.gridLossPercentCuboid = percentCuboidGridC4;
            obj.acRegion.cuboid.c4.gridVolumeLoss = nDefectGridC4 * voxelVolume;

            obj.acRegion.cuboid.outside.defectGridIdx = idxDefectGridOutside;
            obj.acRegion.cuboid.outside.nGridPoints = nDefectGridOutside;
            obj.acRegion.cuboid.outside.gridLossPercentDefect = percentDefectGridOutside;
            obj.acRegion.cuboid.outside.gridVolumeLoss = nDefectGridOutside * voxelVolume;

            disp(['defect grid points/vertices assigned to acetabulum sub-cuboids: pelvis defect ', num2str(pelvisNum)]);

        end

        %% Acetabulum region: sectors for acetabulum region
        function obj = acSector(obj, limits)

            % Centre
            centre = limits.centre;
            x0 = centre(1);
            y0 = centre(2);
            z0 = centre(3);
            % x-depth
            xmin = limits.xmin;
            xmax = limits.xmax;
            % radii
            rInner   = limits.rInner;
            rOuter23 = limits.rOuter23;
            rOuter1  = limits.rOuter1;

            % Save common sector settings
            %obj.acRegion.sector.limits = limits;
            %obj.acRegion.sector.centre = centre;

            % Sector 4: inner full cylinder
            obj.acRegion.sector.s4.type = 'cylinder';
            obj.acRegion.sector.s4.centre = centre;
            obj.acRegion.sector.s4.xmin = xmin;
            obj.acRegion.sector.s4.xmax = xmax;
            obj.acRegion.sector.s4.rmin = 0;
            obj.acRegion.sector.s4.rmax = rInner;
            obj.acRegion.sector.s4.thetaMinDeg = limits.s4ThetaMinDeg;
            obj.acRegion.sector.s4.thetaMaxDeg = limits.s4ThetaMaxDeg;

            % Sector 2: annular sector on negative y side
            obj.acRegion.sector.s2.type = 'annularSector';
            obj.acRegion.sector.s2.centre = centre;
            obj.acRegion.sector.s2.xmin = xmin;
            obj.acRegion.sector.s2.xmax = xmax;
            obj.acRegion.sector.s2.rmin = rInner;
            obj.acRegion.sector.s2.rmax = rOuter23;
            obj.acRegion.sector.s2.thetaMinDeg = limits.s2ThetaMinDeg;
            obj.acRegion.sector.s2.thetaMaxDeg = limits.s2ThetaMaxDeg;

            % Sector 3: annular sector on positive y side
            obj.acRegion.sector.s3.type = 'annularSector';
            obj.acRegion.sector.s3.centre = centre;
            obj.acRegion.sector.s3.xmin = xmin;
            obj.acRegion.sector.s3.xmax = xmax;
            obj.acRegion.sector.s3.rmin = rInner;
            obj.acRegion.sector.s3.rmax = rOuter23;
            obj.acRegion.sector.s3.thetaMinDeg = limits.s3ThetaMinDeg;
            obj.acRegion.sector.s3.thetaMaxDeg = limits.s3ThetaMaxDeg;

            % Sector 1: superior annular sector
            obj.acRegion.sector.s1.type = 'annularSector';
            obj.acRegion.sector.s1.centre = centre;
            obj.acRegion.sector.s1.xmin = xmin;
            obj.acRegion.sector.s1.xmax = xmax;
            obj.acRegion.sector.s1.rmin = rInner;
            obj.acRegion.sector.s1.rmax = rOuter1;
            obj.acRegion.sector.s1.thetaMinDeg = limits.s1ThetaMinDeg;
            obj.acRegion.sector.s1.thetaMaxDeg = limits.s1ThetaMaxDeg;

            disp('acetabulum sector region created')

        end

        %% Acetabulum region: points / vertices inside sectors
        function obj = acSectorInside(obj, insidePoints, insidePointsIdx, vertices)

            % helper
            centre = obj.acRegion.sector.limits.centre;
            x0 = centre(1);
            y0 = centre(2);
            z0 = centre(3);
            % grid points
            x = insidePoints(:,1);
            y = insidePoints(:,2);
            z = insidePoints(:,3);

            rYZ = sqrt((y - y0).^2 + (z - z0).^2);
            thetaDeg = atan2d(y - y0, z - z0);   % 0° at +z

            % Sector 1
            s1 = obj.acRegion.sector.s1;
            isInS1 = x >= s1.xmin & x <= s1.xmax & ...
                rYZ >= s1.rmin & rYZ <= s1.rmax & ...
                thetaDeg >= s1.thetaMinDeg & thetaDeg <= s1.thetaMaxDeg;
            % Sector 2
            s2 = obj.acRegion.sector.s2;
            isInS2 = x >= s2.xmin & x <= s2.xmax & ...
                rYZ >= s2.rmin & rYZ <= s2.rmax & ...
                thetaDeg >= s2.thetaMinDeg & thetaDeg <= s2.thetaMaxDeg;
            % Sector 3
            s3 = obj.acRegion.sector.s3;
            isInS3 = x >= s3.xmin & x <= s3.xmax & ...
                rYZ >= s3.rmin & rYZ <= s3.rmax & ...
                thetaDeg >= s3.thetaMinDeg & thetaDeg <= s3.thetaMaxDeg;
            % Sector 4
            s4 = obj.acRegion.sector.s4;
            isInS4 = x >= s4.xmin & x <= s4.xmax & ...
                rYZ >= s4.rmin & rYZ <= s4.rmax;

            % global grid indices
            obj.acRegion.sector.s1.gridIdx = insidePointsIdx(isInS1);
            obj.acRegion.sector.s2.gridIdx = insidePointsIdx(isInS2);
            obj.acRegion.sector.s3.gridIdx = insidePointsIdx(isInS3);
            obj.acRegion.sector.s4.gridIdx = insidePointsIdx(isInS4);

            % vertices
            xv = vertices(:,1);
            yv = vertices(:,2);
            zv = vertices(:,3);
            rYZv = sqrt((yv - y0).^2 + (zv - z0).^2);
            thetaDegV = atan2d(yv - y0, zv - z0);

            % Sector 1
            isInS1v = xv >= s1.xmin & xv <= s1.xmax & ...
                rYZv >= s1.rmin & rYZv <= s1.rmax & ...
                thetaDegV >= s1.thetaMinDeg & thetaDegV <= s1.thetaMaxDeg;
            % Sector 2
            isInS2v = xv >= s2.xmin & xv <= s2.xmax & ...
                rYZv >= s2.rmin & rYZv <= s2.rmax & ...
                thetaDegV >= s2.thetaMinDeg & thetaDegV <= s2.thetaMaxDeg;
            % Sector 3
            isInS3v = xv >= s3.xmin & xv <= s3.xmax & ...
                rYZv >= s3.rmin & rYZv <= s3.rmax & ...
                thetaDegV >= s3.thetaMinDeg & thetaDegV <= s3.thetaMaxDeg;
            % Sector 4
            isInS4v = xv >= s4.xmin & xv <= s4.xmax & ...
                rYZv >= s4.rmin & rYZv <= s4.rmax;

            % global vertices indices
            obj.acRegion.sector.s1.vertexIdx = find(isInS1v);
            obj.acRegion.sector.s2.vertexIdx = find(isInS2v);
            obj.acRegion.sector.s3.vertexIdx = find(isInS3v);
            obj.acRegion.sector.s4.vertexIdx = find(isInS4v);

            disp('grid points / vertices assigned to acetabulum sectors')

        end

        %% Acetabulum region: defect points inside sectors
        function obj = acSectorDefect(obj, isPointInsideDefect, refSector, voxelSize, pelvisNum)

            % Reference sectors of pelvis
            s1.gridIdx = refSector.s1.gridIdx;
            s2.gridIdx = refSector.s2.gridIdx;
            s3.gridIdx = refSector.s3.gridIdx;
            s4.gridIdx = refSector.s4.gridIdx;

            % Number of reference grid points per sector
            nRefGridS1 = numel(s1.gridIdx);
            nRefGridS2 = numel(s2.gridIdx);
            nRefGridS3 = numel(s3.gridIdx);
            nRefGridS4 = numel(s4.gridIdx);

            % Voxel size [mm]
            voxelVolume = prod(voxelSize);

            % Global indices of defect grid points
            idxDefectGrid = find(isPointInsideDefect);

            % Total number of defect grid points
            nDefectGrid_total = numel(idxDefectGrid);

            % Defect grid points assigned to each sector
            idxDefectGridS1 = intersect(idxDefectGrid, s1.gridIdx);
            idxDefectGridS2 = intersect(idxDefectGrid, s2.gridIdx);
            idxDefectGridS3 = intersect(idxDefectGrid, s3.gridIdx);
            idxDefectGridS4 = intersect(idxDefectGrid, s4.gridIdx);

            % Defect grid points outside all sectors
            idxDefectGridInSectors = [idxDefectGridS1; idxDefectGridS2; idxDefectGridS3; idxDefectGridS4];
            idxDefectGridOutside = setdiff(idxDefectGrid, idxDefectGridInSectors);

            % Number of defect grid points per sector
            nDefectGridS1 = numel(idxDefectGridS1);
            nDefectGridS2 = numel(idxDefectGridS2);
            nDefectGridS3 = numel(idxDefectGridS3);
            nDefectGridS4 = numel(idxDefectGridS4);
            nDefectGridOutside = numel(idxDefectGridOutside);

            % Proportion of total defect grid points
            percentDefectGridS1 = 100 * nDefectGridS1 / nDefectGrid_total;
            percentDefectGridS2 = 100 * nDefectGridS2 / nDefectGrid_total;
            percentDefectGridS3 = 100 * nDefectGridS3 / nDefectGrid_total;
            percentDefectGridS4 = 100 * nDefectGridS4 / nDefectGrid_total;
            percentDefectGridOutside = 100 * nDefectGridOutside / nDefectGrid_total;

            % Loss percent relative to reference pelvis grid points inside each sector
            percentSectorGridS1 = 100 * nDefectGridS1 / nRefGridS1;
            percentSectorGridS2 = 100 * nDefectGridS2 / nRefGridS2;
            percentSectorGridS3 = 100 * nDefectGridS3 / nRefGridS3;
            percentSectorGridS4 = 100 * nDefectGridS4 / nRefGridS4;

            % Save
            obj.acRegion.sector.totalDefect.defectGridIdx = idxDefectGrid; % total defect, not all sub-cuboids (defect can also be outside sectors)
            obj.acRegion.sector.totalDefect.nGridPoints = nDefectGrid_total;
            obj.acRegion.sector.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume;

            obj.acRegion.sector.s1.defectGridIdx = idxDefectGridS1; % grid points (idx) in sector
            obj.acRegion.sector.s1.nGridPoints = nDefectGridS1; % number of grid points in sector
            obj.acRegion.sector.s1.gridLossPercentDefect = percentDefectGridS1; % volume loss (percent) per defect
            obj.acRegion.sector.s1.gridLossPercentSector = percentSectorGridS1; % volume loss (percent) per sector
            obj.acRegion.sector.s1.gridVolumeLoss = nDefectGridS1 * voxelVolume; % volume loss (mm3)

            obj.acRegion.sector.s2.defectGridIdx = idxDefectGridS2;
            obj.acRegion.sector.s2.nGridPoints = nDefectGridS2;
            obj.acRegion.sector.s2.gridLossPercentDefect = percentDefectGridS2;
            obj.acRegion.sector.s2.gridLossPercentSector = percentSectorGridS2;
            obj.acRegion.sector.s2.gridVolumeLoss = nDefectGridS2 * voxelVolume;

            obj.acRegion.sector.s3.defectGridIdx = idxDefectGridS3;
            obj.acRegion.sector.s3.nGridPoints = nDefectGridS3;
            obj.acRegion.sector.s3.gridLossPercentDefect = percentDefectGridS3;
            obj.acRegion.sector.s3.gridLossPercentSector = percentSectorGridS3;
            obj.acRegion.sector.s3.gridVolumeLoss = nDefectGridS3 * voxelVolume;

            obj.acRegion.sector.s4.defectGridIdx = idxDefectGridS4;
            obj.acRegion.sector.s4.nGridPoints = nDefectGridS4;
            obj.acRegion.sector.s4.gridLossPercentDefect = percentDefectGridS4;
            obj.acRegion.sector.s4.gridLossPercentSector = percentSectorGridS4;
            obj.acRegion.sector.s4.gridVolumeLoss = nDefectGridS4 * voxelVolume;

            obj.acRegion.sector.outside.defectGridIdx = idxDefectGridOutside;
            obj.acRegion.sector.outside.nGridPoints = nDefectGridOutside;
            obj.acRegion.sector.outside.gridLossPercentDefect = percentDefectGridOutside;
            obj.acRegion.sector.outside.gridVolumeLoss = nDefectGridOutside * voxelVolume;

            disp(['defect grid points assigned to acetabulum sectors: pelvis defect ', num2str(pelvisNum)]);

        end

%% Function to plot sectors
        function h = plotSectorSurface(obj, sector, faceColor, faceAlpha, nTheta)

            if nargin < 5
                nTheta = 80;
            end
            c = sector.centre;
            y0 = c(2);
            z0 = c(3);
            xmin = sector.xmin;
            xmax = sector.xmax;
            rmin = sector.rmin;
            rmax = sector.rmax;
            th1 = deg2rad(sector.thetaMinDeg);
            th2 = deg2rad(sector.thetaMaxDeg);
            theta = linspace(th1, th2, nTheta);
            % Outer arc
            Yout = y0 + rmax * sin(theta);
            Zout = z0 + rmax * cos(theta);
            % Inner arc
            Yin = y0 + rmin * sin(theta);
            Zin = z0 + rmin * cos(theta);
            % Outer cylindrical surface -> use this handle for legend
            Xo = [xmin * ones(1,nTheta); xmax * ones(1,nTheta)];
            Yo = [Yout; Yout];
            Zo = [Zout; Zout];
            h = surf(Xo, Yo, Zo, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
                'EdgeColor', faceColor, 'LineWidth', 1.2); %
            % Inner cylindrical surface
            if rmin > 0
                Xi = [xmin * ones(1,nTheta); xmax * ones(1,nTheta)];
                Yi = [Yin; Yin];
                Zi = [Zin; Zin];
                surf(Xi, Yi, Zi, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
                    'EdgeColor', faceColor, 'LineWidth', 1.2);
            end
            % Radial side surfaces
            for k = [1 nTheta]
                Xr = [xmin xmax; xmin xmax];
                Yr = [Yin(k) Yin(k); Yout(k) Yout(k)];
                Zr = [Zin(k) Zin(k); Zout(k) Zout(k)];
                surf(Xr, Yr, Zr, 'FaceColor', faceColor, 'FaceAlpha', faceAlpha, ...
                    'EdgeColor', faceColor, 'LineWidth', 1.2);
            end
            % Front / back faces
            if rmin > 0
                Yf = [Yout fliplr(Yin)];
                Zf = [Zout fliplr(Zin)];
            else
                Yf = [Yout y0];
                Zf = [Zout z0];
            end
            Xf = xmin * ones(size(Yf));
            patch(Xf, Yf, Zf, faceColor, 'FaceAlpha', faceAlpha, ...
                'EdgeColor', faceColor, 'LineWidth', 1.2);
            Xb = xmax * ones(size(Yf));
            patch(Xb, Yf, Zf, faceColor, 'FaceAlpha', faceAlpha, ...
                'EdgeColor', faceColor, 'LineWidth', 1.2);
        end

        %% Acetabulum region: hybrid = cuboids with central sector
        function obj = acHybridInside(obj, cuboidRef, sectorRef, insidePoints, insidePointsIdx, vertices)

            % Copy geometry/limits from cuboid
            obj.acRegion.hybrid.limits = cuboidRef.all.limits;

            % Keep cuboid geometry for display
            obj.acRegion.hybrid.mc1.vertices = cuboidRef.c1.vertices;
            obj.acRegion.hybrid.mc1.faces    = cuboidRef.c1.faces;
            obj.acRegion.hybrid.mc1.limits   = cuboidRef.c1.limits;

            obj.acRegion.hybrid.mc2.vertices = cuboidRef.c2.vertices;
            obj.acRegion.hybrid.mc2.faces    = cuboidRef.c2.faces;
            obj.acRegion.hybrid.mc2.limits   = cuboidRef.c2.limits;

            obj.acRegion.hybrid.mc3.vertices = cuboidRef.c3.vertices;
            obj.acRegion.hybrid.mc3.faces    = cuboidRef.c3.faces;
            obj.acRegion.hybrid.mc3.limits   = cuboidRef.c3.limits;

            obj.acRegion.hybrid.mc4.vertices = cuboidRef.c4.vertices;
            obj.acRegion.hybrid.mc4.faces    = cuboidRef.c4.faces;
            obj.acRegion.hybrid.mc4.limits   = cuboidRef.c4.limits;

            % Central sector mc5 geometry
            obj.acRegion.hybrid.mc5 = sectorRef.s4;
            % Adapt central sector depth to cuboid depth, e.g. 2.2 * acentreR
            obj.acRegion.hybrid.mc5.xmin = cuboidRef.c1.limits.xmin;
            obj.acRegion.hybrid.mc5.xmax = cuboidRef.c1.limits.xmax;

            % s4 indices to remove from cuboids
            %mc5GridIdx   = sectorRef.s4.gridIdx;
            %mc5VertexIdx = sectorRef.s4.vertexIdx;
            % Recalculate mc5 indices with adapted x-depth
            centre = sectorRef.limits.centre;
            x0 = centre(1);
            y0 = centre(2);
            z0 = centre(3);
            mc5 = obj.acRegion.hybrid.mc5;
            % Grid points
            x = insidePoints(:,1);
            y = insidePoints(:,2);
            z = insidePoints(:,3);
            rYZ = sqrt((y - y0).^2 + (z - z0).^2);
            isInMc5 = x >= mc5.xmin & x <= mc5.xmax & ...
                rYZ >= mc5.rmin & rYZ <= mc5.rmax;
            mc5GridIdx = insidePointsIdx(isInMc5);
            % Vertices
            xv = vertices(:,1);
            yv = vertices(:,2);
            zv = vertices(:,3);
            rYZv = sqrt((yv - y0).^2 + (zv - z0).^2);
            isInMc5v = xv >= mc5.xmin & xv <= mc5.xmax & ...
                rYZv >= mc5.rmin & rYZv <= mc5.rmax;
            mc5VertexIdx = find(isInMc5v);

            % Grid indices: cuboid minus mc5/s4
            obj.acRegion.hybrid.mc1.gridIdx = setdiff(cuboidRef.c1.gridIdx, mc5GridIdx);
            obj.acRegion.hybrid.mc2.gridIdx = setdiff(cuboidRef.c2.gridIdx, mc5GridIdx);
            obj.acRegion.hybrid.mc3.gridIdx = setdiff(cuboidRef.c3.gridIdx, mc5GridIdx);
            obj.acRegion.hybrid.mc4.gridIdx = setdiff(cuboidRef.c4.gridIdx, mc5GridIdx);

            % Vertex indices: cuboid minus mc5/s4
            obj.acRegion.hybrid.mc1.vertexIdx = setdiff(cuboidRef.c1.vertexIdx, mc5VertexIdx);
            obj.acRegion.hybrid.mc2.vertexIdx = setdiff(cuboidRef.c2.vertexIdx, mc5VertexIdx);
            obj.acRegion.hybrid.mc3.vertexIdx = setdiff(cuboidRef.c3.vertexIdx, mc5VertexIdx);
            obj.acRegion.hybrid.mc4.vertexIdx = setdiff(cuboidRef.c4.vertexIdx, mc5VertexIdx);

            % Central sector remains as own region
            obj.acRegion.hybrid.mc5.gridIdx   = mc5GridIdx;
            obj.acRegion.hybrid.mc5.vertexIdx = mc5VertexIdx;

            disp('grid points / vertices assigned to acetabulum hybrid regions')
        end

        %% Acetabulum region: defect points inside hybrid regions
        function obj = acHybridDefect(obj, isPointInsideDefect, refHybrid, voxelSize, pelvisNum)

            % Hybrid regions of reference pelvis
            mc1.gridIdx = refHybrid.mc1.gridIdx;
            mc2.gridIdx = refHybrid.mc2.gridIdx;
            mc3.gridIdx = refHybrid.mc3.gridIdx;
            mc4.gridIdx = refHybrid.mc4.gridIdx;
            mc5.gridIdx = refHybrid.mc5.gridIdx;

            % Number of reference grid points per region
            nRefMC1 = numel(mc1.gridIdx);
            nRefMC2 = numel(mc2.gridIdx);
            nRefMC3 = numel(mc3.gridIdx);
            nRefMC4 = numel(mc4.gridIdx);
            nRefMC5 = numel(mc5.gridIdx); %%%

            % Voxel size [mm]
            voxelVolume = prod(voxelSize);

            % Global indices of defect grid points
            idxDefectGrid = find(isPointInsideDefect);

            % Total number of defect grid points
            nDefectGrid_total = numel(idxDefectGrid);

            % Defect grid points assigned to each hybrid region
            idxDefectMC1 = intersect(idxDefectGrid, mc1.gridIdx);
            idxDefectMC2 = intersect(idxDefectGrid, mc2.gridIdx);
            idxDefectMC3 = intersect(idxDefectGrid, mc3.gridIdx);
            idxDefectMC4 = intersect(idxDefectGrid, mc4.gridIdx);
            idxDefectMC5 = intersect(idxDefectGrid, mc5.gridIdx);

            % Defect grid points outside all hybrid regions
            idxDefectAssigned = [idxDefectMC1; idxDefectMC2; idxDefectMC3; idxDefectMC4; idxDefectMC5];
            idxDefectOutside = setdiff(idxDefectGrid, idxDefectAssigned);

            % Number of defect grid points per region
            nDefectMC1 = numel(idxDefectMC1);
            nDefectMC2 = numel(idxDefectMC2);
            nDefectMC3 = numel(idxDefectMC3);
            nDefectMC4 = numel(idxDefectMC4);
            nDefectMC5 = numel(idxDefectMC5);
            nDefectOutside = numel(idxDefectOutside);

            % Proportion of total defect grid points
            percentDefectMC1 = 100 * nDefectMC1 / nDefectGrid_total;
            percentDefectMC2 = 100 * nDefectMC2 / nDefectGrid_total;
            percentDefectMC3 = 100 * nDefectMC3 / nDefectGrid_total;
            percentDefectMC4 = 100 * nDefectMC4 / nDefectGrid_total;
            percentDefectMC5 = 100 * nDefectMC5 / nDefectGrid_total;
            percentDefectOutside = 100 * nDefectOutside / nDefectGrid_total;

            % Loss percent relative to reference pelvis grid points inside each region
            percentRegionMC1 = 100 * nDefectMC1 / nRefMC1;
            percentRegionMC2 = 100 * nDefectMC2 / nRefMC2;
            percentRegionMC3 = 100 * nDefectMC3 / nRefMC3;
            percentRegionMC4 = 100 * nDefectMC4 / nRefMC4;
            percentRegionMC5 = 100 * nDefectMC5 / nRefMC5;

            % Save
            obj.acRegion.hybrid.totalDefect.defectGridIdx = idxDefectGrid;
            obj.acRegion.hybrid.totalDefect.nGridPoints = nDefectGrid_total;
            obj.acRegion.hybrid.totalDefect.gridVolumeLoss = nDefectGrid_total * voxelVolume;

            obj.acRegion.hybrid.mc1.defectGridIdx = idxDefectMC1;
            obj.acRegion.hybrid.mc1.nGridPoints = nDefectMC1;
            obj.acRegion.hybrid.mc1.gridLossPercentDefect = percentDefectMC1;
            obj.acRegion.hybrid.mc1.gridLossPercentRegion = percentRegionMC1;
            obj.acRegion.hybrid.mc1.gridVolumeLoss = nDefectMC1 * voxelVolume;

            obj.acRegion.hybrid.mc2.defectGridIdx = idxDefectMC2;
            obj.acRegion.hybrid.mc2.nGridPoints = nDefectMC2;
            obj.acRegion.hybrid.mc2.gridLossPercentDefect = percentDefectMC2;
            obj.acRegion.hybrid.mc2.gridLossPercentRegion = percentRegionMC2;
            obj.acRegion.hybrid.mc2.gridVolumeLoss = nDefectMC2 * voxelVolume;

            obj.acRegion.hybrid.mc3.defectGridIdx = idxDefectMC3;
            obj.acRegion.hybrid.mc3.nGridPoints = nDefectMC3;
            obj.acRegion.hybrid.mc3.gridLossPercentDefect = percentDefectMC3;
            obj.acRegion.hybrid.mc3.gridLossPercentRegion = percentRegionMC3;
            obj.acRegion.hybrid.mc3.gridVolumeLoss = nDefectMC3 * voxelVolume;

            obj.acRegion.hybrid.mc4.defectGridIdx = idxDefectMC4;
            obj.acRegion.hybrid.mc4.nGridPoints = nDefectMC4;
            obj.acRegion.hybrid.mc4.gridLossPercentDefect = percentDefectMC4;
            obj.acRegion.hybrid.mc4.gridLossPercentRegion = percentRegionMC4;
            obj.acRegion.hybrid.mc4.gridVolumeLoss = nDefectMC4 * voxelVolume;

            obj.acRegion.hybrid.mc5.defectGridIdx = idxDefectMC5;
            obj.acRegion.hybrid.mc5.nGridPoints = nDefectMC5;
            obj.acRegion.hybrid.mc5.gridLossPercentDefect = percentDefectMC5;
            obj.acRegion.hybrid.mc5.gridLossPercentRegion = percentRegionMC5;
            obj.acRegion.hybrid.mc5.gridVolumeLoss = nDefectMC5 * voxelVolume;

            obj.acRegion.hybrid.outside.defectGridIdx = idxDefectOutside;
            obj.acRegion.hybrid.outside.nGridPoints = nDefectOutside;
            obj.acRegion.hybrid.outside.gridLossPercentDefect = percentDefectOutside;
            obj.acRegion.hybrid.outside.gridVolumeLoss = nDefectOutside * voxelVolume;

            disp(['defect grid points assigned to hybrid acetabulum regions: pelvis defect ', num2str(pelvisNum)]);
        end

        %% Display radar plot per defect
       function obj = plotRadar(obj, values, plotTitle, labels)
           plotColor = [0, 101/255, 189/255]; % TUMcolors.blue300

           values = values(:)';   % row vector
           labels = labels(:)';   % row vector
           nAxes = numel(values);

           if numel(labels) ~= nAxes
               error('plotRadar:SizeMismatch', ...
                   'Number of labels (%d) must match number of values (%d).', ...
                   numel(labels), nAxes);
           end

           rMax = 100;
           rTicks = [0 25 50 75 100];

           % Dynamic axis positions:
           % first axis at top, then clockwise
           thetaDeg = 90 - (0:nAxes-1) * 360/nAxes;
           theta = deg2rad(thetaDeg);

           % Close polygon
           thetaClosed = [theta theta(1)];
           valuesClosed = [values values(1)];

           % Cartesian coordinates
           xPoly = valuesClosed .* cos(thetaClosed);
           yPoly = valuesClosed .* sin(thetaClosed);

           hold on
           axis equal
           axis off

           % Grid polygons
           for r = rTicks(2:end)
               xg = r * cos(thetaClosed);
               yg = r * sin(thetaClosed);
               plot(xg, yg, '-', ...
                   'Color', [0.75 0.75 0.75], ...
                   'LineWidth', 1);
           end

           % Axes
           for k = 1:nAxes
               plot([0 rMax*cos(theta(k))], ...
                   [0 rMax*sin(theta(k))], '-', ...
                   'Color', [0.6 0.6 0.6], ...
                   'LineWidth', 1);
           end

           % Filled polygon
           patch(xPoly, yPoly, plotColor, ...
               'FaceAlpha', 0.25, ...
               'EdgeColor', plotColor, ...
               'LineWidth', 2);

           % Markers
           plot(xPoly, yPoly, 'o', ...
               'Color', plotColor, ...
               'MarkerFaceColor', plotColor, ...
               'MarkerSize', 5);

           % Labels
           labelRadius = rMax * 1.12;
           for k = 1:nAxes
               xl = labelRadius * cos(theta(k));
               yl = labelRadius * sin(theta(k));

               text(xl, yl, labels{k}, ...
                   'HorizontalAlignment', 'center', ...
                   'VerticalAlignment', 'middle', ...
                   'FontWeight', 'bold', ...
                   'FontSize', 12);
           end

           % Tick labels near top
           for r = rTicks
               text(-8, r, num2str(r), ...
                   'HorizontalAlignment', 'right', ...
                   'VerticalAlignment', 'middle', ...
                   'FontWeight', 'bold', ...
                   'FontSize', 10, ...
                   'Color', [0.1 0.1 0.1]);
           end

           xlim([-1.25*rMax, 1.25*rMax]);
           ylim([-1.15*rMax, 1.15*rMax]);

           title(plotTitle, 'FontWeight', 'bold');
           hold off
       end
    end
end
