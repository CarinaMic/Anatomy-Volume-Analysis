%% Shrink Wrap with alphaShape (with user-input)

% Input:    pelvisNum: Numeric identifier used only for logging
%           patches: Container describing patches
%           edge: Edge/loop geometry (for visualization only)
%           type: 'all','acentre','acetabulum'
%           aRadius: Initial alpha radius %
%           numData: Number of rows in the original defect vertex array
%           importData:  Input geometry (vertices,comVertices,faces,acentre)
%           addPoints: Extra points appended to the point set

% Output:   shrink: Results struct with dynamic field named by 'type'

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shrinkAlphaUI(patches,edge,pelvisNum,type,aRadius,numData,importData,addPoints)

if isfield(importData, 'comVertices') && isfield(importData, 'comFaces')
    verticesImport = importData.comVertices;
    facesImport = importData.comFaces;
else
    verticesImport = importData.vertices;
    facesImport = importData.faces;
end
if nargin < 8
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
    for k = 1:patches.all.numComponents
        plot3(edge.verticesLoops{k}(:,1),edge.verticesLoops{k}(:,2),...
            edge.verticesLoops{k}(:,3),...
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
    if nargin >= 8
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
shrink.(type).radius = aRadius;
shrink.(type).alpha = alphaShape(allVerticesPoints,aRadius);
% ShrinkWrap: vertices -> not all vertices are used
shrinkFaces = boundaryFacets(shrink.(type).alpha);
shrink.(type).faces = shrinkFaces;
shrink.(type).allVerticesPoints = allVerticesPoints; % all vertices + acentre
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