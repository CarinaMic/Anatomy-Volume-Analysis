%% Check for orientation of the triangulated patch (user-input)

% Input:    pelvisNum: Numeric identifier used only for logging
%           edgeNum: Index of the processed edge loop
%           importData: Struct with triangulated mesh (vertices,faces)
%           patchVertices: vertices of the triangulated patch to inspect
%           patchFaces: faces of the patch
%           patchNormals: per-face normals of the patch
%           patchCentreFaces: per-face centroids of the patch
%           type: 'all','acentre','acetabulum'

% Output:   patches: struct with user decision
%               (0 = keep orientation, 1 = flip requested)

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [patches] = flipNormals(pelvisNum,edgeNum,importData,patchVertices,patchFaces,patchNormals,patchCentreFaces,type)

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
title('pelvis ' num2str(pelvisNum));
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

patches.(type).flip{edgeNum,1} = 0;

while confirmButton.Value == 0
    %drawnow;
    pause(0.05);
end

% Plot is closed
close(fig)

% Functions of push buttons
    function flipCallback(src,~)
        patches.(type).flip{edgeNum,1} = src.Value;
    end

disp(['checked orientation of the triangulation (patch): pelvis defect ',num2str(pelvisNum)])

end