%% Surface reconstruction of point cloud with crust function

% Input:    pelvisNum: Numeric identifier used only for logging
%           importData: with comVertices as point cloud baseline
%           addPoints: Additional points to augment the cloud

% Output:   shrink: Struct with fields added under shrink.crust

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = crust(pelvisNum, importData, addPoints)

verticesImport = importData.comVertices;
vertices = [verticesImport; addPoints];

% Matlab File Exchange: Surface Reconstruction From Scattered Points Cloud
% https://de.mathworks.com/matlabcentral/fileexchange/63730-surface-reconstruction-from-scattered-points-cloud
[faces,norm] = MyRobustCrust(vertices);
shrink.crust.faces = faces;
shrink.crust.normals = norm;
shrink.crust.vertices = vertices;

figure(1)
hold on
title('Output Triangulation','fontsize',14)
axis equal
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
view(3);
axis vis3d

disp(['surface reconstruction: pelvis defect ',num2str(pelvisNum)])
end
