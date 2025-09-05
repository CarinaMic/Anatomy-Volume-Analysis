%% Shrink Wrap with boundary function over vertices

% Input:    pelvisNum: Numeric identifier used only for logging
%           acentre: coordinates of acentre
%           vertices: input points to be wrapped

% Output:   shrink
%           shrink.bound.faces: triangle indices of the boundary surface
%           shrink.bound.vertices: combined point set [vertices; acentre]

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function [shrink] = shrinkBound(pelvisNum,acentre,vertices)

% Add the additional point acentre
combinedVertices = vertices;
combinedVertices = [combinedVertices; acentre];

% ShrinkWrap with boundary function
% https://de.mathworks.com/help/matlab/ref/boundary.html#buh3c7k-1-s
% Shrink factor is between 0 and 1. 0 gives the convex hull,
% and 1 gives a compact boundary that envelops the points. (default 0.5)
shrink.bound.faces = boundary(combinedVertices,0.8); % define shrink factor
shrink.bound.vertices = combinedVertices;

disp(['shrinked mesh (boundary): pelvis defect ',num2str(pelvisNum)])

end