
function [C, T, N] = discrete_contour_mesh_patch(V, mode)
%% discrete_contour_mesh_patch : function to mesh one discrete 2D
% or 3D contour composed of sorted or disordered 3D points.
%
% Author and support : nicolas.douillet (at) free.fr, 2021.
%
%
% Syntax
%
% [C, T] = discrete_contour_mesh_patch(V);
% [C, T] = discrete_contour_mesh_patch(V, mode);
% [C, T, N] = discrete_contour_mesh_patch(V);
%
%
% Description
%
% [C, T] = discrete_contour_mesh_patch(V) meshes the discrete contour
% defined by the vertex set V and store the resulting triangulation in T
% after having sorted -default behaviour- V elements by their angle
% relatively to their isobarycentre, and whom result is stored in C.
%
% [C, T] = discrete_contour_mesh_patch(V, mode) uses the mode parameter to
% sort (mode = 'raw') or not (mode = 'sorted') the vertex set V.
%
% [C, T, N] = discrete_contour_mesh_patch(V) also computes and returns the
% vertex normals set N.
%
%
% See also : MESH, TRIMESH
%
%
% Input arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the point set, size(V) = [nb_vertices,3].
%       [ |  | | ]
%
% - mode : character string in the set {'raw','sorted'},
%          the type of contour considered. Case insensitive.
%
%
% Output argument
%
%       [ |  |  |]
% - C = [Cx Cy Cz], real matrix double, the output sorted point set, size(C) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [i1 i2 i3], positive integer matrix double, the resulting triangulation, size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - N : [Nx Ny Nz], real matrix double, the vertex normal vectors, size(N) = [nb_vertices,3].
%       [ |  |  |]
%% Body
tic;
assert(nargin > 0,'Not enought input argument.');
if nargin < 2
   
    mode = 'raw'; % default behaviour
    
end
if strcmpi(mode,'raw') % reorder / angular sort if necessary
    
    G = mean(V,1);
    U = V - repmat(G,[size(V,1),1]);
    ref_vect1 = repmat(U(1,:),[size(U,1),1]);
    bov = cross(U(1,:),U(2,:),2);
    angl = atan2(sign(dot(cross(ref_vect1,U,2),repmat(bov,[size(U,1),1]),2)).*sqrt(sum(cross(ref_vect1,U,2).^2,2)),dot(ref_vect1,U,2));
    [~,idx] = sort(angl);
    C = V(idx,:);    
    
elseif strcmpi(mode,'sorted')
    
    C = V;    
    
else
    
    error('Unrecognized mode.');
    
end
boundary = 1:size(C,1);
% Outward oriented normalized vertex normals
bound_nb_vtx = numel(boundary);
B = repelem(boundary',cat(2,1,2*ones(1,bound_nb_vtx-2),1));
B = reshape(B,[2,bound_nb_vtx-1])';
B = cat(1,B,[B(end,end) B(1,1)]);
Tf = C(B(:,2),:) - C(B(:,1),:);
Tf = Tf ./ sqrt(sum(Tf.^2,2));
Tb = -Tf;
N = -Tf - circshift(Tb,1,1);
G = mean(C,1);
orientation = sign(dot(N,C-repmat(G,[bound_nb_vtx,1]),2));
if ~isequal(orientation,ones(bound_nb_vtx,1))

    N = N.*orientation;

end
N = N ./ sqrt(sum(N.^2,2));
cross_prod = @(boundary_backward,boundary,boundary_forward) cross(C(boundary_forward,:)-C(boundary,:),C(boundary_backward,:)-C(boundary,:),2);
dot_prod   = @(boundary_backward,boundary,boundary_forward)   dot(C(boundary_forward,:)-C(boundary,:),C(boundary_backward,:)-C(boundary,:),2);
sgn = @(bov,cross_prod) sign(dot(cross_prod,repmat(bov,[size(cross_prod,1),1]),2));
edg_angle = @(sgn,cross_prod,dot_prod) atan2(sgn.*sqrt(sum(cross_prod.^2,2)),dot_prod);
% Initialization
T = zeros(0,3);
nb_added_tgl = 0;
bov = compute_boundary_orientation_vector(boundary,C);
while bound_nb_vtx > 2
    
    boundary_forward  = circshift(boundary,-1);
    boundary_backward = circshift(boundary,1);        
    
    c_prod = cross_prod(boundary_backward,boundary,boundary_forward);
    d_prod = dot_prod(boundary_backward,boundary,boundary_forward);
    sgnv = sgn(bov,c_prod);
    edg_angl_list = edg_angle(sgnv,c_prod,d_prod);    
    
    if isequal(sign(edg_angl_list),-ones(size(edg_angl_list,1),size(edg_angl_list,2)))
        
        boundary = fliplr(boundary);
        boundary_forward = circshift(boundary,-1);
        boundary_backward = circshift(boundary,1);
        
        c_prod = cross_prod(boundary_backward,boundary,boundary_forward);
        d_prod = dot_prod(boundary_backward,boundary,boundary_forward);
        sgnv = sgn(bov,c_prod);
        edg_angl_list = edg_angle(sgnv,c_prod,d_prod);        
        
    end          
    
    criterion = edg_angl_list;    
    min_pos = min(criterion(criterion > 0));
    
    if ~isempty(min_pos)
        
        min_angle_idx = find(criterion == min_pos(1,1));
        
    else
        
        min_angle_idx = 1; % special possible case of the last triangle; index of min doesn't matter
        
    end
    
    if min_angle_idx; min_angle_idx = min_angle_idx(1,1); end
    
    new_triangle = [boundary_forward(min_angle_idx), boundary(min_angle_idx), boundary_backward(min_angle_idx)];
    T = add_triangles(new_triangle,T,size(C,1));        
    
    nb_added_tgl = nb_added_tgl + 1;
    boundary(min_angle_idx) = [];
    bound_nb_vtx = numel(boundary);
    
end
end % discrete_contour_mesh_patch
% Subfunction
function bov = compute_boundary_orientation_vector(boundary,V)
nb_edg = numel(boundary); 
bov = cross(mean(V(boundary(1,2:floor(0.5*nb_edg)),:)-repmat(V(boundary(1,1),:),[floor(0.5*nb_edg)-1,1]),1),...
            mean(V(boundary(1,ceil(0.5*nb_edg):end),:)-repmat(V(boundary(1,1),:),[nb_edg-ceil(0.5*nb_edg)+1,1]),1),2);
end % compute_boundary_orientation_vector
% Subfunction

function T_out = add_triangles(T_set, T_in, nb_vtx)
%% add_triangles : function to add new triangles to the triangle set.
%
% Author & support : nicolas.douillet (at) free.fr, 2021.
%
%
% Input arguments
%
%           [|  |  | ]
% - T_set = [i1 i2 i3], positive integer matrix double, the triangle set to add, size(T_set) = [nb_new_triangles,3].
%           [|  |  | ]
%
%          [  |     |     |  ]
% - T_in = [i1_in i2_in i3_in], positive integer matrix double, the input triangulation, size(T_in) = [nb_input_triangles,3].
%          [  |     |     |  ]
%
% - nb_vtx : positive integer scalar double, the number of vertices in the input vertex set. nb_vtx = size(V_in,1).
%
%
% Output arguments
%
%           [  |      |      |   ]
% - T_out = [i1_out i2_out i3_out], positive integer matrix double, the output triangulation, size(T_out) = [nb_output_triangles,3],
%           [  |      |      |   ]
%
%           with nb_output_triangles = nb_input_triangles + nb_new_triangles.
%% Body
% tic;
% Check if vertices id are valid (0 < positive integers <= nb_vtx)
if T_set > 0 & T_set <= nb_vtx & isreal(T_set) & rem(T_set,1) == 0 & floor(T_set) == T_set
    
    % Check if some triangles of T_set already part of T_in
    dpl_tgl_idx = ismember(T_set,T_in,'rows');
    
    if ~nnz(dpl_tgl_idx) 
        
        T_out = cat(1,T_in,T_set);
        % fprintf('%d triangles added in %d seconds.',size(T_set,1),toc);
        
    else % if nnz(dpl_tgl_idx) 
       
        warning('One or more triangle from this set already exist in the current triangulation. Duplicated triangles have been ignored.\n');
        
        % Suppress duplicated triangles
        T_set = T_set(~dpl_tgl_idx,:);
        T_out = cat(1,T_in,T_set);        
        
    end
    
else
    
    error('Unable to perform triangles addition because T_set contains some invalid vertex indices. Vertex indices must be positive integers in the range |[1; %d]|.\n',nb_vtx);
    
end
end % add_triangles

