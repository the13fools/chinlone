function [Mesh] = parse_sample_mesh()

a = load_obj("tet.obj");
M_ret.vertex = a.vertex(:,[1:3])';
M_ret.faces = zeros(3, size(a.faces, 1));

for i = 1:size(a.faces, 1)
    c = a.faces(i);
    open_cell = c{1};
    M_ret.faces(:,i) = open_cell(:,1);
end


Mesh.vertexPoss = M_ret.vertex;

% Mesh.faceVIds = reshape([ids1;ids2;ids4;ids2;ids3;ids4], [3, 2*n_rad*(n_sample-1)]);
Mesh.faceVIds = M_ret.faces;
%
p1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
p2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
p3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
%
e12 = p1 - p2;
e13 = p1 - p3;
Mesh.faceNors = cross(e12, e13);
sqrNorms = sqrt(sum(Mesh.faceNors.*Mesh.faceNors));
Mesh.faceNors = Mesh.faceNors./(kron(ones(3,1), sqrNorms));


Mesh.faceVecs = importdata("tet_single_axis_aligned_field_0.txt")';
% Mesh.faceVecs = faceVecs;
