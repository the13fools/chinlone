function [Mesh] = cylinder_generation(n_rad, n_sample, angle)
% Compute a cylinder with resolution n_rad x n_sample
thetas_rad = 2*pi*(0:(n_rad-1))/n_rad;
d = 2*sin(pi/n_rad);
VertexPoss_x = (kron(0:(n_sample-1), ones(1,n_rad))-(n_sample-1)/2)*d;
VertexPoss_y = kron(ones(1,n_sample), cos(thetas_rad));
VertexPoss_z = kron(ones(1,n_sample), sin(thetas_rad));
Mesh.vertexPoss = [VertexPoss_x; VertexPoss_y; VertexPoss_z];
%
IDX = reshape(1:(n_rad*n_sample), [n_rad, n_sample]);
IDX = [IDX;IDX(1,:)];
%
ids1 = reshape(IDX(1:n_rad, 1:(n_sample-1)), [1, n_rad*(n_sample-1)]);
ids2 = reshape(IDX(2:(n_rad+1), 1:(n_sample-1)), [1, n_rad*(n_sample-1)]);
ids3 = reshape(IDX(2:(n_rad+1), 2:n_sample), [1, n_rad*(n_sample-1)]);
ids4 = reshape(IDX(1:n_rad, 2:n_sample), [1, n_rad*(n_sample-1)]);
%
Mesh.faceVIds = reshape([ids1;ids2;ids4;ids2;ids3;ids4], [3, 2*n_rad*(n_sample-1)]);
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
numF = 2*n_rad*(n_sample-1);
axis_x = kron([1,0,0]', ones(1,numF));
inner = sum(axis_x.*Mesh.faceNors);
axis_x = axis_x - Mesh.faceNors.*(ones(3,1)*inner);
sqrtNorm = sqrt(sum(axis_x.*axis_x));
axis_x = axis_x./(ones(3,1)*sqrtNorm);
axis_y = cross(Mesh.faceNors, axis_x);
faceVecs = cos(angle)*axis_x + sin(angle)*axis_y;
Mesh.faceVecs = faceVecs;