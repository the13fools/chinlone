function [func] = func_from_vectorfield(Mesh, faceVecs)
%
numV = size(Mesh.vertexPoss, 2);
numF = size(Mesh.faceVIds, 2);
% Rotate the vectors
faceVecs = cross(Mesh.faceNors, faceVecs);
% Extract edges and compute the two (sometimes one? adjacent faces
[edges, adjFaces] = mesh_topology(Mesh);

% Estimate (not globally normalized up to a scale) the local scales per
% face
[rows, cols, difVec] = optimize_edge_difference(Mesh.vertexPoss,...
    Mesh.faceVIds, edges, adjFaces, faceVecs);

fprintf('Perform parameter search...\n');
[eigenVec, alpha_opt] = bin_search(rows, cols, difVec,...
    exp(-2), exp(3), 51);
fprintf('Alpha_opt = %f.\n', alpha_opt);

[eigenVec, alpha_opt] = bin_search(rows, cols, difVec,...
    alpha_opt*exp(-0.1), alpha_opt*exp(0.1), 51);
fprintf('Alpha_opt = %f.\n', alpha_opt);

[eigenVec, alpha_opt] = bin_search(rows, cols, difVec,...
    alpha_opt*exp(-0.002), alpha_opt*exp(0.002), 21);
fprintf('Alpha_opt = %f.\n', alpha_opt);

% Write the output
a = eigenVec(1:2:(2*numV));
b = eigenVec(2:2:(2*numV));
s = sqrt(a.*a + b.*b);
costheta = a./s;
sintheta = b./s;
func = acos(costheta);
ids = find(sintheta < 0);
func(ids) = -func(ids);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform bin search to find the optimal alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eigenVec, alpha_opt] = bin_search(rows, cols, difVec,...
    alpha_min, alpha_max, num)
%
t = 0:(1/(num-1)):1;
alphas = exp(log(alpha_min)*(1-t)+log(alpha_max)*t);
residuals = zeros(1, length(alphas));
for i = 1:length(alphas)
    [d, A] = connection_adjacency(rows, cols, difVec, alphas(i));
    residuals(i) = objective_function(d, A);
end
[s,id] = min(residuals);
alpha_opt = alphas(id);

[d, A] = connection_adjacency(rows, cols, difVec, alpha_opt);
numV = length(d);
L = sparse(1:(2*numV), 1:(2*numV),kron(d, ones(1,2))) - A;
[eigenVec, e_cur] = eigs(L, 1, 1e-10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% curvature fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [k] = curvature_fitting(x, y)
A = [x'.*x', x', ones(3,1)];
coeff = inv(A'*A)*(y*A)';
k = coeff(1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the connection adjacency matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d, A] = connection_adjacency(rows, cols, difVec, alpha)
% A function that generates a connection laplacian from the relative
% rotations angles along the mesh edges
% d: a vector that encodes the vertex degrees
% A: a generalized adjacency matrix that encodes pair-wise 
numE = length(rows);
numV = max(rows);
TP = sparse(rows, cols, ones(1,numE), numV, numV);
d = sum(TP + TP');
rowsA(1,:) = 2*kron(rows, ones(1,2))-1;
rowsA(2,:) = rowsA(1,:)  + 1;
colsA = zeros(2, 2*numE);
colsA(:,1:2:(2*numE)) = 2*kron(cols, ones(2,1))-1;
colsA(:,2:2:(2*numE)) = colsA(:,1:2:(2*numE)) + 1;
valsA = zeros(2, 2*numE);
valsA(1,1:2:(2*numE)) = cos(difVec*alpha);
valsA(2,1:2:(2*numE)) = sin(difVec*alpha);
valsA(2,2:2:(2*numE)) = valsA(1,1:2:(2*numE));
valsA(1,2:2:(2*numE)) = -valsA(2,1:2:(2*numE));
A = sparse(rowsA, colsA, valsA, 2*numV, 2*numV);
A = A+ A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the deriavtive of the connection Laplacian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L_prime] = connection_laplacian_derivative(rows, cols, difVec, alpha)
% A function that generates a connection laplacian from the relative
% rotations angles along the mesh edges
% d: a vector that encodes the vertex degrees
% A: a generalized adjacency matrix that encodes pair-wise 
numE = length(rows);
numV = max(rows);
rowsA(1,:) = 2*kron(rows, ones(1,2))-1;
rowsA(2,:) = rowsA(1,:)  + 1;
colsA = zeros(2, 2*numE);
colsA(:,1:2:(2*numE)) = 2*kron(cols, ones(2,1))-1;
colsA(:,2:2:(2*numE)) = colsA(:,1:2:(2*numE)) + 1;
valsA = zeros(2, 2*numE);
valsA(1,1:2:(2*numE)) = sin(difVec*alpha);
valsA(2,1:2:(2*numE)) = -cos(difVec*alpha);
valsA(2,2:2:(2*numE)) = valsA(1,1:2:(2*numE));
valsA(1,2:2:(2*numE)) = -valsA(2,1:2:(2*numE));
A = sparse(rowsA, colsA, valsA, 2*numV, 2*numV);
L_prime = A+ A';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the objective value of the connection Laplacian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eigenVal] = objective_function(d, A)
%
numV = length(d);
L = sparse(1:(2*numV), 1:(2*numV),kron(d, ones(1,2))) - A;
%
eigenVal = eigs(L, 1, 1e-10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimize the edge differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rows, cols, difVec] = optimize_edge_difference(vertexPoss,...
    faceVIds, edges, adjFaces, faceVecs)
% Estimate s_f
numF = size(faceVecs, 2);
validIds = find(min(adjFaces) > 0);
edgeVecs = vertexPoss(:,edges(1,validIds))-...
    vertexPoss(:,edges(2,validIds));
de_left = sum(edgeVecs.*faceVecs(:, adjFaces(1,validIds)));
de_right = sum(edgeVecs.*faceVecs(:, adjFaces(2,validIds)));
%Perform eigen-decomposition to find the global scale
rowsJ = ones(2,1)*(1:length(validIds));
valsJ = [de_left; -de_right];
J = sparse(rowsJ, adjFaces(:,validIds), double(valsJ),...
    length(validIds), numF);
L = J'*J;
[eigenVec, eigenVal] = eigs(L, 1, -1e-10);
if sum(eigenVec) < 0
    eigenVec = -eigenVec;
end
eigenVec = eigenVec';
relativeScales = eigenVec/mean(eigenVec);
%
faceVecs = faceVecs.*(ones(3,1)*relativeScales);
p1 = vertexPoss(:, faceVIds(1,:));
p2 = vertexPoss(:, faceVIds(2,:));
p3 = vertexPoss(:, faceVIds(3,:));
e12 = p1 - p2;
e23 = p2 - p3;
e31 = p3 - p1;
d12 = sum(e12.*faceVecs);
d23 = sum(e23.*faceVecs);
d31 = sum(e31.*faceVecs);
rows = [faceVIds(1,:), faceVIds(2,:), faceVIds(3,:)];
cols = [faceVIds(2,:), faceVIds(3,:), faceVIds(1,:)];
difVec = [d12,d23,d31];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the topology of the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [edges, adjFaces] = mesh_topology(Mesh)
%
numV = size(Mesh.vertexPoss, 2);
numF = size(Mesh.faceVIds, 2);
v1Ids = Mesh.faceVIds(1, :);
v2Ids = Mesh.faceVIds(2, :);
v3Ids = Mesh.faceVIds(3, :);
%
rows_A_vv = [v1Ids, v2Ids, v3Ids];
cols_A_vv = [v2Ids, v3Ids, v1Ids];
vals_A_vv = ones(1, 3*numF);
A_vv = sparse(rows_A_vv, cols_A_vv, vals_A_vv, numV, numV);
A_vv = A_vv + A_vv';
[rows, cols, vals] = find(A_vv);
ids = find(rows < cols);
edges = [rows(ids), cols(ids)]';
numE = size(edges, 2);
%
cols_A_vf = ones(3,1)*(1:numF);
rows_A_vf = Mesh.faceVIds;
vals_A_vf = ones(3, numF);
A_vf = sparse(rows_A_vf, cols_A_vf, vals_A_vf, numV, numF);
for vId = 1 : numV
    nFIds{vId} = find(A_vf(vId,:));
end
%
adjFaces = zeros(2, numE);
for eId = 1 : numE
    v1Id = edges(1, eId);
    v2Id = edges(2, eId);
    stats = [];
    for j = 1 : length(nFIds{v1Id})
        fId = nFIds{v1Id}(j);
        flag = is_face(Mesh.faceVIds(:, fId), [v1Id, v2Id]);
        if flag == 1
            stats = [stats, fId];
        end
    end
    adjFaces(1:length(stats), eId) = stats;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = is_face(vids, edge)
flag = 0;
if min(vids(1:2)) == edge(1) && max(vids(1:2)) == edge(2)
    flag = 1;
    return;
end
if min(vids([1,3])) == edge(1) && max(vids([1,3])) == edge(2)
    flag = 1;
    return;
end
if min(vids([2,3])) == edge(1) && max(vids([2,3])) == edge(2)
    flag = 1;
    return;
end

