Mesh.centriods = zeros(3, size(Mesh.faceVIds, 2));

for i = 1:size(Mesh.faceVIds, 2)

Mesh.centriods(:, i) =  Mesh.vertexPoss(:, Mesh.faceVIds(1, i))
                        + Mesh.vertexPoss(:, Mesh.faceVIds(2, i))
                        + Mesh.vertexPoss(:, Mesh.faceVIds(3, i));
end

Mesh.centriods = Mesh.centriods * 1./3.;

cent_temp = Mesh.centriods';

scatter3(cent_temp(:,1),cent_temp(:,2),cent_temp(:,3),'r')
