function egi = compute_egi(vertex,face)

% compute_egi - compute the extended gaussian image of a mesh model
%
%   egi = compute_egi(vertex,face);
%
%   egi is a Nfaces x 3 matrix of vectors in 3-D.  Each vector is the 
% normal to a face in the mesh model, and its length is proportional to
% area of the face.
%
%   based on code (c) 2004 Gabriel Peyré
% 


[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,2);
nvert = size(vertex,2);
normal = zeros(3, nvert);

% unit normals to the faces
normalf = crossp( vertex(:,face(2,:))-vertex(:,face(1,:)), ...
                  vertex(:,face(3,:))-vertex(:,face(1,:)) );
orignormalf = normalf;
d = sqrt( sum(normalf.^2,1) ); d(d<eps)=1;  
%note this d is twice the area of the triangular face
normalf = normalf ./ repmat( d, 3,1 );

%unit normal to the vertex
 normal = zeros(3,nvert);
 for i=1:nface
     f = face(:,i);
    for j=1:3
        normal(:,f(j)) = normal(:,f(j)) + normalf(:,i);
    end
end
% normalize
d = sqrt( sum(normal.^2,1) ); d(d<eps)=1;
normal = normal ./ repmat( d, 3,1 );

% enforce that the normal are outward
v = vertex - repmat(mean(vertex,1), 3,1);
s = sum( v.*normal, 2 );
if sum(s>0)<sum(s<0)
    % flip
    normal = -normal;
    normalf = -normalf;
    orignormalf = -orignormalf;
end

egi = 0.5*orignormalf';  % the 0.5 corrects for the area of parallelogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(1,:) = x(2,:).*y(3,:) - x(3,:).*y(2,:);
z(2,:) = x(3,:).*y(1,:) - x(1,:).*y(3,:);
z(3,:) = x(1,:).*y(2,:) - x(2,:).*y(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vertex,face] = check_face_vertex(vertex,face, options)

% check_face_vertex - check that vertices and faces have the correct size
%
%   [vertex,face] = check_face_vertex(vertex,face);
%
%   Copyright (c) 2007 Gabriel Peyre

vertex = check_size(vertex,2,4);
face = check_size(face,3,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = check_size(a,vmin,vmax)
if isempty(a)
    return;
end
if size(a,1)>size(a,2)
    a = a';
end
if size(a,1)<3 && size(a,2)==3
    a = a';
end
if size(a,1)<=3 && size(a,2)>=3 && sum(abs(a(:,3)))==0
    % for flat triangles
    a = a';
end
if size(a,1)<vmin ||  size(a,1)>vmax
    error('face or vertex is not of correct size');
end

function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) 2003 Gabriel Peyré


fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);



[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;
% read Face 1  1088 480 1022
[A,cnt] = fscanf(fid,'%d %d %d %d\n', 4*nface);
if cnt~=4*nface
    warning('Problem in reading faces.');
end
A = reshape(A, 4, cnt/4);
face = A(2:4,:)+1;


fclose(fid);

