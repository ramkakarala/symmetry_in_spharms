function [vertices,densePoints,tri] = readshapeinoff(fname,NDENSE)
% usage
%        pts = readshapeinoff(fname)
% reads the 3d coordinates of a shape in OFF format
if (nargin < 1)
    % rad a shape in .off file
    [fname,pname]=uigetfile('*.off','Get shape in .off file');
    fname = [pname fname];
end;
if (nargin < 2)
    NDENSE = 100000;
end;

fid = fopen(fname,'rt');
%%
line=fgetl(fid);
line=fgetl(fid);
par = sscanf(line,'%d %d %d');
Npts = par(1);
numFaces = par(2);
%% read points
vertices = zeros(Npts,3);
pts = 1;
while (~feof(fid) && pts <= Npts)
    line = fgetl(fid);
    par = sscanf(line,'%f %f %f');
    if ~isempty(par)
        vertices(pts,:)=par(1:3);
        % disp(par);
        pts=pts+1;
    end;
    if (rem(pts,1000)==0)
        disp(pts);
    end;
end;
%% Triangulate the faces

tri = [];
% Go over each face
for i=1:numFaces
    line = fgetl(fid);
    var = sscanf(line,'%f %f %f');
    faces(i,:) = var;
    if(faces(i,1) == 3)
        tri = [tri;faces(i,2:end)];
    else
        %Triangulate individual polygons and store triangles
    end;
    
end;
tri = tri+1;
figure,
subplot(1,2,1);
trisurf(tri, vertices(:,1), vertices(:,2), vertices(:,3));
axis equal
axis off;
daspect([1 1 1]);
view(3);
title('Mesh');

fclose(fid);


%% Dense Sampling of points

fprintf('Calculating the dense points ...\n');

areas = zeros(length(tri), 1);

zaxis = [0,0,1];
%Find area of each triangle
for ii = 1:length(tri)
    p0 = vertices(tri(ii,1),:);
    p1 = vertices(tri(ii,2),:);
    p2 = vertices(tri(ii,3),:);
    T = [p0;p1;p2];
    %Find normal to the triangle
    normal = cross(p1-p0, p2-p0);
    
    %Find rotation matrix to rotate normal to zaxis i.e., make triangle
    %parallel to the X-Y plane
    R = RotMatV1V2(normal, zaxis);
    
    %Rotate Triangle
    RotatedT = T*R;
    
    X = RotatedT(:,1);
    Y = RotatedT(:,2);
    areas(ii) = polyarea(X, Y);
end;

totalArea = sum(areas);
normArea = areas./totalArea;

%Sample points proportional to the area
P = [];
densePoints = []; %Store the dense sampled points
zaxis = [0,0,1];
try
    for ii = 1:length(tri)
        p0 = vertices(tri(ii,1),:);
        p1 = vertices(tri(ii,2),:);
        p2 = vertices(tri(ii,3),:);
        T = [p0;p1;p2];
        
        %Find normal to the triangle
        normal = cross(p1-p0, p2-p0);
        
        %Find rotation matrix to rotate normal to zaxis i.e., make triangle
        %parallel to the X-Y plane
        R = RotMatV1V2(normal, zaxis);
        RotatedT = T*R;
        
        %Number of dense points in tri i is proportional to its area.
        N = ceil(NDENSE*normArea(ii));
        P1 = RotatedT(1,1:2);
        P2 = RotatedT(2,1:2);
        P3 = RotatedT(3,1:2);
        
        meanZ = mean(RotatedT(:,3));
        
        t = sqrt(rand(N,1));
        s = rand(N,1);
        P = (1-t)*P1 + (t.*(1-s))*P2+(t.*s)*P3;
        Z = ones(size(P,1), 1).*meanZ;
        newP = [P, Z];
        
        %Rotate the sampled points back to the original plane
        densePoints = [densePoints; newP*R'];
    end;
catch
    ii
end;
subplot(1,2,2); plot3(densePoints(:,1),densePoints(:,2),densePoints(:,3),'kx'),axis equal off;
title('point cloud');

end
%
%----------------------------------------------
%
function rotationMatrix = RotMatV1V2(x, y)

%Calculates the rotation matrix that rotates a 3D vector x to a 3D vector y

rotationAxis=cross(x,y);
if(sum(rotationAxis) == 0)
    rotationMatrix = eye(3);
else
    rotationAngle=acosd( dot(x,y)/norm(x)/norm(y) );
    
    rotationMatrix=R3d(rotationAngle, rotationAxis);
end;
end

function R=R3d(deg,u)
%R3D - 3D Rotation matrix counter-clockwise about an axis.
%
%R=R3d(deg,u)
%
% deg: The counter-clockwise rotation about the axis u in degrees.
%


R=eye(3);
u=u(:)/norm(u);
x=deg; %abbreviation

for ii=1:3
    
    v=R(:,ii);
    
    R(:,ii)=v*cosd(x) + cross(u,v)*sind(x) + (u.'*v)*(1-cosd(x))*u;
    %Rodrigues' formula
    
end

end
