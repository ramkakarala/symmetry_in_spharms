%% This script recreates figure 1
clear;
close all;
%% now read in the world map
F = readS2Kitfdata('mapofworld.dat');
BW = max(size(F))/2;
%% compute transform
[x,y,z]=sphere(2*BW-1);
Npts = 2*BW;
jkvec = 0:2*BW-1;
thetavec = pi/(4*BW)*(2*jkvec+1);
varphivec = 2*pi*jkvec/(2*BW);
[varphi,theta]=meshgrid([varphivec 2*pi],[0 thetavec pi]);
% note the addition of 0, pi, and 2pi above makes the spherical
% triangulation complete. Those points are NOT sampled by S2kit
% for purposes of computing Spherical harmonic transform
tri=delaunay(theta,varphi);
x=cos(varphi).*sin(theta);
y=sin(varphi).*sin(theta);
z=cos(theta);
%% extrapolate
Fe = zeros(size(theta));
Fe(2:2*BW+1,1:2*BW) = F;
% the last column of Fe is at varphi = 2pi, same as 0
Fe(:,end) = Fe(:,1);
% first row of Fe is where theta = 0; extrapolate it by adding deriv
linearexp = F(1,:)+F(1,:)-F(2,:);
Fe(1,:) = mean(linearexp)*ones(1,size(theta,2));
% last row of Fe is where theta = pi; extrapolate it
linearexp = F(end,:)+F(end,:)-F(end-1,:);
Fe(end,:) = mean(linearexp)*ones(1,size(theta,2));
%%
mF = min(Fe(:));
MF = max(Fe(:));
Fs = Fe - min(Fe(:))+1.0;  % now a function value of zero means it lies on the sphere
%Fs = 1 + (Fe-mF)/(MF-mF);
figure
trisurf(tri,x.*Fs,y.*Fs,z.*Fs,Fs);  % plot F(P)*P in Foley's notation
title('World');
axis equal
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%% Fourier expansion movie
figure;

L=BW-1;
clear M;
ll = L;

sigma=0.0005; %important value, chosen after some experiments
[Fsmoothworld,fourier_coeffsworld]= SPHARMsmoothscalarcomplexp(Fe,ll,sigma,theta,varphi,'normal');
if (norm(imag(Fsmoothworld))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
Fsmoothr=real(Fsmoothworld);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%% replace with their magnitudes
L = BW-1;
Fmag = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    fmag = [sqrt(fvec*fvec'),zeros(1,flen-1)];
    % coeffs are [0,+1,+2,..+ell,-1,-2...-ell]
    % must unwrap them into [-ell,..,-1,0,1,...ell];
   % funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
   % frotvec = funwrap*Dmat{k+1};
    % now rewrap them to store;
   % frotwrap = [frotvec(k+1),frotvec(k+2:end),fliplr(frotvec(1:k))];
    Fmag(k+1,1:flen) = fmag; %rotwrap;
end;
%% reconstruct the rotated function from coefficients only
fmagrecon = SPHARMreconscalarcomplex(Fmag,L,sigma,theta,varphi,'normal');
if (norm(imag(fmagrecon))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
Fsmoothr=real(fmagrecon);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('rotated original, nagnitude only recon L=%d',ll));
%%-----------------------------------------------------------
%% reconstruct the function from real part of coefficients only
frealrecon = SPHARMreconscalarcomplex(real(fourier_coeffsworld),L,sigma,theta,varphi,'normal');
if (norm(imag(frealrecon))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
Fsmoothr=real(frealrecon);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('rotated original, nagnitude only recon L=%d',ll));
%%-----------------------------------------------------------
%% replace with the phase
L = BW-1;
Fpha = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    fpha = (1/flen)*fvec/norm(fvec);  
    % note 1/flen is a non-informative magnitude that comes from the delta
    % function
    Fpha(k+1,1:flen) = fpha; 
end;
%% reconstruct the rotated function from phase only
fpharecon = SPHARMreconscalarcomplex(Fpha,L,sigma,theta,varphi,'normal');
if (norm(imag(fpharecon))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
Fsmoothr=real(fpharecon);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF); 
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('phase only recon L=%d',ll));
%%------------------------------------------------
%%
I=imread('cameraman.tif');
figure,imshow(I);
IF2=fft2(I);
Ip2=IF2./abs(IF2);
I2=zeros(size(I));
I2(1:4,1:4)=1;
I2F=fft2(I2);
Ip2randmag=Ip2.*abs(I2F);
ip2=real(ifft2(Ip2randmag));
ip2m=ip2/max(ip2(:))*255;
figure,imshow(uint8(ip2m)),