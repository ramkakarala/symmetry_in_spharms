%% This script is provided for review purposes of the CVPR submission
% No 690, "3-D symmetry estimation in the phase domain"
%%
clear;
close all;
%% now read in the shape in off format
[V,DenseV,triobject]=readshapeinoffandplotmesh;
NptsV=max(size(V));
NptsDV=max(size(DenseV));
VcentD = mean(DenseV);
VmuO = V-ones(NptsV,1)*mean(V);
VmuD = DenseV-ones(NptsDV,1)*mean(DenseV);
%% pick a random rotation
Ro = euleranglestorotmatrix(2*pi*rand,pi*rand,2*pi*rand);
%% rotate original shape
VmuOR = VmuO*Ro;
VmuOR = VmuOR-ones(size(VmuOR,1),1)*mean(VmuOR);  % subtract off any mean
%% reorient by PA -- this is not necessary to the method, but it saves us 
% from running many different rotations of each shape for experiments
CO = cov(VmuOR);
[VC,EC]=eig(CO);
VmuORC = VmuOR*VC;
[phi,thetat,rVO] = cart2sph(VmuORC(:,1),VmuORC(:,2),VmuORC(:,3));
theta=pi/2-thetat;
%% Spherical harmonics of original shape vertices
L=16;
sigma=0.0005;
[Fsmoothls,FrotO] = SPHARMsmoothscalarcomplexp(rVO,L,sigma,theta,phi,'normal');
%% rotate densified shape
VmuDR = VmuD*Ro;
VmuDR = VmuDR-ones(size(VmuDR,1),1)*mean(VmuDR);  % subtract off any mean
%% reorient densified
CD = cov(VmuDR);
[VD,ED]=eig(CD);
VmuDRC = VmuDR*VD;
[phi,thetat,rVD] = cart2sph(VmuDRC(:,1),VmuDRC(:,2),VmuDRC(:,3));
theta=pi/2-thetat;
%% SPH of densified
L=16;
sigma=0.0005;
[Fsmoothls,FrotD] = SPHARMsmoothscalarcomplexp(rVD,L,sigma,theta,phi,'normal');
%% compute EGI mapping from triangles in shape
egiRC=compute_egi(VmuORC,triobject);
egiRCmu = egiRC-ones(max(size(egiRC)),1)*mean(egiRC);
%% compute polar coordinates
[phi,thetat,rVE] = cart2sph(egiRCmu(:,1),egiRCmu(:,2),egiRCmu(:,3));
theta=pi/2-thetat;
%% SPH of EGI
L=16;
sigma=0.0005;
[Fsmoothls,FrotE] = SPHARMsmoothscalarcomplexp(rVE,L,sigma,theta,phi,'normal');
%% covariance matrix estimate as baseline
tic;
[Rb3d,ratioco,bestalphac,bestbetac,bestgammac]=symshcovar(VmuORC,FrotO,L);
toc;
%% nested grid search two variable
tic;
[Rb2,ratio2o,bestalpha2,bestbeta2,bestgamma2]=optimizesymplanec(FrotO,L,30,5,bestalphac,bestbetac);
toc;
%% search using ISA, a form of APF
tic;
[Rboa,ratioapfo,bestalpha2,bestbeta2,bestgamma2]=...
    optimizelinearphasefit_apf_v2(FrotO,L,10,5);
toc;
%% how did we do?
fprintf(1,'Original Vertex results: N-GRID = %f,SH-ISA = %f SH-COVAR=%f\n',...
    ratio2o,ratioapfo,ratioco);
%% repeat with dense points
tic;
[Rb3d,ratiocd,bestalphac,bestbetac,bestgammac]=symshcovar(VmuDRC,FrotD,L);
toc;
%% two variable grid search
tic;
[Rb2,ratio2d,bestalpha2,bestbeta2,bestgamma2]=optimizesymplanec(FrotD,L,30,5,bestalphac,bestbetac);
toc;
%% search with ISA
tic;
[Rb2,ratioapfd,bestalpha2,bestbeta2,bestgamma2]=...
    optimizelinearphasefit_apf_v2(FrotD,L,10,5);
toc;
%% how did we do now?
fprintf(1,'Dense vertex results: N-GRID = %f, SH-ISA = %f SH-COVAR=%f\n',...
    ratio2d,ratioapfd,ratiocd);
%% repeat above with egi
%%
tic;
[Rbce,ratioce,bestalphac,bestbetac,bestgammac]=symshcovar(egiRCmu,FrotE,L);
toc;
%% two variable
tic;
[Rb2e,ratio2e,bestalpha2,bestbeta2,bestgamma2]=optimizesymplanec(FrotE,L,30,5,bestalphac,bestbetac);
toc;
%% search using APF
tic;
[Rbeapf,ratioapfe,bestalpha2,bestbeta2,bestgamma2]=...
    optimizelinearphasefit_apf_v2(FrotE,L,10,5);
toc;

%% how did EGI do?
fprintf(1,'EGI results: N-GRID = %f,SH-ISA = %f SH-COVAR=%f\n',...
    ratio2e,ratioapfe,ratioce);
%% plot results using densified vertices with ISA result
figure
trisurf(triobject, VmuORC(:,1), VmuORC(:,2), VmuORC(:,3));
colormap(copper);
axis equal
axis off;
daspect([1 1 1]);
view(3);
if (ratioapfd < ratioapfo)
title(sprintf('E=%4.2f\n',ratioapfd));
hold on;
Pyr=[0 0 0 [0 0 1]*Rb2 [1 0 0]*Rb2];
drawPlane3d(Pyr,'g');  % from geom3d toolbox
hold off;
else
    title(sprintf('E=%4.2f\n',ratioapfo));
hold on;
Pyr=[0 0 0 [0 0 1]*Rboa [1 0 0]*Rboa];
drawPlane3d(Pyr,'g');  % from geom3d toolbox
hold off;
end;
%% plot results for egi
figure
trisurf(triobject, VmuORC(:,1), VmuORC(:,2), VmuORC(:,3));
colormap(copper);
axis equal
axis off;
daspect([1 1 1]);
view(3);
title('Plane from EGI shown on mesh');
hold on;
%cd geom3d
Pyr=[0 0 0 [0 0 1]*Rbeapf [1 0 0]*Rbeapf];
drawPlane3d(Pyr,'g');
hold off;
 