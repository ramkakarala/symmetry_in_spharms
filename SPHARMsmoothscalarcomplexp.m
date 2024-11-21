function [frecon,coeff] = SPHARMsmoothscalarcomplexp(f,L,sigma,theta,varphi,opmode)
%---------------------------------------------------------------------------------------
%[frecon, coeff]=SPHARMsmoothscalarcomplexp(f,L,sigma,theta,varphi)
%
%
%    f            : 2D matrix of function values one at each theta, varphi
%
%    L            : The maximal degree of weighted-SPHARM representation.
%                       Read the paper below to find it optimally.
%
%  sigma          : bandwith of weighted-SPHARM representation
%                       It is the bandwidth of heat kernel on unit sphere.
%                       When sigma=0, it is the traditional SPHARM representation. 
%                       range beween 0.0001 and 0.01 will be sufficient for cortical surfaces.
%
% f_smooth: The weighted-SPHARM result. The dimension is identical to coord.
%
% coeff   : The estimated SPHARM coefficients (Fourier coefficients) given as a structured array
%                        containg coeff.x, coeff.y, coeff.z
%                        coeff.x is the SPHARM coefficients of x-cooridinate given as (L+1)*(2L+1)*n_subject 
%                        coeff.x(3,:,10) is the 2nd degree SPHARM coefficients of all order for the 10th subject.  
%
% opmode: one of the following 
%         'magnitude only'
%         'phase only'
%          none
%
% (C) Moo K. Chung, 2006, 2007
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mchung@stat.wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to understand the notations and the algorithm.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
%
% Update history Sept 19 2006; July 5, 2007
%----------------------------------------------------------------------------------------

%reshape for use with wspharm conventions
magnitudeonly = 0;
phaseonly = 0;
normalizedc = 0;
if nargin > 5
    if (~isempty(strfind(lower(opmode),'magnitude')))
        magnitudeonly = 1;
    end;
    if (~isempty(strfind(lower(opmode),'phase')))
        phaseonly = 1;
    end;
    if (~isempty(strfind(lower(opmode),'normal')))
       normalizedc = 1;
    end;  
end;
n_vertex = prod(size(f));   % the number of vertices in a mesh. 
[R,C]=size(f);
f = reshape(f,1,n_vertex);
latitude = reshape(theta,1,n_vertex);
longitude = reshape(varphi,1,n_vertex);

cestimate=zeros(n_vertex,1); 
beta=zeros(L+1,2*L+1,1);


% 0-degree SPHARM coefficients. Step 2 in the iterative resiual fitting (IRF) algorithm. See the paper.

Y=Y_l(0,latitude,longitude)'; % Y is now a col vector
Ycommon=pinv(Y'*Y)*Y';
    
betal=(Ycommon*f');
if (magnitudeonly)
    betal = abs(betal);
end;
fsmooth = Y*betal;
cestimate = fsmooth;
beta(1,1) = betal';
nsum = 0;
for l=1:L
  fprintf(1,'%d,',l); 
  if (rem(l,20)==0)
      fprintf(1,'\n');
  end;
    
       
        % Step 4: residual. See the paper for detail 
        f_j = f'-cestimate; %(:,i);
        
        temp = Y_l(l,latitude,longitude);
        Y=temp';
        % real(Y_l) gives negative harmonics imag(Y_l) gives positive harmonics 
        clear temp;
        pfact = diag((-1).^(1:l));
        Yc = conj(Y(:,2:(l+1)))*pfact;
        Ya = [Y Yc];
        %Ya =[real(Y)   imag(Y(:,2:(l+1)))];
        Ycommon=pinv(Ya'*Ya)*Ya';
        
        % Step 5: refitting the residual. See the paper for detail
        betal = Ycommon*f_j;
       
        if magnitudeonly
            betalmag = zeros(2*l+1,1);
            betalmag(1)=norm(betal);  % first element has the magnitude
            betal = betalmag;
        elseif phaseonly
            betal = betal/norm(betal);
        end;
        nsum = nsum+norm(betal);
        beta(l+1,1:2*l+1)=betal';
        fsmooth=Ya*betal;
        cestimate=cestimate+ exp(-l*(l+1)*sigma)*fsmooth;
        
  
end;

fprintf(1,'done\n');
if normalizedc
    beta(1,1)=beta(1,1)+nsum;
    Y = Y_l(0,latitude,longitude)'; % Y is now a col vector
    fsmooth = Y*nsum;
    cestimate = cestimate+fsmooth;
end;
coeff=beta;
frecon=reshape(cestimate,R,C);
