function R = euleranglestorotmatrix(alpha,beta,gamma)
% usage
%        R = euleranglestorotmatrix(alpha,beta,gamma)
%
% Provides the 3x3 rotation matrix R that provides the rotation
% defined by the Euler angles. Note: it is assumed that R is applied
% to row vectors x, i.e., xR = y.  The euler angles are 
% alpha = rotation around z axis
% beta = rotation around y axis following alpha (note sign)
% sign matches wigner matrix Dl{k+1} = exp(-j*alpha)*DC{k+1}*exp(-j*gamma)
% where DC{k+1} = dmatrixbeta(L+1,betao);
% finally
% gamma = rotation around z axis following alpha, beta

Ra = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
Rb = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
Rg = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
R  = Ra*Rb*Rg;

