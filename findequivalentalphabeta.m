function [alpha,beta,Rbe]=findequivalentalphabeta(ntest)
% given a normal vector, find the two angles that are needed to map
% the Y axis to it.  
beta = atan2(ntest(3),-ntest(1));
if beta < 0
    beta = beta+pi;
end;
alpha = acos(ntest(2));
if (abs(ntest(1))>abs(ntest(3)))
if (sign(ntest(1)) ~= sign(-cos(beta)*sin(alpha)))
    alpha = -alpha;
end;
else
    if (sign(ntest(3)) ~= sign(sin(beta)*sin(alpha)))
    alpha = -alpha;
end;
end;
Rbe=euleranglestorotmatrix(alpha,beta,0);