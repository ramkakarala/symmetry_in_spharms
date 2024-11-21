function [Rb,ratio,bestalpha,bestbeta,bestgamma]=symshcovar(X,F,L)
% usage
%       [Rb,ratio,bestalpha,bestbeta,bestgamma]=symshcovar(X,F,L);
% X is the Nx3 data, F is its spherical harmonics, and L is the order
% estimates symmetry using the covariance matrix to provide candidates
load littledfor16;
if (L>16)
    error('Cannot use precomputed wigner D');
end;
CovD=cov(X);
[Vc,Dc]=eig(CovD);
Vj=Vc';
for i=1:3
    [fit(i),alphav(i),betav(i),resid(i)]=evalfitnessofvect(F,L,Vj(i,:),DC,DCt);
end;
[bestvfit,ind]=max(fit);
bestalpha=alphav(ind); bestbeta=betav(ind); bestfitv=fit(ind);
ratio=resid(ind);
Rb=euleranglestorotmatrix(bestalpha,bestbeta,0);
bestgamma = 0;
fprintf(1,'Covariance residual=%f\n',ratio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fitvalue,alpha,beta,pctresidual]=evalfitnessofvect(F,L,vect,DC,DCt)
% usage
%        fitval=evalfitnessofvect(F,L,vect)
% F is coefficients; L is order; vector normal vector of symmetry plane
%
[alpha,beta]=findequivalentalphabeta(vect);
D1=MakeEquivalentWignerMatrices(L,DC,DCt, alpha,beta,0);
fitvalue = 0;
for k=2:L
    Fl = F(k,1:2*k-1);
    funwrap = [fliplr(Fl(k+1:end)),Fl(1),Fl(2:k)];
    Dt = D1{k}'*conj(D1{k});
    lvalue = real(funwrap*Dt*conj(funwrap'));
    fitvalue = fitvalue+lvalue;
end;
sumnorm=0;
for k=2:L
      Fl=F(k,1:2*k-1);
      sumnorm = sumnorm+Fl*Fl';
end;
pctresidual= 100*(sumnorm-fitvalue)/(sumnorm);
