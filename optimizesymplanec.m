function [Rb,pctresidual,bestalpha,bestbeta,bestgamma]=...
    optimizesymplanec(F,L,N,Rpt,prevbestalpha,prevbestbeta)
% usage
% [Rb,pctresidual,bestalpha,bestbeta,bestgamma]=
%  optimizesymplanec(F,L,N,Rpt,prevbestalpha,prevbestbeta)
%
% Finds the Euler angles of the best linear phase fit to spherical
% harmonic coefficients F using nested grid search (N-GRID).  Parameters
% F is the original data
% L maximum order of the coeffs. F is Lx(2L+1)
% N number of grid elements in each angle
% Rpt number of repeats of the grid search
% prevbestalpa,beta can be set to outptus of SH-COVAR
% outputs are best rotation matrix and best alpha,beta,gamma in the Z-Y-Z eulera angles that
% provide the best linear phase fit, and the model fit residual (closer to 0
% means better)
%

if (nargin<5)
    prevbestalpha = 0;
    prevbestbeta = 0; 
end;
load littledfor16;
if (L>16)
    error('Cannot use precomputed wigner D');
end;
% calculate upper bound on fit 
sumnorm=0;
for k=2:L
    Fl=F(k,1:2*k-1);
    Funwrap{k}=[fliplr(Fl(k+1:end)),Fl(1),Fl(2:k)];
    Funwrapt{k}=conj(Funwrap{k}');
    sumnorm = sumnorm+Fl*Fl';
end;

% grid
alpha = (0:2:2*N-1)/N*(2*pi);
beta = (0:2:2*N-1)'/N*pi/2;
Npts = length(alpha)*length(beta)+2;
as=zeros(Npts,1);
bs=as;
as(1)=prevbestalpha;
bs(1)=prevbestbeta;
bs(end)=pi;
as(end)=0;
ind = 1;
for c=2:N-1
    for r=1:N
        ind=ind+1;
        as(ind)=alpha(r); bs(ind)=beta(c);
    end;
end;


% coarse grid search
bestfitvalue = 0;
bestalpha = prevbestalpha;
bestbeta = prevbestbeta;
bestgamma = 0;
Ltest = L;
Nangles= N;
Nrepeats= Rpt;
for n=1:Npts
    fitvalue = 0;
    Dl = MakeEquivalentWignerMatrices(Ltest,DC,DCt, as(n),bs(n),0);
    for k=2:Ltest
         Dt = Dl{k}'*conj(Dl{k});
        lvalue = real(Funwrap{k}*Dt*Funwrapt{k});
        fitvalue = fitvalue+lvalue;
    end;
    if fitvalue > bestfitvalue
        bestfitvalue = fitvalue;
        bestalpha = as(n);
        bestgamma = 0;
        bestbeta = bs(n);
    end;
end;
fprintf(1,'coarse search done using %d points\n',Npts);
pctresidual= 100*(sumnorm-bestfitvalue)/(sumnorm);
fprintf(1,'Pct residual after coarse search = %f\n',pctresidual);

% repeat using local grid search
steps=10;
rangea = pi;
rangeb = pi/2;


for t=1:Nrepeats
    oldbestfit = bestfitvalue;
    rangea = rangea/steps;
    rangeb = rangeb/steps;
    asv=linspace(bestalpha-rangea,bestalpha+rangea,steps);
    bsv=linspace(bestbeta-rangeb,bestbeta+rangeb,steps);
    [asm,bsm]=meshgrid(asv,bsv);
    as=asm(:);
    bs=bsm(:);
    Npts2 = length(as);
    for n=1:Npts2
        fitvalue = 0;
        Dl = MakeEquivalentWignerMatrices(Ltest,DC,DCt, as(n),bs(n),0);
        for k=2:Ltest
            Dt = Dl{k}'*conj(Dl{k});
            lvalue = real(Funwrap{k}*Dt*Funwrapt{k});
            fitvalue = fitvalue+lvalue;
        end;
        if fitvalue > bestfitvalue
            bestfitvalue = fitvalue;
            bestalpha = as(n);
            bestgamma = 0;
            bestbeta = bs(n);
        end;
    end;
    pctresidual= 100*(sumnorm-bestfitvalue)/(sumnorm);
    fprintf(1,'Pct residual after repeat %d = %f\n',t,pctresidual);
    fractchange=(bestfitvalue-oldbestfit)/oldbestfit;
    if (fractchange<0.01)
        break;
    end;

end;

Rb = euleranglestorotmatrix(bestalpha,bestbeta,0);