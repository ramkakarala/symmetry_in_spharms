function Dl = MakeEquivalentWignerMatrices(L, DC, DCt, alpha, beta, gamma)
% use dcomposition to speed up computation of wigner matrices
a1 = alpha-pi/2;
a2 = gamma+pi/2;
Dl = cell(L,1);
for k=0:L-1
    ind = -k:k;
    D = (exp(-1i*ind'*a1)*exp(-1i*ind*beta) .* DCt{k+1})*(ones(length(ind),1)*exp(-1i*ind*a2) .* DC{k+1});
    t = real(trace(D));
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else
            ch = 1;
        end;
    end;
    if (abs (t-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    Dl{k+1}=D;
end;