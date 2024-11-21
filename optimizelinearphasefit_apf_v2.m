function [Rb,pctresidual,bestalpha,bestbeta,bestgamma]=...
    optimizelinearphasefit_apf(F, L, N, M, b, alpha)
%
% optimize symmetry plane fit using Interacting Simulated Annealing
% a form of annealed particle filter
%
    load littledfor16;  
    sumnorm=0;
    for k=2:L
          Fl=F(k,1:2*k-1);
          sumnorm = sumnorm+Fl*Fl';
    end;
    fun = @(x)optfun(x, L, F, DC, DCt);
    if nargin < 5
        N = 100;
        M = 8;
        b = 0.25;
        alpha = 5;
    end
    [vec1, val1] = loptimize_apf(fun, N, M, b, alpha);
    [vec2, val2] = loptimize_apf(fun, N, M, b, alpha);
    if val1 < val2
        vec = vec2;
        val = val2;
    else
        vec = vec1;
        val = val1;
    end
    [bestalpha, bestbeta, bestgamma] = qtoea(vec);
    Rb = euleranglestorotmatrix(bestalpha,bestbeta,bestgamma);
    pctresidual= 100*(sumnorm-val)/(sumnorm);
end

function fitvalue = optfun(unitr3vec, Ltest, F, DC, DCt)
    assert(abs(norm(unitr3vec)-1) < 1e-12);
	fitvalue = 0;
    [alpha, beta, gamma] = qtoea(unitr3vec);
    %Dl = makewignermatrices(Ltest,DC,alpha,gamma);
    Dl = MakeEquivalentWignerMatrices(Ltest,DC,DCt,alpha,beta,gamma);
    for k=2:Ltest
        Fl = F(k,1:2*k-1);
        funwrap = [fliplr(Fl(k+1:end)),Fl(1),Fl(2:k)];
        Dt = Dl{k}'*conj(Dl{k});
        lvalue = real(funwrap*Dt*conj(funwrap'));
        fitvalue = fitvalue+lvalue;
    end;
end
function [psi, theta, phi] = qtoea(Q)
    psi=atan2((Q(:,2).*Q(:,3)-Q(:,1).*Q(:,4)),(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)));
    theta=acos(Q(:,4).^2-Q(:,1).^2-Q(:,2).^2+Q(:,3).^2);
    phi=atan2((Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)));
end

function [vec, val] = loptimize_apf(fun, N, M, alpha, gamma)
    start = [0, 1, 0, 0];
    samples = repmat(start, N, 1);
    likely = zeros(N, 1);
    historyl = zeros(N, M);
    history = cell(1, M);
    res = 2;
    draw = 0;
    kappa = 0;
    if draw ~= 0; h = figure; hold on; end
    for i=1:M
        for j = 1:N
            samples(j, :) = randvonMisesFisherm(4, 1, kappa, samples(j, :));
            likely(j) = fun(samples(j, :));
        end
        if res == 2; history{i} = samples; historyl(1:N, i) = likely; end
        if draw ~= 0; cla; plot3(samples(:, 1), samples(:, 2), samples(:, 3), '.r');xlim([-1, 1]);ylim([-1, 1]);zlim([-1, 1]); end

        beta = getbeta(likely, 1, alpha, 0.01);
        anpdf = getpdf(beta*likely);
        selind = isanresample(anpdf, N);
        samples = samples(selind, :);
        if draw ~= 0; plot3(samples(:, 1), samples(:, 2), samples(:, 3), '.b'); pause(1); end
        if i == 1
            kappa = 1;
        end
        kappa = kappa*gamma;
    end
    if draw ~= 0; close(h); end
    
    if res == 1
        vec = sum(samples, 1);
        Rbar = norm(vec);
        vec = vec ./ Rbar;    
    elseif res == 2
        [j, i] = find(historyl == max(max(historyl)));
        vec = history{i}(j, :);
    end
    val = fun(vec);
end

function pdf = getpdf(log_unnormalized_weights)
    pdf = exp(log_unnormalized_weights - max(log_unnormalized_weights));
    if isnumeric(sum(pdf)) && sum(pdf) ~= 0
        pdf = pdf ./ sum(pdf);
    else
        pdf = zeros(size(log_unnormalized_weights));
        [~, ind] = max(log_unnormalized_weights);
        pdf(ind) = 1;
    end
end
function ind = isanresample(anpdf, N)
    mind = (anpdf >= rand(N, 1)*max(anpdf));
    ind = [find(mind); mresample(anpdf, N - sum(mind))];
    assert(length(ind) == N);
end
function [selind] = mresample(pdf, N)
    while(1)
        cdf = cumsum([0; pdf(:)]);
        samples = rand(1, N);
        selind = sum(repmat(cdf, 1, N) <= repmat(samples, size(cdf, 1), 1), 1)';
        if all(selind <= numel(pdf))
            break;
        end
        fprintf('mresamples looping due to numerical instability\n');
    end
end
function beta = getbeta(lnlikelihood, beta, alphadesired, ~)
    global ln_sample_likelihood;
    global adesired;
    assert(alphadesired >= 1/length(lnlikelihood), 'Cooling too quick');
    ln_sample_likelihood = lnlikelihood;
    adesired = alphadesired ;
    beta = fzero(@getalpha, beta);
end

function alpha = getalpha(beta)
    global ln_sample_likelihood;
    global adesired;
    if isinf(beta)
        alpha = 1/length(ln_sample_likelihood);
    elseif beta <= 0
        alpha = 1;
    else
        func = beta * ln_sample_likelihood;
        func = func - max(func);
        func = exp(func);
        if isfinite(sum(func)) && sum(func) ~= 0
            func = func ./ sum(func);
            invD = sum(func .^ 2);
            assert(invD > 0);
            alpha = 1 / (invD * length(func));
        else
            alpha = 1/length(ln_sample_likelihood);
        end
    end
    alpha = alpha - adesired;
end
function [ X ] = randvonMisesFisherm(m, n, kappa, mu)
% RANDVONMISESFISHERM Random number generation from von Mises Fisher
% distribution.
% X = randvonMisesFisherm(m, n, kappa) returns n samples of random unit 
% directions in m dimensional space, with concentration parameter kappa,
% and the direction parameter mu = e_m
% X = randvonMisesFisherm(m, n, kappa, mu) with direction parameter mu
% (m-dimensional column unit vector)
%
% Sungkyu Jung, Feb 3, 2010.

if nargin < 3, help randvonMisesFisher3, return, end
if nargin == 3,  muflag = false;
else muflag = true;
end

if m < 2; 
    disp('Message from randvonMisesFisherm.m: dimension m must be > 2'); 
    disp('Message from randvonMisesFisherm.m: Set m to be 2'); 
    m = 2;
end

if kappa < 0; 
    disp('Message from randvonMisesFisherm.m: kappa must be >= 0'); 
    disp('Message from randvonMisesFisherm.m: Set kappa to be 0'); 
    kappa = 0;
end

%
% the following algorithm is following the modified Ulrich's algorithm 
% discussed by Andrew T.A. Wood in "SIMULATION OF THE VON MISES FISHER 
% DISTRIBUTION", COMMUN. STATIST 23(1), 1994.

% step 0 : initialize
b = (-2*kappa + sqrt(4*kappa^2 + (m-1)^2))/(m-1);
x0 = (1-b)/(1+b);
c = kappa*x0 + (m-1)*log(1-x0^2);

% step 1 & step 2
nnow = n; w = [];
%cnt = 0;
while(true)
    ntrial = max(round(nnow*1.2),nnow+10) ;
    Z = betarnd((m-1)/2,(m-1)/2,ntrial,1);
    U = rand(ntrial,1);
    W = (1-(1+b)*Z)./(1-(1-b)*Z);
    
    indicator = kappa*W + (m-1)*log(1-x0*W) - c >= log(U);
    if sum(indicator) >= nnow
        w1 = W(indicator);
        w = [w ;w1(1:nnow)];
        break;
    else
        w = [w ; W(indicator)];
        nnow = nnow-sum(indicator);
        %cnt = cnt+1;disp(['retrial' num2str(cnt) '.' num2str(sum(indicator))]);
    end
end

% step 3
V = UNIFORMdirections(m-1,n);
X = [repmat(sqrt(1-w'.^2),m-1,1).*V ;w'];

if muflag
    mu = mu / norm(mu);
    X = rotMat(mu)'*X;
end
end


function V = UNIFORMdirections(m,n)
% generate n uniformly distributed m dim'l random directions
% Using the logic: "directions of Normal distribution are uniform on sphere"

V = zeros(m,n);
nr = randn(m,n); %Normal random 
for i=1:n
    while 1
        ni=nr(:,i)'*nr(:,i); % length of ith vector
        % exclude too small values to avoid numerical discretization
        if ni<1e-10 
            % so repeat random generation
             nr(:,i)=randn(m,1);
        else
             V(:,i)=nr(:,i)/sqrt(ni);
            break;
        end
    end
end

end

function rot = rotMat(b,a,alpha)
% ROTMAT returns a rotation matrix that rotates unit vector b to a
%
%   rot = rotMat(b) returns a d x d rotation matrix that rotate
%   unit vector b to the north pole (0,0,...,0,1)
%
%   rot = rotMat(b,a ) returns a d x d rotation matrix that rotate
%   unit vector b to a
%
%   rot = rotMat(b,a,alpha) returns a d x d rotation matrix that rotate
%   unit vector b towards a by alpha (in radian)
%
%    See also .

% Last updated Nov 7, 2009
% Sungkyu Jung


[s1 s2]=size(b);
d = max(s1,s2);
b= b/norm(b);
if min(s1,s2) ~= 1 || nargin==0 , help rotMat, return, end  

if s1<=s2;    b = b'; end

if nargin == 1;
    a = [zeros(d-1,1); 1];
    alpha = acos(a'*b);
end

if nargin == 2;
    alpha = acos(a'*b);
end
if abs(a'*b - 1) < 1e-15; rot = eye(d); return, end
if abs(a'*b + 1) < 1e-15; rot = -eye(d); return, end

c = b - a * (a'*b); c = c / norm(c);
A = a*c' - c*a' ;

rot = eye(d) + sin(alpha)*A + (cos(alpha) - 1)*(a*a' +c*c');
end
