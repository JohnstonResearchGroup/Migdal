function out = get_mu(Z,X,ek,WN,nk,beta,fill)
%GET_MU Find the chemical potential using the bisection method.
%   The bisection method is guaranteed to converge if DOS(w) >= 0.
erreps = 1e-9;

%set the initial values of the chemical potential.
muL = min(min(ek));
muR = max(max(ek));

fillL = get_filling(Z,X,ek,WN,nk,beta,muL);
fillR = get_filling(Z,X,ek,WN,nk,beta,muR);

% In some rare cases, the self-energy X is so large that filling is not in 
% interval [fillL fillR]. In those cases, we expand the interval [muL muR] 
% using the while loop below.
iter = 0;
wbnd = muR - muL;

while ( (fill-fillL)*(fill-fillR)>0 )
    iter = iter + 1;
    if (fill > fillR)
        muR = muR + wbnd;
        fillR = get_filling(Z,X,ek,WN,nk,beta,muR);
    elseif (fill < fillL)
        muL = muL - wbnd;
        fillL = get_filling(Z,X,ek,WN,nk,beta,muL);
    end
    if (iter > 10)
        error('Self-energy too large, cannot find the chemical potential!')
    end
end

%Now we find mu using Matlab fzero method (super-linear convergence). 
%fzero uses the algorithm that was originated by T. Dekker, which uses a
%combination of bisection, secant, and inverse quadratic interpolation
%methods. See Press et al., Numerical Recipes, Chapter 9.3 Van
%Wijngaarden-Dekker-Brent Method.
options = optimset;
options.TolX = erreps;
out = fzero(@(xmu) get_filling(Z,X,ek,WN,nk,beta,xmu)-fill,[muL muR],options);
