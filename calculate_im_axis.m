fprintf('\n')
fprintf('**************************************************************\n')
fprintf('  filling0 = %10.6f, T = %10.6f [t], beta = %10.6f [1/t]\n',filling0,T,beta)
fprintf('  w(n=0) = %10.6f [t], number of wn>0 numwi =%6d\n',pi/beta,numwi)
fprintf('Begin imaginary axis calculation.\n')

pre = lam0/(beta*4*nk*nk);
%define temporary array for Fourier transformed quantities
den = zeros(2*nk,2*nk);
gZ = zeros(2*nk,2*nk);
gX = zeros(2*nk,2*nk);
GZ = zeros(2*nk,2*nk,2*numwi);
GX = zeros(2*nk,2*nk,2*numwi);
gP1 = zeros(2*nk,2*nk);
gP2 = zeros(2*nk,2*nk);
GP1 = zeros(2*nk,2*nk,2*numwi+1);
GP2 = zeros(2*nk,2*nk,2*numwi+1);


loop = 1;
iter = 0;
fprintf(['  Iter.  _______dZ_______  _______dX_______'...
         '  _______dP1______  _______dP2______',...
         '  _______mu_______  _______<n>______',...
         '\n'])
while (loop == 1)
    iter = iter + 1;

    %mix new and old self-energy to speed up convergence
    if iter>1
        P1 = P1*wgt + (1-wgt)*P1old;
        P2 = P2*wgt + (1-wgt)*P2old;
        Z = Z*wgt + (1-wgt)*Zold;
        X = X*wgt + (1-wgt)*Xold;
    end

    %store self-energies from the previous iteration (after mixing)
    Zold = Z;
    Xold = X;
    P1old = P1;
    P2old = P2;

    %set \xi_k = \epsilon_k - mu
    EK = EK0 - mu;
    %compute the momentum grid FFTs
    for n = 1:length(NU)
        if (n ~= length(NU))
            den(:,:) = Zold(:,:,n).^2 + (EK(:,:)+Xold(:,:,n)).^2;
            gZ(:,:) = Zold(:,:,n)./den;
            gX(:,:) = (EK(:,:)+Xold(:,:,n))./den;
            GZ(:,:,n) = fft2(gZ);
            GX(:,:,n) = fft2(gX);
        end
        den(:,:) = P1old(:,:,n).^2 + P2old(:,:,n).^2;
        gP1(:,:) = P1old(:,:,n)./den;
        gP2(:,:) = P2old(:,:,n)./den;
        GP1(:,:,n) = fft2(gP1);
        GP2(:,:,n) = fft2(gP2);
    end

    %calculate the phonon self-energy
    P1 = zeros(2*nk,2*nk,2*numwi+1);
    P2 = zeros(2*nk,2*nk,2*numwi+1);
    for n = 1:length(NU)
        P1(:,:,n) = 1 + (NU(n)/wph)^2;
    end
    for nu = -numwi:numwi
        for n = -numwi:numwi-1
            m = n + nu;
            %\omega_{m'} = \omega_m + \nu_n
            if (m < numwi && m >= -numwi)
                P1(:,:,nu+numwi+1) = P1(:,:,nu+numwi+1) ...
                    + 2*pre*ifft2(GX(:,:,n+numwi+1).*(GX(:,:,m+numwi+1))) ...
                    - 2*pre*ifft2(GZ(:,:,n+numwi+1).*(GZ(:,:,m+numwi+1)));
                P2(:,:,nu+numwi+1) = P2(:,:,nu+numwi+1) ...
                    + 2*pre*ifft2(GX(:,:,n+numwi+1).*(GZ(:,:,m+numwi+1))) ...
                    + 2*pre*ifft2(GZ(:,:,n+numwi+1).*(GX(:,:,m+numwi+1)));
            end
        end
    end

    %calculate the electron self-energy
    for n = 1:2*numwi
        Z(:,:,n) = WN(n);
    end
    X = zeros(2*nk,2*nk,2*numwi);
    for n = -numwi:numwi-1
        wn = WN(n+numwi+1);
        for m = -numwi:numwi-1
            nu = n - m;
            %\nu_n = \omega_m - \omega_{m'}
            if (nu >= -numwi & nu <= numwi)
                Z(:,:,n+numwi+1) = Z(:,:,n+numwi+1) ...
                    + pre*ifft2(GZ(:,:,m+numwi+1).*GP1(:,:,nu+numwi+1)) ...
                    - pre*ifft2(GX(:,:,m+numwi+1).*GP2(:,:,nu+numwi+1));
                X(:,:,n+numwi+1) = X(:,:,n+numwi+1) ...
                    - pre*ifft2(GX(:,:,m+numwi+1).*GP1(:,:,nu+numwi+1)) ...
                    - pre*ifft2(GZ(:,:,m+numwi+1).*GP2(:,:,nu+numwi+1));
            end
        end
    end

    %remove small imaginary part from Z, X, P1, P2 due to fft
    Z = real(Z);
    X = real(X);
    P1 = real(P1);
    P2 = real(P2);

    %find the chemical potential for the desired filling0
    mu = get_mu(Z,X,EK0,WN,nk,beta,filling0);
    %recomput the filling as a check
    filling = get_filling(Z,X,EK0,WN,nk,beta,mu);

    tmp = abs(P1-P1old);    diffp1 = max(tmp(:));
    tmp = abs(P2-P2old);    diffp2 = max(tmp(:));
    tmp = abs(Z-Zold);      diffz  = max(tmp(:));
    tmp = abs(X-Xold);      diffx  = max(tmp(:));
    if (max([diffp1 diffp2 diffz diffx]) < eps), loop = 0; end
    fprintf('  %5d  %16.12f  %16.12f  %16.12f  %16.12f',iter,diffz,diffx,diffp1,diffp2)
    fprintf('  %16.12f  %16.12f\n',mu,filling)
    %decide if we exit
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check self-energy convergence!\n',maxiter)
        %continue
    end
end %self-energy loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%compute the RPA CDW susceptibility \chi^\text{CDW}(q=(\pi,\pi),\omega=0)
pol = (wph*wph)*(P1(:,:,numwi+1) - 1 + 1i*P2(:,:,numwi+1));
pol = real(pol);
chi0 = -pol/(alpha*alpha);
chi = chi0./(1-lam0*chi0);

%compute the RPA pairing susceptibility \chi^\text{SC}(q=0,\omega=0)
%(i) construct the Green's function
Polarization = zeros(2*nk,2*nk);
Sigma = zeros(2*nk,2*nk);
G = zeros(2*nk,2*nk,2*numwi);
D = zeros(2*nk,2*nk,2*numwi+1);
%update \xi_k = \epsilon_k - mu with converged mu
EK = EK0 - mu;
for n = -numwi:numwi
    Polarization(:,:) = wph*wph*(-1 - (NU(n+numwi+1)/wph)^2 ...
        + P1(:,:,n+numwi+1) + 1i*P2(:,:,n+numwi+1));
    D(:,:,n+numwi+1) = 1./(-wph*wph - NU(n+numwi+1)^2 - Polarization(:,:));
    if (n ~= numwi)
        Sigma = 1i*(WN(n+numwi+1)-Z(:,:,n+numwi+1)) + X(:,:,n+numwi+1);
        G(:,:,n+numwi+1) = 1./(1i*WN(n+numwi+1) - EK(:,:) - Sigma);
    end
end

%(ii) construct F0 = G(k,i\omega_n)G(-k,-i\omega_n)
F0 = G(:,:,:).*G(:,:,end:-1:1);
%solve vertex equation by iterations
fprintf('\n')
fprintf('Begin pairing vertex calculation.\n')
fprintf(['  Iter.  _______dTT______', '\n'])
iter = 0;
loop = 1;
mat1 = zeros(2*nk,2*nk);
mat2 = zeros(2*nk,2*nk);
TT = ones(2*nk,2*nk,2*numwi);
TTold = zeros(size(TT));
while (loop == 1)
    iter = iter + 1;
    TTold = TT;
    TT = ones(2*nk,2*nk,2*numwi);
    for n = -numwi:numwi-1
        for np = -numwi:numwi-1
            nnu = n-np;
            if(abs(nnu) <= numwi)
                mat1(:,:) = fft2(D(:,:,nnu+numwi+1));
                mat2(:,:) = fft2(TTold(:,:,np+numwi+1).*F0(:,:,np+numwi+1));
                TT(:,:,n+numwi+1) = TT(:,:,n+numwi+1) ...
                    - (alpha*alpha/(beta*4*nk*nk))*ifft2(mat1.*mat2);
            end
        end
    end
    tmp = abs(TT-TTold);   diffT = max(tmp(:));
    fprintf('  %5d  %16.12f\n',iter,diffT)
    %decide if we exit
    if (diffT < eps), loop = 0; end
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check pairing vertex convergence!\n',maxiter)
        %continue
    end
end %vertex equation loop

%(iii) compute the pair-field susceptibility
tmp = F0.*TT;
chisc = sum(real(tmp(:)))/(4*nk*nk*beta);   %q=(0,0), \omega=0

