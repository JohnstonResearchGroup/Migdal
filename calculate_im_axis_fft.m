fprintf('\n')
fprintf('**************************************************************\n')
fprintf('  filling0 = %10.6f, T = %10.6f [t], beta = %10.6f [1/t]\n',filling0,T,beta)
fprintf('  w(n=0) = %10.6f [t], number of wn>0 numwi =%6d\n',pi/beta,numwi)
fprintf('Begin imaginary axis calculation.\n')

if (2*numwi) ~= numel(WN)
    error('Wrong Matsubara frequency grid!')
end
%define the freqency grid
WNU = pi*2*(-numwi:numwi-1)/beta;
fmWt = getNCwt(beta,numwi,NCorder,1);   %fermion freq weight
bsWt = getNCwt(beta,numwi,NCorder,0);   %boson freq weight

pre = 1/(nktot*beta);
filling = get_filling(Z,X,EK0,WN,nk,beta,mu);
fprintf('  filling  = %12.8f, mu = %12.8f [t]\n',filling,mu)

%clear the arrays changing size (numwi) at different temperatures
clear('Gn','chic0','Vnm','invD0') 
%array for Fourier transforms
Gn(1:2*nk,1:2*nk,1:2*numwi) = 0;
chic0(1:2*nk,1:2*nk,1:2*numwi) = 0;
Vnm(1:2*nk,1:2*nk,1:2*numwi) = 0;
%array for end point
Gbeta(1:2*nk,1:2*nk) = 0;
%inverse of phonon Green's function D0
invD0 = reshape(repmat(1 + (WNU/wph).^2, [4*nk*nk, 1]), [2*nk, 2*nk, 2*numwi]);
%extend matrix fqph to bosonic freqeuncy dimension
fqphw = repmat(fqph, [1, 1, 2*numwi]);
    

loop = 1;
iter = 0;
fprintf(['  Iter.  _______dZ_______  _______dX_______',...
         '  _______mu_______  _______<n>______\n'])
maxchi_cutoff = 0.9999;
%Anderson acceleration
if (mixing == 1) && (mix_method == 2)
    nptmatZ = 2*nk*2*nk*2*numwi;
    xx(1:nptmatZ,1) = 0;    %store Gnold in one column
    gg(1:nptmatZ,1) = 0;    %store Gn in one column
    mMax = 5;               %mMax>0, 5
    % maximum number of stored residuals (non-negative integer); should not
    % exceed number of rows of xx vectors (so we have a tall/thin matrix in
    % QR decompostion)
    droptol = 1.e10;
    % tolerance for dropping stored residual vectors to improve
    % conditioning: If droptol > 0, drop residuals if the condition number
    % exceeds droptol; if droptol <= 0, do not drop residuals.    
    % dampbt = @(x)0.2*(x<20) + 0.2*(x>=20 && x<25) + 0.5*(x>=25); %0.2
    % dampbt = @(x)1*(x<30) + 0.8*(x>=30);
    dampbt = 1;            
    % damping factor: If dampbt > 0 (and dampbt ~= 1), then the step is damped
    % by dampbt; otherwise, the step is not damped. NOTE: dampbt can be a
    % function handle; form dampbt(iter), where iter is the iteration number
    % and 0 < dampbt(iter) <= 1.
    if (useGqCoupling == 1) && (filling0 >= 0.8)
        AAstart = 20;
    else
        AAstart = 1;    %AAstart>=1, 5
    end
    % acceleration delay factor: If AAstart > 0, start acceleration when
    % iter = AAstart.
    res_hist = [];
    % residual history matrix (iteration numbers and residual norms).
    DG = [];
    % Storage of g-value differences.
    R = []; Q = [];
    mAA = 0;
    % Initialize the number of stored residuals.
end

while (loop == 1)
    iter = iter + 1;

    %updata old self-energies (with mixing)
    if (mixing == 1) && (iter>1)
        if (mix_method == 0)
            Zold = wgt*Z + (1-wgt)*Zold;
            Xold = wgt*X + (1-wgt)*Xold;
        elseif mix_method == 1
            %mix Green's function Gn
            %(a) update Gn with pre-mixed self-energies
            Gn = 1./(1i*Z - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]) - X);
            %(b) mix the Green's function
            Gnold = wgt*Gn + (1-wgt)*Gnold;
            %(c) find the self-energies from the mixed Green's function
            denom = 1./Gnold;
            Zold = imag(denom);
            Xold = -real(denom) - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]);
        else
            if mix_method ~= 2
                error('mix_method = %d is not supported!',mix_method)
            end
        end
        %Recalculate chemical potential mu after mixing. 
        %Note: this may cause convergence problem if the mixed Green's
        %function or self-energies overshoot too much, so it has been 
        %commented out now.
        %mu = get_mu(Z,X,EK0,WN,nk,beta,filling0);
    else
        Zold = Z;
        Xold = X;
    end

    %update Green's function G = Gn, updata old normal self-energy: S
    Gn = 1./(1i*Zold - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]) - Xold);
    if (mix_method == 1) && (mixing == 1) && (iter>1)
        tmp = abs(Gn - Gnold);
        err1 = max(tmp(:));
        if err1 > 1e-12
            fprintf('    %g\n',err1)
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save old Green's function Gn(k,i\omega_n) for mixing
    if (mix_method == 1 || mix_method == 2) && (mixing == 1) && (iter==1)
        Gnold = Gn;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    %compute the self-energy ----------------------------------------------
    %forward Fourier transform Green's function (k,i\omega_n) to (r,\tau)
    Gn = fourier_fwd(Gn,EK0-mu,beta,Nk,numwi,1,1,useSymmetry,runInSerial);   
    Gbeta = -Gn(:,:,1);
    Gbeta(1,1) = Gbeta(1,1) - 1;
    chic0(:,:,1) = Gn(:,:,1).*Gbeta;
    for mm = 2:2*numwi
        chic0(:,:,mm) = Gn(:,:,mm).*Gn(:,:,2*numwi+2-mm);
    end
    %backward Fourier transform susceptibility (r,\tau) to (q,i\omega_n)
    chic0 = real(fourier_bwd(chic0,0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial));
    if (lam0 ~= 0)
        rechic0 = chic0(:,:,numwi+1);
        if useGqCoupling == 1
            maxchi = 2*lam0*max(rechic0(:).*fqph(:));     %2 from spin sum
        else
            maxchi = 2*lam0*max(rechic0(:));     %2 from spin sum
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% User input required here.  The following 'if' statements serve the purpose
%%% of adjusting the denominator cutoff to be incrementally closer to 1 the
%%% temperature approaches the Tc value for the given parameter set.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (T > 0.8)       						%**USER INPUT**
            maxchi_cutoff = 0.995;				%**USER INPUT**
	    elseif (T>0.6)							%**USER INPUT**
            maxchi_cutoff = 0.9995;  			%**USER INPUT**
        elseif (T>0.10)    						%**USER INPUT**
            maxchi_cutoff = 0.99995;			%**USER INPUT**
        else 								
            maxchi_cutoff = 0.999999;		    %**USER INPUT**
        end
        if (maxchi > maxchi_cutoff)
            fprintf('  Warning: max lam0*(2*chic0) = %g,  decrease lam0*(2*chi_max) = %g\n',maxchi,maxchi_cutoff)
            if useGqCoupling == 1
                idx_mask = rechic0.*fqph > maxchi_cutoff/(2*lam0);
                %idx_mask is a logic matrix; it excludes elements for fqph==0
                rechic0(idx_mask) = (maxchi_cutoff/(2*lam0))./fqph(idx_mask);
            else
                rechic0( rechic0 > maxchi_cutoff/( 2*lam0 ) ) = maxchi_cutoff/(2*lam0);
            end
            chic0(:,:,numwi+1) = rechic0;
            
            %fprintf('  Warning: max lam0*(2*chic0) = %g\n',maxchi)
        end     
    end    
    %compute the effective interaction 
    %Vnm = (g_0^2*D_0) / (1 + g_0^2*D_0*2\chi_0^c) 
    %    = (-lam0) / [1 + (\omega_\nu/wph)^2 + (-lam0*2)\chi_0^c]
    %Note: we define an effective interaction Vnm --> -Vnm, where the minus
    %sign is from the interaction term in the Feynman diagram for the
    %self-energy so that the self-energy S = Vnm.*Gn without a minus sign.
    if useGqCoupling == 1
        Vnm = lam0*fqphw./(invD0 + (-lam0*2)*fqphw.*chic0);
    else
        Vnm = lam0./(invD0 + (-lam0*2)*chic0);
    end
    %forward Fourier transform effective interaction (q,i\omega_n) to (r,\tau)
    Vnm = fourier_fwd(Vnm,[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
    S = Vnm.*Gn;
    %backward Fourier transform self-energies (r,\tau) to (k,i\omega_n)
    S = fourier_bwd(S,-Vnm(1,1,1),NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial);
    
    %calculate Z,X; symmetrize +/- frequency; even-k symmetry used
    for mm = (numwi+1):(2*numwi)
        Z(:,:,mm) = WN(mm);
        Z(:,:,mm) = Z(:,:,mm) - imag(S(:,:,mm));
        X(:,:,mm) = real(S(:,:,mm));
    end    
    for mm = (numwi+1):(2*numwi)
        Z(:,:,2*numwi+1-mm) =-Z(:,:,mm);
        X(:,:,2*numwi+1-mm) = X(:,:,mm);
    end
    if checkFreqSymm == 1
        tmp = imag(S);
        errZ = min(min(min(abs(tmp(:,:,1:end) + tmp(:,:,end:-1:1)))));
        tmp = real(S);
        errX = min(min(min(abs(tmp(:,:,1:end) - tmp(:,:,end:-1:1)))));
        if errZ>1e-16 || errX>1e-16
            error('  Error of frequency symmetry is too large! [errZ errX errP] =[%g %g]',...
                errZ,errX)
        end
    end
    %find the chemical potential for the desired filling0
    mu = get_mu(Z,X,EK0,WN,nk,beta,filling0);
    %recomput the filling as a check
    filling = get_filling(Z,X,EK0,WN,nk,beta,mu);  

    % Decide if we have to exit
    tmp = abs(Z-Zold);      diffz  = max(tmp(:));
    tmp = abs(X-Xold);      diffx  = max(tmp(:));       
    if (diffz < eps && diffx < eps)
        loop = 0;
        chi_converge(1) = 1;
    end
    fprintf('  %5d  %16.12f  %16.12f',iter,diffz,diffx)
    fprintf('  %16.12f  %16.12f\n',mu,filling)
    
    %Decide if we need to exit the loop.
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check convergence!\n',maxiter)
        chi_converge(1) = 0;
    end
    
    %Anderson acceleration
    if (loop ~= 0) && (mixing == 1) && (mix_method == 2)
        % Compute the current residual norm.
        %xx = real([Zold(:); Xold(:); Pold(:)]);
        %gval = real([Z(:); X(:); P(:)]);
        
        Gn = 1./(1i*Z - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]) - X);
        xx = Gnold(:);
        gval = Gn(:);

        fval = gval - xx;
        res_norm = max(abs(fval));   % = norm(A,inf)
        res_hist = [res_hist;[iter,res_norm]];

        if iter < AAstart %|| mMax == 0 %note we set mMax>0
            % Without acceleration, update x <- g(x) to obtain the next
            % approximate solution.
            Gnold = wgt*Gn + (1-wgt)*Gnold;
            %find the self-energies from the mixed Green's function
            denom = 1./Gnold;
            Zold = imag(denom);
            Xold = -real(denom) - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]);           
        else
            % Apply Anderson acceleration.
            % Update the df vector and the DG array.
            if iter > AAstart
                df = fval-f_old;
                if mAA < mMax
                    DG = [DG gval-g_old];
                else
                    DG = [DG(:,2:mAA) gval-g_old];
                end
                mAA = mAA + 1;
            end
            f_old = fval;
            g_old = gval;

            if mAA == 0     %iter == AAstart
                % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
                fprintf('  %d, start Anderson acceleration and store %d steps of self-energies\n', iter, mMax);
                Gnold = Gn;
                %find the self-energies from the mixed Green's function
                denom = 1./Gnold;
                Zold = imag(denom);
                Xold = -real(denom) - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]);
            else
                % If mAA > 0, solve the least-squares problem and update the solution.
                if mAA == 1
                    % If mAA == 1, form the initial QR decomposition.
                    R(1,1) = norm(df);
                    Q = R(1,1)\df;
                else
                    % If mAA > 1, update the QR decomposition.
                    if mAA > mMax
                        % If the column dimension of Q is mMax, delete the first column and
                        % update the decomposition.
                        [Q,R] = qrdelete(Q,R,1);
                        mAA = mAA - 1;
                        if size(R,1) ~= size(R,2),
                            Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                        end
                    end
                    % Now update the QR decomposition to incorporate the new column.
                    for j = 1:mAA - 1
                        R(j,mAA) = Q(:,j)'*df;
                        df = df - R(j,mAA)*Q(:,j);
                    end
                    R(mAA,mAA) = norm(df);
                    Q = [Q, R(mAA,mAA)\df];
                end
                if droptol > 0
                    % Drop residuals to improve conditioning if necessary.
                    condDF = cond(R);
                    while condDF > droptol && mAA > 1
                        fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                        [Q,R] = qrdelete(Q,R,1);
                        DG = DG(:,2:mAA);
                        mAA = mAA - 1;
                        % The following treats the qrdelete quirk described above.
                        if size(R,1) ~= size(R,2)
                            Q = Q(:,1:mAA); R = R(1:mAA,:);
                        end
                        condDF = cond(R);
                    end
                end
                % Solve the least-squares problem.
                gamma = R\(Q'*fval);
                % Update the approximate solution.
                xx = gval - DG*gamma;
                % Apply damping if dampbt is a function handle or if dampbt > 0
                % (and dampbt ~= 1).
                if isa(dampbt,'function_handle')
                    xx = xx - (1-dampbt(iter))*(fval - Q*R*gamma);
                else
                    if dampbt > 0 && dampbt ~= 1
                        xx = xx - (1-dampbt)*(fval - Q*R*gamma);
                    end
                end
                
                Gnold(:) = xx(1:nptmatZ);
                %find the self-energies from the mixed Green's function
                denom = 1./Gnold;
                Zold = imag(denom);
                Xold = -real(denom) - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]);
            end
        end %if iter < AAstart
    end %if (mixing == 1) && (mix_method == 2)
end %self-consistency loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now that we have determined the self-energies, we compute the RPA CDW 
%susceptibility \chi^\text{CDW}(q=(\pi,\pi),\omega=0)
%and the RPA pairing susceptibility \chi^\text{SC}(q=0,\omega=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct the electron Green's function
Gn = 1./(1i*Z - repmat(EK0(:,:)-mu, [1, 1, 2*numwi]) - X);
%construct F0 = G(k,i\omega_n)G(-k,-i\omega_n)
F0 = abs(Gn).^2;            %real values

%forward Fourier transform Green's function (k,i\omega_n) to (r,\tau)
Gn = fourier_fwd(Gn,EK0-mu,beta,Nk,numwi,1,1,useSymmetry,runInSerial);
Gbeta = -Gn(:,:,1);
Gbeta(1,1) = Gbeta(1,1) - 1;
chic0(:,:,1) = Gn(:,:,1).*Gbeta;
chic0(:,:,2:end) = Gn(:,:,2:end).*Gn(:,:,end:-1:2);

%backward Fourier transform susceptibility (r,\tau) to (q,i\omega_n)
chic0 = real(fourier_bwd(chic0,0,NCorder,bsWt,beta,Nk,numwi,0,useSymmetry,runInSerial));
if useGqCoupling == 1
    P1 = invD0 + (-lam0*2)*fqphw.*chic0;
else
    P1 = invD0 + (-lam0*2)*chic0;
end
P2 = zeros(size(P1));

%compute the effective interaction
if useGqCoupling == 1
    Vnm = lam0*fqphw./(invD0 + (-lam0*2)*fqphw.*chic0);
else
    Vnm = lam0./(invD0 + (-lam0*2)*chic0);
end

%forward Fourier transform effective interaction (q,i\omega_n) to (r,\tau)
Vnm = fourier_fwd(Vnm,[],beta,Nk,numwi,0,0,useSymmetry,runInSerial);
if useGqCoupling == 1
    chi = 2*chic0(:,:,numwi+1)./(1 + (-lam0*2)*fqph.*chic0(:,:,numwi+1));
else
    chi = 2*chic0(:,:,numwi+1)./(1 + (-lam0*2)*chic0(:,:,numwi+1));
end

fprintf('\n')
fprintf('Begin pairing vertex calculation.\n')
fprintf(['  Iter.  _______dTT______', '\n'])
iter = 0;
loop = 1;
TT = ones(2*nk,2*nk,2*numwi);
TTold = zeros(size(TT));

%Anderson acceleration
if (mixing == 1) && (mix_method == 2)
    nptmatZ = 2*nk*2*nk*2*numwi;
    xx(1:nptmatZ,1) = 0;    %store Gnold in one column
    gg(1:nptmatZ,1) = 0;    %store Gn in one column
    mMax = 5;               %mMax>0, 5
    droptol = 1.e10;
    dampbt = 1;
    % damping factor: If dampbt > 0 (and dampbt ~= 1), then the step is damped
    % by dampbt; otherwise, the step is not damped. NOTE: dampbt can be a
    % function handle; form dampbt(iter), where iter is the iteration number
    % and 0 < dampbt(iter) <= 1.
    AAstart = 1;    %AAstart>=1, 5
    % acceleration delay factor: If AAstart > 0, start acceleration when
    % iter = AAstart.
    res_hist = [];
    % residual history matrix (iteration numbers and residual norms).
    DG = [];
    % Storage of g-value differences.
    R = []; Q = [];
    mAA = 0;
    % Initialize the number of stored residuals.
end

while (loop == 1)
    iter = iter + 1;
           
    %updata vertex function (with mixing)
    if (mixing == 1) && (iter>1)
        if (mix_method == 0) || (mix_method == 1)            
            TTold = wgt*TT + (1-wgt)*TTold;
        else
            if mix_method ~= 2
                error('mix_method = %d is not supported!',mix_method)
            end
        end
    else
        TTold = TT;
    end    
    
    %forward Fourier transform F_0.*\Lambda (k,i\omega_n) to (r,\tau)
    tmp = fourier_fwd(F0.*TTold,[],beta,Nk,numwi,0,1,useSymmetry,runInSerial);
    %backward Fourier transform Vnm.*F_0.*\Lambda (r,\tau) to (k,i\omega_n)
    tmp = real(fourier_bwd(Vnm.*tmp,0,NCorder,fmWt,beta,Nk,numwi,1,useSymmetry,runInSerial));
    TT = 1 + tmp;
    tmp = abs(TT-TTold);   diffT = max(tmp(:));
    fprintf('  %5d  %16.12f\n',iter,diffT)
    %decide if we exit
    if (diffT < eps)
        loop = 0; 
        chi_converge(2) = 1;
    end
    if (iter >= maxiter)
        loop = 0;
        fprintf('  Maximal iterations =%5d reached. Check pairing vertex convergence!\n',maxiter)
        chi_converge(2) = 0;
    end
    
    %Anderson acceleration
    if (loop ~= 0) && (mixing == 1) && (mix_method == 2)
        % Compute the current residual norm.
        xx = TTold(:);
        gval = TT(:);

        fval = gval - xx;
        res_norm = max(abs(fval));   % = norm(A,inf)
        res_hist = [res_hist;[iter,res_norm]];

        if iter < AAstart %|| mMax == 0 %note we set mMax>0
            % Without acceleration, update x <- g(x) to obtain the next
            % approximate solution.
            TTold = wgt*TT + (1-wgt)*TTold;
        else
            % Apply Anderson acceleration.
            % Update the df vector and the DG array.
            if iter > AAstart
                df = fval-f_old;
                if mAA < mMax
                    DG = [DG gval-g_old];
                else
                    DG = [DG(:,2:mAA) gval-g_old];
                end
                mAA = mAA + 1;
            end
            f_old = fval;
            g_old = gval;

            if mAA == 0     %iter == AAstart
                % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
                fprintf('  %d, start Anderson acceleration and store %d steps\n', iter, mMax);
                TTold = TT;
            else
                % If mAA > 0, solve the least-squares problem and update the solution.
                if mAA == 1
                    % If mAA == 1, form the initial QR decomposition.
                    R(1,1) = norm(df);
                    Q = R(1,1)\df;
                else
                    % If mAA > 1, update the QR decomposition.
                    if mAA > mMax
                        % If the column dimension of Q is mMax, delete the first column and
                        % update the decomposition.
                        [Q,R] = qrdelete(Q,R,1);
                        mAA = mAA - 1;
                        if size(R,1) ~= size(R,2),
                            Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                        end
                    end
                    % Now update the QR decomposition to incorporate the new column.
                    for j = 1:mAA - 1
                        R(j,mAA) = Q(:,j)'*df;
                        df = df - R(j,mAA)*Q(:,j);
                    end
                    R(mAA,mAA) = norm(df);
                    Q = [Q, R(mAA,mAA)\df];
                end
                if droptol > 0
                    % Drop residuals to improve conditioning if necessary.
                    condDF = cond(R);
                    while condDF > droptol && mAA > 1
                        fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, mAA-1);
                        [Q,R] = qrdelete(Q,R,1);
                        DG = DG(:,2:mAA);
                        mAA = mAA - 1;
                        % The following treats the qrdelete quirk described above.
                        if size(R,1) ~= size(R,2)
                            Q = Q(:,1:mAA); R = R(1:mAA,:);
                        end
                        condDF = cond(R);
                    end
                end
                % Solve the least-squares problem.
                gamma = R\(Q'*fval);
                % Update the approximate solution.
                xx = gval - DG*gamma;
                % Apply damping if dampbt is a function handle or if dampbt > 0
                % (and dampbt ~= 1).
                if isa(dampbt,'function_handle')
                    xx = xx - (1-dampbt(iter))*(fval - Q*R*gamma);
                else
                    if dampbt > 0 && dampbt ~= 1
                        xx = xx - (1-dampbt)*(fval - Q*R*gamma);
                    end
                end

                TTold(:) = xx(1:nptmatZ);
            end
        end %if iter < AAstart
    end %if (mixing == 1) && (mix_method == 2)
end %vertex equation loop

%(iii) compute the pair-field susceptibility
tmp = F0.*TT;
chisc = sum(real(tmp(:)))/(4*nk*nk*beta);   %q=(0,0), \omega=0

