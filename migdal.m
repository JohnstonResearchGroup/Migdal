% function migdal(nk,lam0,wph)
function migdal(nk,lam0,wph,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% migdal.m  
%--------------------------------------------------------------------------
% Revision History
%--------------------------------------------------------------------------
% 23 October 2018. 
% Version 1.0 
% Initial release of the code, developed by Yan Wang, Phil Dee, and 
% Steven Johnston. The code is maintained by Steven Johnston. Any bug 
% fixes can be reported to sjohn145@utk.edu. 
% 
% Documentation for the code can be found in the preprint paper 
% arXiv:XXXX.XXXXX (2018) or in (INSERT REFERENCE ON PUBLICATION). 
% If you use this code for your research, please cite this paper.
%
% This code and its its documentation is released “as is” without any 
% warranty, and thus with no additional responsibility or liability. 
% 
% All lines ending with "**USER INPUT**" set control or model parameters 
% that should be supplied by the user. 
%--------------------------------------------------------------------------
% Input: nk - number of momentum points in the k_x direction.  The number 
%             of point in the k_y direction is taken to be the same. 
%        lam0 - The bare value of the dimensionless electron-phonon
%               coupling. 
%        wph - The bare value of the (dispersionless)  Phonon mode. 
%        varargin - Input arguements to determine the behavior of the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc             %Clear the Screen
clear           %Clear memory
close all       %Close all prior Figure windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the momentum dependence of the electron-phonon coupling
gqnum = 0;				 %**USER INPUT**
if gqnum == 0            %Holstein isotropic coupling
 varargin = [];
elseif gqnum == 1        %Buckling mode g(q)~cos(q/2)^2 [CDW q ~ (0,0)]
 varargin = {1,'bk'};
elseif gqnum == 2        %Breathing mode g(q) ~sin(q/2)^2 [CDW q ~(pi,pi)]
 varargin = {1,'br'};    %               (more difficult to converge)
elseif gqnum == 3        %Forward scattering, Exponential form g(q) ~ exp(-q/q0)
 varargin = {1,'exp'};
elseif gqnum == 5        %Forward scattering, Gaussian form g(q) ~ exp(-(q/q0)^2)
 varargin = {1,'gauss'};
elseif gqnum == 5        %Forward scattering, Kronecker-delta function
 varargin = {1,'kd'};
elseif gqnum == 6        %Forward scattering, box function
 varargin = {1,'window'}; 
else
 error('The value of gqnum does not correspond to an implemented type.');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize parameters (set 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 1;			      % NN-hopping  **USER INPUT**		
tp = 0.0; 		      % NNN-hopping parameter **USER INPUT** 
mu = 0;			      % Initial guess for chemical potential **USER INPUT**		
Wband = 8*t;          % 2D Electronic Bandwidth t,t' model: 8t **USER INPUT**	
wph = 0.3*t;  	      % Phonon frequency **USER INPUT**	
nk = 8;		          % Half of linear dimension # of k-points **USER INPUT**	
lam0 = 0.3*Wband; 	  % bare e-ph coupling (Marsiglio's, [Energy]) **USER INPUT**	
wc = 80*Wband;		  % Matsubara Frequency cutoff.  wc>> Wband **USER INPUT**
Nk = [nk nk];     	  % One quadrant of K-mesh **USER INPUT**	
nktot = 4*Nk(1)*Nk(2);% Total # of k-points 				

%%% Note on usage of the code in loops:
% When running over a list of different lam0 values, uncomment the lines 
% below and specify the range of values to use. You also need to uncomment 
% the end statement for the for loop (tagged with "% Lambda Loop end")

% vlam0 = linspace(0.1*Wband, 0.3*Wband, 10); %0.3*Wband;
% nlam0 = numel(vlam0);
%%% Lambda Loop start %%%
% for nl = 10:10   %nlam0 %**USER INPUT**	
% lam0 = vlam0(nl);		  %**USER INPUT**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments on the relationship of couplings
% g_0^2 = \hbar*\alpha^2/(M*2*\omega_E) = \lambda_0*(\hbar*\omega_E)/2, units: [E]^2.
% \lambda_0 = \alpha^2/(M*\omega_E^2) = 2*g_0^2/(\hbar*\omega_E), units: [E].
% K = M*\omega_E^2, units: [E]/[L]^2. \alpha, units: [E]/[L].
% in the code and in Marsiglio's paper:
% \hbar\alpha/\sqrt(M) --> alpha [E]^{3/2}; \hbar*\omega_E --> wph; \lambda_0 --> lam0
% g_0^2 = alpha^2/(2*wph) = (lam0*wph)/2; lam0 = alpha^2/wph^2 = 2*g_0^2/wph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g0 = sqrt(lam0*wph/2);							%**USER INPUT**
alpha = wph*sqrt(lam0);							%**USER INPUT**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfill = 0.86;								                    %**USER INPUT**
% Comment the above line and uncomment next 5 lines to create a loop over 
% filling values   
%nstart = 0.80;								                    %**USER INPUT**
%nfinal = 0.90;   % left off lam=0.21111 and nstart = 0.68		%**USER INPUT**
%deltanf = 0.02;							                    %**USER INPUT**
%nfpts  = round((nfinal-nstart)/deltanf,0) + 1 ;			    %**USER INPUT**
%nfill = linspace(nstart, nfinal, nfpts ); %linspace(0.04, 1.0, 49); 	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(varargin)
    useGqCoupling = 0;
else
    if isnumeric(varargin{1})
        useGqCoupling = varargin{1};
        if useGqCoupling == 1
            vertex_type = varargin{2};
        end
    else
        error('Wrong Input!');
    end
end
% For the following parameters 1 = true, 0 = false. 
restart = 1;		%**USER INPUT** Is this a restart of a prior run?
FFT_freq = 1;       %**USER INPUT** Use the FFT rout.	
flySelfE = 1;       %**USER INPUT** load self-energy calculated on the fly		
runInSerial = 1;    %**USER INPUT** switch between for-loop and parfor-loop		
useSymmetry = 1;    %**USER INPUT** C_2(inversion) symmetry is assumed
NCorder = 2;		%**USER INPUT** 
Power2NumFreq = 0;	%**USER INPUT** Is the number of Matsubara freq equal to a power of 2?

% Set the mixing method. 
%mix_method: 0, mix self-energy; 
%            1, mix Green's function; 
%            2, mix Green's function and use Anderson acceleration
mixing = 1.0;		%**USER INPUT** 
mix_method = 2;		%**USER INPUT** 


if useGqCoupling == 1
    maxiter = 300;	%**USER INPUT** Maximum number of iterations
    wgt = 0.05;		%**USER INPUT** Weight on mixing
else
    maxiter = 300; 	%**USER INPUT** Maximum number of iterations
    wgt = 0.6;		%**USER INPUT** Weight on mixing
end

eps = 1e-8;			%**USER INPUT** Convergence Criterion

checkFreqSymm = 0;  %check frequency even/odd symmetry
if useGqCoupling == 1 && FFT_freq ~= 1
    error('Momentum dependent coupling not added non-fft code calculate_im_axis!')
end

printchiq = 0;      %**USER INPUT** Set to 1 if you want to print Xcdw(q).  		
printGfunc = 0;     %**USER INPUT** Set to 1 if you want to print Re[-G(r=0,tau=beta/2)]. 
checkFreqSymm = 0;  %**USER INPUT** check frequency even/odd symmetry			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Use arrays to store the magnitude of temperature increments, maximum, 
%%% and minimum range values for each dTi(i). This storage method is used 
%%% because certain temperature ranges may require finer temperature steps. 
dTi =  [1.00 0.1 0.01 0.001 0.0005 0.0001 0.00001 0.000001]; %**USER INPUT** 
Tmax = [2.00 0.9 0.29 0.199 0.0995 0.0699 0.03949 0.000009]; %**USER INPUT** 
Tmin = [1.00 0.3 0.20 0.100 0.0700 0.0495 0.00001 0.000001]; %**USER INPUT** 
  
vT = [];
for nt = 1:length(dTi)
  vT = [vT,[Tmax(nt):-dTi(nt):Tmin(nt)]];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the the file IO locations  **USER INPUT**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileDir = './';
filenamestr = ['_nkEq' num2str(nk) '_lam0Eq' num2str(lam0) ...
               '_wphEq' num2str(wph) '_nEq' num2str(nfill) ...
               '_Phase_g_q_' num2str(gqnum) '_tpEq' num2str(abs(tp)) ...
               '_meff.dat'];
fileout = ['out' filenamestr];
filechi = ['chi' filenamestr];
filechiq = ['chiq' filenamestr];
fileGfunc = ['Gfunc' filenamestr];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the momentum grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = (0:1:2*nk-1)*pi/nk;
[KX,KY] = meshgrid(K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the bare band dispersion (without the chemical potential).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EK0 = -2*t*(cos(KX)+cos(KY)) - 4*tp*cos(KX).*cos(KY);
EK = zeros(size(EK0));

%set up the q-dependence fqph of coupling vertex |g(q)|^2 = g_0^2*fqph =
%(0.5*lam0*wph)*fqph
if useGqCoupling == 1  
    vertexinput = {vertex_type,nk};
    fqph = vertexq(KX,KY,vertexinput{:});
else
    vertexinput = [];
    fqph = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define empty vectors to collect the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vpol = [];
vchi = [];
vchisc = [];
vmu = [];

if restart == 1
    %open the output files to append and write; create if necessary
    fid = fopen([fileDir,filechi],'a');
    fidq = fopen([fileDir,filechiq],'a');
    fidG = fopen([fileDir,fileGfunc],'a');
else
    %open to write; discard any old data; create if necessary
    fid = fopen([fileDir,filechi],'w');
    fidq = fopen([fileDir,filechiq],'w');
    fidG = fopen([fileDir,fileGfunc],'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start outputting to the screen and output file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidout = fopen([fileDir, fileout],'a');
fclose(fidout);
diary([fileDir, fileout]);
fprintf('\n')
fprintf('Susceptibility calculation\n')
fprintf('  grid nk = %4d, convergence criterion = %g\n',nk,eps)
fprintf('  lambda_0 = %g [t], g_0 = %g [t], w_ph = %g [t],\n',lam0,g0,wph)
fprintf('  mixing = %d, weight = %g, mix_method = %d\n',mixing,wgt,mix_method)
fprintf('  maxiter = %d\n',maxiter)
fprintf('  useSymmetry = %1d\n',useSymmetry)
fprintf('  NCorder = %1d\n',NCorder)
fprintf('  restart = %1d\n',restart)
fprintf('  flySelfE = %1d\n',flySelfE)
fprintf('  FFT_freq = %1d\n',FFT_freq)

if useGqCoupling == 1
    fprintf('  useGqCoupling = %1d,',useGqCoupling)
    if strcmp(vertexinput{1},'exp') || strcmp(vertexinput{1},'gauss')
        fprintf('  vertex f(q) = ''%s'', q0 = %g\n',vertexinput{1},q0)
    else
        fprintf('  vertex f(q) = ''%s''\n',vertexinput{1})
    end
end

start_time = datestr(now);
fill_done = [];
T_done = [];
if restart == 1
    fprintf(['Restart calculation ' start_time '\n']);
    sfile = dir([fileDir,filechi]);
    if sfile.bytes == 0
        tmp = [];
        this_fill_done = 0;
    else
        tmp = dlmread([fileDir,filechi]);
        fill_done = tmp(:,1);
        T_done = tmp(:,2);
        vthis_fill_done = tmp(:,end);
    end
else
    this_fill_done = 0;
    fprintf(['Start new calculation ' start_time '\n'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vfill = nfill; 

npt_T = numel(vT);
npt_fill = numel(vfill);
npt_this_run = 0;
vtm = [];
npt_to_run = npt_T*npt_fill;
chi_converge = [1 1];
chi_negative = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start Filling anf Temperature loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nnf = 1:npt_fill
    filling0 = vfill(nnf);
    for nnT = 1:npt_T
        tic;
        npt_to_run = npt_to_run - 1;
        T = vT(nnT);
        beta = 1/T;        
        if (restart == 1) && (npt_this_run == 0) && ~isempty(fill_done)
            %check if data point (filling0, T) has been done
            idx1 = find(abs(filling0 - fill_done)<2e-6);
            idx2 = find(abs(T - T_done)<2e-6);
            this_fill_done = ~isempty(find(vthis_fill_done(idx1)==1));            
            point_done = ~isempty(intersect(idx1,idx2));
            if this_fill_done
                if (nnT == npt_T)
                    fprintf('  filling0 = %10.6f done. Go to next filling.\n',filling0)
                end
                continue
            else
                if point_done
                    fprintf('  filling0 = %10.6f, T = %10.6f [t], beta = %10.6f [1/t] done. Go to next T.\n',filling0,T,beta)
                    continue
                else
                    fprintf('  filling0 = %10.6f, T = %10.6f [t], beta = %10.6f [1/t] calculation restarting.\n',filling0,T,beta)
                    npt_this_run = 1;
                    chi_converge = [1 1];
                    chi_negative = [0 0];
                    this_fill_done = 0;
                end
            end
        else
            if (nnf == 1) && (nnT == 1)
                npt_this_run = 1;
            end
        end

        if this_fill_done && (nnT>1)
            vtm(npt_this_run) = toc;
            if (nnT == npt_T)
                fprintf('  filling0 = %10.6f done. Go to next filling.\n',filling0)
                fprintf('Time = %.2f s. Total Time = %.2f s.\n',vtm(npt_this_run),sum(vtm))
                fprintf('Done %6d points. To go %6d points.\n',npt_this_run,npt_to_run)
                fprintf('**************************************************************\n')
            end
            npt_this_run = npt_this_run + 1;
            continue
        end        
        
        if (flySelfE == 1) && (nnT > 1) && (npt_this_run > 1)
            WNold = WN;
            WNUold = WNU;
            NUold = NU;
        end
                                
        %Now we define the frequency grids.
        numwi = round(wc/(pi/beta));
        numwi = max(numwi,8);
        WN = (pi/beta)*(2*(-numwi:numwi-1)+1);
        WNU = pi*2*(-numwi:numwi-1)/beta;
        NU = (pi/beta)*2*(-numwi:numwi);
        %Note: WNU and NU are both bosonic Matsubara frequencies, but the
        %number of WNU is the same as the WN (for fermion). When
        %frequency-time FFT is used (FFT_freq=1), WNU should be used for
        %boson grid so that the time space grid is same for boson and
        %fermion.
        
        % Z stores the value of w_n*Z(k,w_n).
        % X stores the value of  chi(k,w_n).
        % P1 stores the real part of Pi(k,w_n).
        % P2 stores the imaginary part of Pi(k,w_n).

        if (flySelfE == 1) && (nnT > 1) && (npt_this_run > 1)
            %do not interpolate wn*Z term, the extrapolate error near the
            %tail can be large            
            for nn = 1:numel(WNold)
                Z(:,:,nn) = Z(:,:,nn)/WNold(nn);
                %\omega_n*Z --> Z used in interpSelfE
            end            
            if FFT_freq == 1
                [Z X P1 P2 mu] = interpSelfE(WN,[],1,Z,X,[],[],WNold,[],mu);                
            else
                %do not interpolate P1 including NU/wph term, the extrapolate
                %error near the tail can be large
                for nn = 1:length(NUold)
                    P1(:,:,nn) = P1(:,:,nn) - (1 + (NUold(nn)/wph)^2);
                end
                P2(:) = 0;
                [Z X P1 P2 mu] = interpSelfE(WN,NU,1,Z,X,P1,P2,WNold,NUold,mu);
                for nn = 1:length(NU)
                    P1(:,:,nn) = P1(:,:,nn) + (1 + (NU(nn)/wph)^2);
                end
                P2(:) = 0;
            end
            for nn = 1:numel(WN)
                Z(:,:,nn) = Z(:,:,nn)*WN(nn);
                %Z --> \omega_n*Z used calculation
            end            
        else
            %define electron self-energy array            
            Z = ones(2*nk,2*nk,2*numwi);
            for nn = 1:2*numwi
                Z(:,:,nn) = WN(nn);
                %Z --> \omega_n*Z used calculation
                S = zeros(2*nk,2*nk,2*numwi);
            end            
            X = zeros(2*nk,2*nk,2*numwi);
            Xold = zeros(2*nk,2*nk,2*numwi);
            Zold = zeros(2*nk,2*nk,2*numwi);
            %define phonon self-energy array
            if FFT_freq == 1
                P1 = [];
                P2 = [];
                P1old = [];
                P2old = [];
            else
                P1 = zeros(2*nk,2*nk,2*numwi+1);        %real part
                P2 = zeros(2*nk,2*nk,2*numwi+1);        %imag part
                P1old = zeros(2*nk,2*nk,2*numwi+1);
                P2old = zeros(2*nk,2*nk,2*numwi+1);
                %initialize the phonon self-energy (real part)
                for nn = 1:length(NU)
                    P1(:,:,nn) = 1 + (NU(nn)/wph)^2;
                end
            end
        end
        
        %find the chemical potential for the desired filling0
        mu = get_mu(Z,X,EK0,WN,nk,beta,filling0);
        
        if FFT_freq == 1            
            calculate_im_axis_fft;            
        else
            calculate_im_axis;
        end

	    %calculate the Fermi-surface average of Z(kx,ky,wn=pi/beta) 
	    Zavg = sum(sum(Z(:,:,numwi+1).*gauss(EK0-mu,0.01)))/sum(sum(gauss(EK0-mu,0.01)));  
	    % calculate LDOS, using the real part of the green's function 
	    ldos = -real(Gn(1,1,numwi+1))/T; 

        %store data in column-vector
        vmu = [vmu; mu];
        [chi_max idx] = max(chi(:));
        qx_max = KX(idx)/pi;
        qy_max = KY(idx)/pi;
        vchi = [vchi; [chi(1,1), chi(nk+1,nk+1), qx_max, qy_max, chi_max] ];              
        %q=(0,0), (pi,pi), q_max \omega=0
        vpol = [vpol; P1(nk+1,nk+1,numwi+1)];
        vchisc = [vchisc; chisc];

        %check if chi_CDW or chi_sc diverges and decide to exit early
        chi_negative = [(min(chi(:))<0) (chisc<0)];
        this_fill_done = (chi_converge(1)==0 || chi_converge(2)==0 || chi_negative(1)==1 || chi_negative(2)==1);  
        
        mushift0 = lam0;
        S_H = -lam0*filling0;
        %format:10 field width + 1 space
        fprintf(fid,'%10.6f ',lam0./Wband, filling0, T, vpol(end));
        %format:(continue from the above) 12 field width + 1 space
        fprintf(fid,'%12.6f ',vchi(end,:), vchisc(end));
        %format:(continue from the above) 10 field width + 1 space
        fprintf(fid,'%10.6f ',vmu(end), S_H, mushift0, (vmu(end) + mushift0 + S_H), (Zavg/(pi*T)-1), ldos);
        %format:(continue from the above) + end of line
        fprintf(fid,'%5.3f %1d %1d %1d\n',wgt,chi_negative(1),chi_negative(2),this_fill_done);
        
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print the q-dependence of  Xcdw 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These lines can be used to print along the (0,0)-->(pi,pi) direction
%         for ii = 1:2*nk            
%              fprintf(fidq,'%5.3f %5.3f %12.6f %12.6f %12.6f \n',vT(nnT), vfill(nnf), KX(1,ii)/pi, KY(ii,1)/pi, chi(ii,ii)); 
%         end
% These lines can be used to print all of the points. Note that this output 
% file with 'chiq' in the filename stores this data, and it can be a very 
% large file for nk>32.  Opening this file directly is not advised.  
        if printchiq == 1
           for jj = 1:2*nk
              for ii = 1:2*nk
            	 fprintf(fidq,'%5.3f %5.3f %12.6f %12.6f %12.6f \n',...
                         vT(nnT), vfill(nnf), KX(1,ii)/pi, KY(jj,1)/pi,...
                         chi(ii,jj)); 
              end
           end    
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To output the Real part of the Matsubara GF, use the routine below.
        ReGn = zeros(2*numwi);
        if printGfunc == 1
            for nn = 1:2*numwi
                ReGn(nn) = real(-Gn(1,1,nn));
            end
            fprintf(fidG,'%5.6f %5.3f %12.6f \n',vT(nnT), vfill(nnf),ReGn(numwi));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		

        vtm(npt_this_run) = toc;
        fprintf('\n')
        fprintf('Done filling0 = %10.6f, T = %10.6f [t], beta = %10.6f [1/t]\n',filling0,T,beta)        
        fprintf('Time = %.2f s. Total Time = %.2f s.\n',vtm(npt_this_run),sum(vtm))
        fprintf('Done %6d points. To go %6d points.\n',npt_this_run,npt_to_run)
        fprintf('**************************************************************\n')
        npt_this_run = npt_this_run + 1;
        if nnT == npt_T
            chi_converge = [1 1];
            chi_negative = [0 0];
            this_fill_done = 0;
        end
    end
end

%end  % Lambda Loop end
fclose(fid);
fclose(fidq);
fprintf('\n');
fprintf('Done susceptibility calculation. Total Time = %.2f s\n',sum(vtm));
diary off;
