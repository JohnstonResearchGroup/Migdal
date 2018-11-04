function [Zut Xut P1ut P2ut mu] = interpSelfE(WNin,NUin,varargin)
%INTERPSELFE Interpolate the self-energy on different Matsubara frequency
%grid WN = pi*(2*(-numwi:numwi-1)+1)/beta.

flySelfE = varargin{1};
if flySelfE == 1
    Z = varargin{2};
    X = varargin{3};
    P1 = varargin{4};
    P2 = varargin{5};
    WN = varargin{6};
    NU = varargin{7};
    mu = varargin{8};
else
    dataFileName = varargin{2};
    load(dataFileName,'Z','X','P1','P2','WN','NU','mu','-mat')
end
[nk1 nk2 nw] = size(Z);
nwin = numel(WNin);
nwbin = numel(NUin); 
if mod(nw,2)~=0 || mod(nwin,2)~=0
    error('Number of Matsubara frequencies must be even!')
end
if nw~=numel(WN)
    error('Wrong grid for self-energy data!')
end

%interpolate electron self-energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = permute(Z,[3 1 2]);          %now Z becomes [nw nk1 nk2]
X = permute(X,[3 1 2]);
Zut(1:nwin,1:nk1,1:nk2) = 0;
Xut(1:nwin,1:nk1,1:nk2) = 0;

%1D interpolation with spline-extrapolated values
Zut = interp1(WN,Z,WNin,'spline','extrap');
Xut = interp1(WN,X,WNin,'spline','extrap');

%Replace spline-extrapolated values by 'nearest' values
idxL = find(WNin<WN(1));
idxR = find(WNin>WN(end));
if ~isempty(idxL)
    idxL = idxL(end);
    Zut(1:idxL,:,:) = Zut(repmat(idxL+1,1,idxL),:,:);
    Xut(1:idxL,:,:) = Xut(repmat(idxL+1,1,idxL),:,:);
end
if ~isempty(idxR)
    idxR = idxR(1);
    Zut(idxR:nwin,:,:) = Zut(repmat(idxR-1,1,nwin-idxR+1),:,:);
    Xut(idxR:nwin,:,:) = Xut(repmat(idxR-1,1,nwin-idxR+1),:,:);
end

Zut = permute(Zut,[2 3 1]);      %now Zut becomes [nk1 nk2 nwin]
Xut = permute(Xut,[2 3 1]);

%interpolate phonon self-energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nwbin == 0
    P1ut = [];      %now P1ut becomes [nk1 nk2 nwbin]
    P2ut = [];
else
    P1 = permute(P1,[3 1 2]);          %now P1 becomes [nw+1 nk1 nk2]
    P2 = permute(P2,[3 1 2]);
    P1ut(1:nwbin,1:nk1,1:nk2) = 0;
    P2ut(1:nwbin,1:nk1,1:nk2) = 0;

    %1D interpolation with spline-extrapolated values
    P1ut = interp1(NU,P1,NUin,'spline','extrap');
    P2ut = interp1(NU,P2,NUin,'spline','extrap');

    %Replace spline-extrapolated values by 'nearest' values
    idxL = find(NUin<NU(1));
    idxR = find(NUin>NU(end));
    if ~isempty(idxL)
        idxL = idxL(end);
        P1ut(1:idxL,:,:) = P1ut(repmat(idxL+1,1,idxL),:,:);
        P2ut(1:idxL,:,:) = P2ut(repmat(idxL+1,1,idxL),:,:);
    end
    if ~isempty(idxR)
        idxR = idxR(1);
        P1ut(idxR:nwbin,:,:) = P1ut(repmat(idxR-1,1,nwbin-idxR+1),:,:);
        P2ut(idxR:nwbin,:,:) = P2ut(repmat(idxR-1,1,nwbin-idxR+1),:,:);
    end

    P1ut = permute(P1ut,[2 3 1]);      %now P1ut becomes [nk1 nk2 nwbin]
    P2ut = permute(P2ut,[2 3 1]);
end
