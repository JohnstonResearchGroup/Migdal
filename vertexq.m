function vq = vertexq(KX,KY,varargin)
%VERTEXQ Generates the momentum dependence of the electron-phonon coupling
%The vertex is assumed to have the from |g(q)|^2 = g_0^2*f(q).
%   f(q) = vq. The normalization is defined using 
%   [1/N_q]\sum_{q_x,q_y} f(q) = 1 in discrete form
%   or [1/(2\pi)^2]\int_{q\in BZ} dq_x dq_y f(q) = 1 in continuous form.
% 
%   For cuprates, see Devereaux/Virosztek/Zawadowski PRB1995;
%   Bulut/Scalapino PRB1996; Sandvik/Scalapino/Bickers PRB2004.
%             O(3)
%             |
%   Cu--O(2)--Cu--O(2)
%   |         |
%   O(3)      O(3)      y
%   |         |         |
%   Cu--O(2)--Cu        |-->x
%   (0) apical oxygen mode: O(4)-Cu, where O(4) is located above each Cu.
%   (1) breathing mode: O(2)-Cu, O(3)-Cu (half breathing: only O(2) or O(3)).
%   (2) buckling mode: Cu-O(2)-Cu,  Cu-O(3)-Cu. A1g: O(2)/O(3) in phase;
%   B1g: O(2)/O(3) out of phase.
%   From Bulut/Scalapino PRB1996:
%   (0) |g(q)|^2 = g_0^2;
%   (1) |g(q)|^2 = g_0^2*(sin(qx/2)^2 + sin(qy/2)^2);
%   (2) |g(q)|^2 = g_0^2*(cos(qx/2)^2 + cos(qy/2)^2);
%   where g_0^2 = (2*)|g|^2/(2M\omega_ph) [factor 2 for case (1) and
%   (2) because of the presence of two oxygens].
%   funtype:
%   'br': vq = sin(qx/2)^2 + sin(qy/2)^2
%   'bk': vq = cos(qx/2)^2 + cos(qy/2)^2
%
%   We also include functionality for forward scatting. 
%   funtype: 'exp', 'gauss', 'kd', 'window' for forward scattering coupling
%   with a peak at q=0, and an optional width q0.
%

funtype = varargin{1};
if strcmp(funtype,'exp')
    q0 = varargin{2}/2;
    vq = 4*(pi^2)*exp(-sqrt(KX.*KX+KY.*KY)/q0)/(2*pi*q0*q0);
elseif strcmp(funtype,'gauss')
    q0 = varargin{2};
    vq = 4*(pi^2)*exp(-(KX.*KX + KY.*KY)/(2*q0*q0))/(2*pi*q0*q0);
elseif strcmp(funtype,'kd')
    nk = varargin{2};
    vq = ((KX.*KX+KY.*KY)==0) * ((2*nk)^2);
elseif strcmp(funtype,'window')
    nk = varargin{2};
    delq = varargin{3};
    vqmat = (KX.*KX+KY.*KY)<=(delq^2);
    Nq = numel(find(vqmat));
    vq = vqmat*((2*nk)^2)/Nq;
elseif strcmp(funtype,'br')
    nk = varargin{2};
    vq = sin(KX/2).^2 + sin(KY/2).^2;
elseif strcmp(funtype,'bk')
    nk = varargin{2};
    vq = cos(KX/2).^2 + cos(KY/2).^2;
else
    error([funtype, ' is not implemented!'])
end
