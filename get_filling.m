function out = get_filling(Z,X,ek,WN,nk,beta,mu)

nktot = 4*nk*nk;
numwi = length(WN)/2;

fill = 0;       %for two spins
fill = fill + 2*sum(sum(fermi(ek-mu,beta)))/nktot;

green(1:nktot,1) = 0;
den(1:nktot,1) = 0;
Z = reshape(Z,nktot,[]);
X = reshape(X,nktot,[]);
ek = reshape(ek-mu,nktot,[]);
for nn = (numwi+1):(numwi*2)
    wn = WN(nn);
    den = Z(:,nn).^2 + (ek+X(:,nn)).^2;
    green = (ek+X(:,nn))./den;
    fill = fill - 4*sum(green)/(nktot*beta);       %4=2(from +/-wn) * 2(from trace<-->spin)

    den =  wn^2 + ek.^2;
    green = ek./den;
    fill = fill + 4*sum(green)/(nktot*beta);
end

out = fill;
