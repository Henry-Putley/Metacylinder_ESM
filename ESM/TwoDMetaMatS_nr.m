function P=TwoDMetaMatS_nr(M,delta,BesselJs,BesselHs,dBesselJs,dBesselHs,S_mp,nr,BesselJsin,dBesselJsin)
% Function describing an infinite array of metacylinders
P=zeros(4*M+4);
tmpt=2*M+2;
mpo=M+1;

for p=-M-1:1:M
    for m=-M-1:1:M
        S_mp(p+mpo+1,m+mpo+1)=(-1)^(m-p)*S_mp(p+mpo+1,m+mpo+1)/BesselHs(m+tmpt+1);
    end
end


P(1:tmpt,1:tmpt)=BesselHs(mpo+1:mpo+tmpt)./BesselHs(mpo+1:mpo+tmpt).*eye(tmpt)+diag(BesselJs(mpo+1:mpo+tmpt))*S_mp; 
P(tmpt+1:end,1:tmpt)=dBesselHs(mpo+1:mpo+tmpt)./BesselHs(mpo+1:mpo+tmpt).*eye(tmpt)+diag(dBesselJs(mpo+1:mpo+tmpt))*S_mp;

for p=-M-1:1:M  
    % Top half
    P(p+mpo+1,tmpt+1:tmpt+mpo)=-(1/2)*exp(-1i*delta*p)*(1i^p).*((-1).^(0:M)).*(BesselJsin(p+tmpt+1:p+tmpt+1+M)+BesselJsin(p+tmpt+1:-1:p+tmpt+1-M));
    P(p+mpo+1,tmpt+mpo+1:end)=-(1/2)*exp(-1i*delta*p)*(1i^p).*((-1).^p).*(BesselJsin(p+tmpt+1:p+tmpt+1+M)+BesselJsin(p+tmpt+1:-1:p+tmpt+1-M));
    % Bottom half
    P(p+tmpt+mpo+1,tmpt+1:tmpt+mpo)=-(1/nr)*(1/2)*exp(-1i*delta*p)*(1i^p).*((-1).^(0:M)).*(dBesselJsin(p+tmpt+1:p+tmpt+1+M)+dBesselJsin(p+tmpt+1:-1:p+tmpt+1-M));
    P(p+tmpt+mpo+1,tmpt+mpo+1:end)=-(1/nr)*(1/2)*exp(-1i*delta*p)*(1i^p).*((-1).^p).*(dBesselJsin(p+tmpt+1:p+tmpt+1+M)+dBesselJsin(p+tmpt+1:-1:p+tmpt+1-M));
end


end