function [S_mp,SY_mp]=LatticeSumMat_faster(M,k_perp,xi,m,mx,A,Qhx,Thx,acc,accx,BesselAccs)
tmpo=2*M+1; tmpt=2*M+2; mpo=M+1; mpt=M+2;


[Sm,SmY]=getlatticesums_faster(k_perp,xi,m,mx,A,Qhx,Thx,acc,accx,BesselAccs);
for p=-M-1:M 
    for m=-M-1:M
        S_mp(p+mpt,m+mpt)=Sm(m-p+tmpt);
        SY_mp(p+mpt,m+mpt)=SmY(m-p+tmpt);
    end
end
end
