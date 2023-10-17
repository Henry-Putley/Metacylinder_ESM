function [Sm,SmY] = getlatticesums_faster(k2,xi,m,mx,A,Qhx,Thx,acc,accx,BesselAccs)
% function for obtaining the lattice sums


Smx = 1./besselj(m+acc,k2.*xi).*( ... 
    -i.^m.*4./A.*sum( ... 
       (k2./Qhx).^accx.* ...
       BesselAccs ... % This is pre-calculated at each kB 
       ./(Qhx.^2-k2.^2).*exp(i.*(Thx.*mx)) ...
                    ,1)...
                                );
                            
 
nn = 1:acc(1);
 Smx(1) = Smx(1) +...
     -1./besselj(acc(1),k2.*xi).*(bessely(acc(1),k2.*xi) + 1/pi.*sum( ...
            factorial(acc(1)-nn)./factorial(nn-1).*(2./k2./xi).^(acc(1)-2.*nn+2) ) ) ;
         
SmY = [conj(fliplr(Smx(2:end))) Smx];
SmJ = [zeros(size(Smx(2:end))) -1 zeros(size(Smx(2:end)))];

Sm = SmJ + i.*SmY;

return
                            
                           
                            
                            
                            
                            
                            
                            










