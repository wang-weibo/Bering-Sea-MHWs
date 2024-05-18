% Ekman Transport with finite depth
% 2023/3/16

function [Ue,Ve]=ekman_finite(lat,lon,u10,v10,dep,sic)

Az=0.01;
rho = 1025; 
h=dep;
Cd=1.25e-3;
% for i=1:m
%     for j=1:n
%         a=v10(i,j);
%         if(a>10)
%             Cd(i,j)=0.49+0.065*a;
%         end
%         if(a>=3 && a<=10)
%             Cd(i,j)=1.14;
%         end
%         if(a<=3)
%             Cd(i,j)=0.62+1.56/a;
%         end
%     end
% end
% Cd=Cd*10^-3;

f=coriolisf(lat);
a=sqrt(f/(2*Az));
[taux,tauy]=windstress(u10,v10,'ci',sic);

ddc1=(tauy./(a*rho*Az)).*(1./(cosh(2*a.*h)+cos(2*a.*h)));
sx1=(ddc1./a).*(0.5*(cosh(2*a.*h)+cos(2*a.*h))-cosh(a.*h).*cos(a.*h));
sy1=(ddc1./a).*sinh(a.*h).*sin(a.*h);

ddc2=(taux./(a*rho*Az)).*(1./(cosh(2*a.*h)+cos(2*a.*h)));
sx2=(ddc2./a).*(0.5*(cosh(2*a.*h)+cos(2*a.*h))-cosh(a.*h).*cos(a.*h));
sy2=(ddc2./a).*sinh(a.*h).*sin(a.*h);

Ue=sx1+sy2;
Ve=sy1-sx2;

end