function [Am,Bm,Qm,Zm,Pm,Am1,Am0]=s_pre_t(bc,mc,R,z01,z02,Qt,Zt)
Ht=Zt-[z01,z02;z01,z02];
A1=bc*ones(2)+mc*Ht;
A1=A1.*Ht;
Am=R*((A1(1,1)+A1(1,2))/2)+(1-R)*((A1(2,1)+A1(2,2))/2);
B1=bc*ones(2)+2*mc*Ht;
Bm=R*((B1(1,1)+B1(1,2))/2)+(1-R)*((B1(2,1)+B1(2,2))/2);
Qm=R*((Qt(1,1)+Qt(1,2))/2)+(1-R)*((Qt(2,1)+Qt(2,2))/2);
Zm=R*((Zt(1,1)+Zt(1,2))/2)+(1-R)*((Zt(2,1)+Zt(2,2))/2);
Am1=(Zm-z02)*(bc+(mc*(Zm-z02)));
Am0=(Zm-z01)*(bc+(mc*(Zm-z01)));
P1=sqrt(mc^2+1)*ones(2);P1=bc*ones(2)+2*P1.*Ht;
Pm=R*((P1(1,1)+P1(1,2))/2)+(1-R)*((P1(2,1)+P1(2,2))/2);
