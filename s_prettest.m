function [Am,Bm,Qm,Pm,Am1,Am0,Vm] = s_prettest ( R , wch , hch , Qt , Z , Bsl )
A = zeros ( 2 );
B = zeros ( 2 );
P = zeros ( 2 );
for i = 1:2
    for j = 1:2
        if Z ( i , j ) <= hch
            B ( i , j ) = wch;
            A ( i , j ) = Z ( i , j ) * B ( i , j );
            P ( i , j ) = B ( i , j ) + Z ( i , j ) * 2;        
        else
            B ( i , j ) = Bsl;
            A ( i , j ) = hch * wch;
            P ( i , j ) = hch * 2 + wch * 2;
        end
    end
end
Am=R*((A(1,1)+A(1,2))/2)+(1-R)*((A(2,1)+A(2,2))/2);
Bm=R*((B(1,1)+B(1,2))/2)+(1-R)*((B(2,1)+B(2,2))/2);
Qm=R*((Qt(1,1)+Qt(1,2))/2)+(1-R)*((Qt(2,1)+Qt(2,2))/2);
Vt = Qt./A;
Vm=R*((Vt(1,1)+Vt(1,2))/2)+(1-R)*((Vt(2,1)+Vt(2,2))/2);
Zm=R*((Z(1,1)+Z(1,2))/2)+(1-R)*((Z(2,1)+Z(2,2))/2);
if Zm <= hch
    Am0 = Zm * wch;
else
    Am0 = hch * wch;
end
Am1 = Am0;
Pm=R*((P(1,1)+P(1,2))/2)+(1-R)*((P(2,1)+P(2,2))/2);