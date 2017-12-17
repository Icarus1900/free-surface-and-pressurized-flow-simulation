function [Am,Bm,Qm,Pm,Am1,Am0] = s_pret ( R , Rt , Qt , Z , Bsl )
A = zeros ( 2 );
B = zeros ( 2 );
P = zeros ( 2 );
for i = 1 : 2
    for j = 1 : 2 
if Z ( i , j ) >= Rt
    A ( i , j ) = Rt^2 * pi / 4;
    B ( i , j ) = Bsl;
    P ( i , j ) = pi * Rt;
elseif Z ( i , j ) >= Rt / 2
    A ( i , j ) = ( ( 2 * pi - 2 * acos ( ( 2 * Z ( i , j ) - Rt ) / Rt ) ) / 2 / pi ) * Rt^2 * pi / 4 + ( Z ( i , j ) - Rt / 2 ) * (  ( Rt / 2 )^2 - ( Z ( i , j ) - Rt / 2 )^2 ) ^0.5;
    B ( i , j ) = 2 * ( ( Rt / 2 )^2 - ( Z ( i , j ) - Rt / 2 )^2 ) ^0.5;
    P ( i , j ) = ( ( 2 * pi - 2 * acos ( ( 2 * Z ( i , j ) - Rt ) / Rt ) ) / 2 / pi ) * pi * Rt;
else
    A ( i , j ) = ( ( 2 * acos ( ( Rt - 2 * Z ( i , j ) ) / Rt ) ) / 2 / pi ) * Rt^2 * pi / 4 - (  Rt / 2 - Z ( i , j ) ) * (  ( Rt / 2 )^2 - ( Z ( i , j ) - Rt / 2 )^2 ) ^0.5;   
    B ( i , j ) = 2 * ( ( Rt / 2 )^2 - (  Rt / 2  - Z ( i , j ) )^2 ) ^0.5;
    P ( i , j ) = ( ( 2 * acos ( ( Rt - 2 * Z ( i , j ) ) / Rt ) ) / 2 / pi ) * pi * Rt;
end
    end
end
Am=R*((A(1,1)+A(1,2))/2)+(1-R)*((A(2,1)+A(2,2))/2);
Bm=R*((B(1,1)+B(1,2))/2)+(1-R)*((B(2,1)+B(2,2))/2);
Qm=R*((Qt(1,1)+Qt(1,2))/2)+(1-R)*((Qt(2,1)+Qt(2,2))/2);
Am1=Am;
Am0=Am;
Pm=R*((P(1,1)+P(1,2))/2)+(1-R)*((P(2,1)+P(2,2))/2);