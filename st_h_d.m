function Ht = st_h_d ( DST , n2 , R , T1 , Q , Hout )
Ht = zeros ( T1 , 1 );
Ht ( T1 ) = Hout;
l = DST;
v = 4 * Q / pi / R^2;
for kk = T1 - 1 : -1 : 1
    Ht ( kk ) = Ht ( kk + 1 ) + l * ( n2^2 * v^2 / ( R / 4 )^( 4 / 3 ) );
end
