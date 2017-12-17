clc;
close all;
num1 = xlsread ( 'data' , 1 );
lp = num1 ( 1 );
i = num1 ( 2 );
n1 = num1 ( 3 );
b = num1 ( 4 );
m = num1 ( 5 );
nl = num1 ( 6 );
zu = lp * i;
zd = -lp * i;
hu = num1 ( 7 );
num2 = xlsread ( 'data' , 2 );
lt = num2 ( 1 );
n2 = num2 ( 2 );
R = num2 ( 3 );
nt = num2 ( 4 );
a = num2 ( 5 );
num3 = xlsread ( 'data' , 3 );
Qini = num3 ( 1 );
Ts = num3 ( 2 );
Te = num3 ( 3 );
Tf = num3 ( 4 );
Dt = num3 ( 5 );
Ts = Ts * 60 / Dt;
Te = Te * 60 / Dt;
Tf = Tf * 60 / Dt;
num4 = xlsread ( 'data' , 4 );
mq = num4 ( 2 );
bq = num4 ( 3 );
he = num4 ( 5 );
num5 = xlsread ( 'data' , 5 );
trigger = num5 ( 1 );
A = pi * R^2 / 4;
Bsl = 9.8 * A / a^2;
DSL = lp / nl;
DST = lt / nt;
Qe = zeros ( Tf , 1 );
Qa = zeros ( Tf , 1 );
He = zeros ( Tf , 1 );
Q = zeros ( Tf , 1 );
H1 = zeros ( Tf , 1 );
H2 = zeros ( Tf , 1 );
H3 = zeros ( Tf , 1 );
H4 = zeros ( Tf , 1 );
Dh = zeros ( Tf , 1 );
hf = zeros ( Tf , 1 );
hi = zeros ( Tf , 1 );
z = zeros ( 1 , 2 * nl + nt );
z ( 1 : nl + 1 ) = zu : -zu / nl : 0;
z ( nl + 1 : nl + nt + 1 ) = 0;
z ( nl + nt + 1 : 2 * nl + nt + 1 ) = 0 : zd / nl : zd;
Q(1:Ts)=Qini;
for step = Ts+1:Tf
    %Q(step)=Qini + sin(pi/100*(step-Ts));
    Q(step)=Qini + 5;
end
switch trigger
    case 1
        for step = 1 : Tf
            if step == 1               
                Hp2 = sf_h_d ( DSL , i , n1 , b , m , nl , Qini , hu );
                Ht = st_h_d ( DST , n2 , R , nt , Qini , Hp2 ( 1 ) );
                Hp1 = sf_h_d ( DSL , i , n1 , b , m , nl + 1 , Qini , Ht ( 1 ) );
                Qall ( 1 : nl + 1 ) = Qini;
                Qall ( nl + 1 : nl + nt + 1 ) = Qini;
                Qall ( nl + nt + 1 : nt + 2 * nl + 1 ) = Qini;
                Zall ( 1 : nl + 1 ) = z ( 1 : nl + 1 ) + Hp1;
                Zall ( nl + 2 : nl + nt + 1 ) = z ( nl + 2 : nl + nt + 1 ) + Ht';
                Zall ( nl + nt + 2 : nt + 2 * nl + 1 ) = Hp2 + z ( nl + nt + 2 : 2 * nl + nt + 1 );
%                 Q(1:Ts)=Qini;
%                 Q(Ts+1:Tf)=Qini - 2;
                Qn = mq * bq * hu * ( 2 * 9.8 * hu )^0.5;
                Qe ( step ) = Qall ( nl + nt + 1 );
                Qa ( step ) = Qall ( nl + 1 );
                H1 ( step ) = Zall ( nl + 1);
                H2 ( step ) = Zall ( nl + nt );
                Ha = Zall ( nl + 1 );
                He = Zall ( nl + nt + 1 );
                na = 2 * nl + nt;
%                 figure ( 4 )
%                 plot ( [0 :lp/nl: lp , lp + lt/nt : lt/nt : lp + lt , lp + lt + lp/nl : lp/nl : 2 * lp + lt] , z ( 1 : na + 1 ) ,  [0 :lp/nl: lp , lp + lt/nt : lt/nt : lp + lt , lp + lt + lp/nl : lp/nl : 2 * lp + lt] , Zall ( 1 : na + 1 ) );
%                 legend('底坡','FontSize',14,'测压管水头线','FontSize',14)
%                 ylabel ( '水位(m)'),'FontSize',14; xlabel('水平距离（m）','FontSize',14) ;
%                 axis( [ 0 2*lp+lt+40 zd-2 Zall(1)+2 ] );
            else
                Qn ( step ) = he * mq * bq * ( Zall ( na ) - zd ) * ( 2 * 9.8 * ( Zall ( na ) - zd ) )^0.5;
                [ Qall , Zall ] = s_Pre1 ( na , lp , n1 , b , m , nl , lt , R , nt , n2 , Q( step ) , Qn ( step ) , Dt , Qall , Zall , z , 0 , Bsl );
                Qe ( step ) = Qall ( nl + nt );
                Qa ( step ) = Qall ( nl + 1 );
                H1 ( step ) = Zall ( nl + 1 );
                H2 ( step ) = Zall ( nl + nt );
                Zp1 = Zall ( 1 : nl );
            end
        end
        for x = 2 : Tf
         hf ( x ) = ( ( n2^2 ) * ( 4 * Qa(x) / pi / R^2 )^2 ) / ( ( R / 4 )^(4/3) ) * lt;
         H3 ( x ) = H1 ( x ) - 1.000*hf ( x );
         hi(x) = (Qa(x)-Qa(x-1))*lt/Dt/3.14/9.8;
         H4(x) = H3(x) - hi(x);
        end       
        x = 20 : Tf;      
        figure ( 1 )
        plot ( x / 60 * Dt , H1 ( 20 : Tf ) , x / 60 * Dt , H2 ( 20 : Tf ) )
        legend({'倒虹吸前水位','倒虹吸后水位'},'FontSize',14);
        ylabel ( '水位(米)','FontSize',14) ; xlabel('时间（分钟)','FontSize',14) 
        figure ( 2 )
        plot ( x / 60 * Dt , Qa ( 20 : Tf ) )
        legend({'倒虹吸内流量'},'FontSize',14);
        ylabel ( '流量(m3/s)','FontSize',14) ; xlabel('时间（分钟)','FontSize',14) 
        
%         figure ( 3 )
        figure(5)        
        plot ( x / 60 * Dt , H3 ( 20 : Tf ), '+', x / 60 * Dt , H2 ( 20 : Tf ) )
        ylabel ( '水位(m)') ; xlabel('时间（分钟)') ;
        legend({'H1-hf','H2'});        
        figure(6)
        plot ( x / 60 * Dt , H4 ( 20 : Tf ), '+', x / 60 * Dt , H2 ( 20 : Tf ) )
        ylabel ( '水位(m)','FontSize',20) ; xlabel('时间（分钟)','FontSize',20) ;
        legend({'H1-hf-hi','H2'},'FontSize',24);

%         figure ( 3 )
%         plot ( x / 60 * Dt , Qe ( 2 : Tf ) )
%         ylabel ( '流量(m3/s)') ; xlabel('时间（h)') ;
%         legend('管道内流量')
    case 2
        lp = 8.8;
        wch = 0.51;
        hch = 0.148;
        i =0;
        n =0.012;
        Bsl = 0.01;
        lt = 0.11;
        nl = lp/lt + 1;
        DT = 0.1;
        Ts = 6.6/DT+1;
        Z = zeros ( nl , 1 );
        Hp = zeros ( nl , 1 );
        Qp = zeros ( nl , 1 );
        Z ( 1 : nl ) = 0;
        H ( 1 ) = 0.1280;
        H ( 2 ) = 0.1290;
        H ( 3 ) = 0.1300;
        H ( 4 ) = 0.1330;
        H ( 5 ) = 0.1360;
        H ( 6 ) = 0.1400;
        H ( 7 ) = 0.1420;
        H ( 8 ) = 0.1440;
        H ( 9 ) = 0.1460;
        H ( 10 ) = 0.1470;
        H ( 11 ) = 0.1500;
        H ( 12 ) = 0.1516;
        H ( 13 ) = 0.1540;
        H ( 14 ) = 0.1580;
        H ( 15 ) = 0.1610;
        H ( 16 ) = 0.1640;
        H ( 17 ) = 0.1688;
        H ( 18 ) = 0.1728;
        H ( 19 ) = 0.1752;
        H ( 20 ) = 0.1778;
        H ( 21 ) = 0.1800;
        H ( 22 ) = 0.1830;
        H ( 23 ) = 0.1870;
        H ( 24 ) = 0.1900;
        H ( 25 ) = 0.1920;
        H ( 26 ) = 0.1920;
        H ( 27 ) = 0.1916;
        H ( 28 ) = 0.1896;
        H ( 29 ) = 0.1880;
        H ( 30 ) = 0.1880;
        H ( 31 ) = 0.1888;
        H ( 32 ) = 0.1908;        
        H ( 33 ) = 0.1948;
        H ( 34 ) = 0.1972;
        xi = 0.1:0.1:6.7;
        H=interp1(0.1:0.2:6.7, H, xi, 'linear');
        Hch(1:Ts+1)=hch;
        plot ( xi , H , xi , Hch )
        legend('上游水位','管顶')
        ylabel ( 'H(m)') ; xlabel('t（s）') ;
        for step = 1 : Ts
            if step == 1
                Hp ( 1 : nl, 1 ) = 0.128;
                Qp ( 1 : nl, 1 ) = 0;
            else
                [ Qp, Hp ] = s_Pretonlytest ( nl , lp , wch , hch , nl , n , H ( step ) , 0.1280 , DT , Qp , Hp , Bsl );
                if step == 16
                    H1 = Hp;
                else
                    if step == 31
                        H2 = Hp;
                    else
                        if step == 46
                            H3 = Hp;
                        else
                            if step == 61
                                H4 = Hp;
                            end
                        end
                    end
                end
            end
        end   
        H11 = [ 0.144852 0.143865 0.142877 0.141749 0.140621 0.13921 0.137659 0.136107 0.134556...
            0.133145 0.131735 0.130748 0.129901 0.129196 0.128773 0.128491 0.128209 0.128068...
            0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128...
            0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128...
            0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128...
            0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128...
            0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128 0.128];
        H22 = [0.163893	0.16347	0.163047 0.162623 0.1622 0.161777 0.161213 0.16079 0.160367...
            0.159944 0.159803 0.159238 0.158392	0.157546 0.157264 0.157828 0.157405	0.146403...
            0.1433 0.141044	0.13921	0.137518 0.135966 0.134556 0.133286	0.132158 0.131171...
            0.130183 0.129619 0.129196 0.128773	0.128491 0.128209 0.128068 0.128 0.128...
            0.128 0.128 0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128...
            0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128...
            0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128...
            0.128 0.128	0.128 0.128	0.128 0.128];
        H33 = [0.188434	0.188011 0.187588 0.187024 0.186601	0.186178 0.185755 0.185331 0.184767...
            0.184344 0.18378 0.183357 0.182793 0.182228	0.181805 0.181241 0.180818 0.180254...
            0.17969	0.179267 0.178702 0.178279 0.177715	0.177151 0.176728 0.176305 0.175882...
            0.175458 0.174894 0.17433 0.173766 0.173343	0.173343 0.173202 0.172638 0.171368...
            0.169817 0.169676 0.17165 0.17165 0.143865 0.136389	0.132581 0.130748 0.12976...
            0.129055 0.128773 0.128491 0.12835 0.128209	0.128068 0.128 0.128 0.128 0.128...
            0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128...
            0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128 0.128	0.128];
        H44 = [0.188575	0.188434 0.188152 0.18787 0.187729 0.187447	0.187165 0.186883 0.186601...
            0.186319 0.186178 0.185896 0.185614	0.185331 0.18519 0.184908 0.184626 0.184485...
            0.184203 0.184062 0.183921 0.183639	0.183498 0.183216 0.183075 0.182934	0.182793...
            0.182652 0.182511 0.18237 0.182228 0.182087	0.181946 0.181946 0.181805 0.181805...
            0.181664 0.181664 0.181523 0.181523	0.181523 0.181523 0.181382 0.181382	0.181382...
            0.181382 0.181382 0.181241 0.181241	0.181241 0.181241 0.181241 0.181241	0.180959...
            0.180818 0.180818 0.181241 0.181664	0.181241 0.179972 0.178843 0.179831	0.182934...
            0.18237	0.143583 0.134274 0.130042 0.128632	0.128068 0.128 0.128 0.128 0.128 0.128...
            0.128 0.128	0.128 0.128	0.128 0.128	0.128];
        x1 = 0 : lt :lp;
        x2 = 0 : 0.125 : 10;
        plot ( x2 , H1 , x2 , H2 , x2 , H3 , x2 , H4 , x2 , H11 , x2 , H22 , x2 , H33 , x2 , H44 )
        axis ( [ 0 , 10 , 0.1 , 0.2 ] )
        ylabel ( 'H(m)') ; xlabel('x（m）') ;
    case 3
        H1 = zeros ( Tf , 1 );
        H2 = zeros ( Tf , 1 );
        H3 = zeros ( Tf , 1 );
        H4 = zeros ( Tf , 1 );
        Dh = zeros ( Tf , 1 );
        hf = zeros ( Tf , 1 );
        z = zeros ( 1 , 2 * nl + nt );
        A = pi * R^2 / 4;
        Bsl = 9.8 * A / a^2;
        DSL = lp / nl;
        DST = lt / nt;
        H = zeros ( 1 , Tf );
        for step = 1 : Tf
            if step == 1
                z ( 1 : nl + 1 ) = zu : -zu / nl : 0;
                z ( nl + 1 : nl + nt + 1 ) = 0;
                z ( nl + nt + 1 : 2 * nl + nt + 1 ) = 0 : zd / nl : zd;
                Hp2 = sf_h_d ( DSL , i , n1 , b , m , nl , Q ( step ) , hu );
                Ht = st_h_d ( DST , n2 , R , nt , Q ( step ) , Hp2 ( 1 ) );
                Hp1 = sf_h_d ( DSL , i , n1 , b , m , nl + 1 , Q ( step ) , Ht ( 1 ) );
                Qall ( 1 : nl + 1 ) = Q ( step );
                Qall ( nl + 1 : nl + nt + 1 ) = Q ( step );
                Qall ( nl + nt + 1 : nt + 2 * nl + 1 ) = Q ( step );
                Zall ( 1 : nl + 1 ) = z ( 1 : nl + 1 ) + Hp1;
                Zall ( nl + 2 : nl + nt + 1 ) = z ( nl + 2 : nl + nt + 1 ) + Ht';
                Zall ( nl + nt + 2 : nt + 2 * nl + 1 ) = Hp2 + z ( nl + nt + 2 : 2 * nl + nt + 1 );
                hp = Zall ( 1 );
                H ( step ) = hp;
                Qn = he * mq * bq * hu * ( 2 * 9.8 * hu )^0.5;
                Qe ( step ) = Qall ( nl + nt + 1 );
                H1 ( step ) = Zall ( nl );
                H2 ( step ) = Zall ( nl + nt );
                na = 2 * nl + nt + 1;
                figure ( 4 )
                plot ( [0 :lp/nl: lp , lp + lt/nt : lt/nt : lp + lt , lp + lt + lp/nl : lp/nl : 2 * lp + lt] , z ( 1 : na ) , [0 :lp/nl: lp , lp + lt/nt : lt/nt : lp + lt , lp + lt + lp/nl : lp/nl : 2 * lp + lt] , Zall ( 1 : na ) , '--' );
                legend('底坡','水面线')
                ylabel ( '水位(m)') ; xlabel('水平距离（m）') ;
                axis( [ 0 2*lp+lt+40 zd-2 Zall(1)+2 ] );
            else
                H ( step ) = hp + 0.1 * sin( ( step - 1 ) * pi / 60 );
                Qn ( step ) = he * mq * bq * ( Zall ( 2 * nl + nt + 1 ) - zd ) * ( 2 * 9.8 * ( Zall ( 2 * nl + nt + 1 ) - zd ) )^0.5;
                for SV_I = 1 : SV
                    [ Qall , Zall ] = s_Pre5 ( na , lp , n1 , b , m , nl , lt , R , nt , n2 , H ( step ) , Qn ( step ) , Dt , Qall , Zall , z , 0 , Bsl );
                end
                Qe ( step ) = Qall ( nl + nt + 1 );
                H1 ( step ) = Zall ( nl );
                H2 ( step ) = Zall ( nl + nt );
                Zp1 = Zall ( 1 : nl );
            end
        end
        for x = 2 : Tf
            hi ( x ) = lt / 9.8 / 3.14 * ( Qe ( x ) - Qe ( x - 1 ) ) / Dt / 60;
            hf ( x ) = ( ( n2^2 ) * ( 4 * Qe(x) / pi / R^2 )^2 ) / ( ( R / 4 )^(4/3) ) * lt;
            H3 ( x ) = H1 ( x ) - 1.000*hf ( x );
            H4 ( x ) = H3 ( x ) - 30 * hi ( x );
        end
        x = 2 : Tf;
        figure ( 1 )
        subplot(1,2,1)
        plot ( x / 60 * Dt , H1 ( 2 : Tf ) )
        legend('管道前水位');
        ylabel ( '水位(m)') ; xlabel('时间（h)')
        subplot(1,2,2)
        plot ( x / 60 * Dt , H2 ( 2 : Tf ) )
        legend('管道后水位');
        ylabel ( '水位(m)') ; xlabel('时间（h)')
        figure ( 2 )
        subplot(1,2,1)
        plot ( x / 60 * Dt , H3 ( 2 : Tf ), x / 60 , H2 ( 2 : Tf ) ,'--')
        ylabel ( '水位(m)') ; xlabel('时间（h)') ;
        legend('H1-hf','H2');
        subplot(1,2,2)
        plot ( x / 60 * Dt , H4 ( 2 : Tf ), x / 60 , H2 ( 2 : Tf ) ,'--')
        ylabel ( '水位(m)') ; xlabel('时间（h)') ;
        legend('H1-hf-k*hi','H2');
        figure ( 3 )
        plot ( x / 60 * Dt , Qe ( 2 : Tf ) )
        ylabel ( '流量(m3/s)') ; xlabel('时间（h)') ;
        legend('管道内流量')
end