function [Q1,Z1] = s_Pre1 ( na , lc , nc , bc , mc , nl , lt , Rt , nt , n2 , P1 , Pn , DT , Q , Z , TZ0 , hj , Bsl )
G = 9.81;
Emax = 0.1;
R = 0.73;
N = na;
Q1 = zeros ( 1 , N );
Z1 = zeros ( 1 , N );
Qt = zeros ( 2 );
Zt = zeros ( 2 );
ZQ = zeros ( 2 * N , 1 );
for I = 1 : N
    ZQ ( 2 * I - 1 ) = Z ( I );                            % 未知量向量
    ZQ ( 2 * I ) = Q ( I );
end
ZQ1 = zeros ( 2 * N , 1 );
ZQ2 = ZQ;
Aeq = zeros ( 8 * N - 6 , 3 );
Beq = zeros ( 2 * N , 1 );
Aeq ( 1 , 1 ) = 1;
Aeq ( 1 , 2 ) = 2;
Aeq ( 1 , 3 ) = 1;
Aeq ( 8 * N - 6 , 1 ) = 2 * N;
Aeq ( 8 * N - 6 , 2 ) = 2 * N;
Aeq ( 8 * N - 6 , 3 ) = 1;
Beq ( 1 ) = P1;
Beq ( 2 * N ) = Pn;
temp = 0;
while norm ( ZQ1 - ZQ2 ) > Emax||temp <= 1                              %当ZQ1与ZQ2足够相似时停止计算
    temp = temp + 1;
    ZQ1 = ZQ2;
    for I = 1 : N - 1                             %计算方格编号
        if I < nl || I > nl + nt         %渠道→管道→渠道                           
            DS = lc / ( nl - 1 );
            Qt ( 1 , 1 ) = ZQ1 ( 2 * I );
            Qt ( 1 , 2 ) = ZQ1 ( 2 * I + 2 );
            Qt ( 2 , 1 ) = Q ( I );
            Qt ( 2 , 2 ) = Q ( I + 1 );
            Zt ( 1 , 1 ) = ZQ1 ( 2 * I - 1 );
            Zt ( 1 , 2 ) = ZQ1 ( 2 * I + 1 );
            Zt ( 2 , 1 ) = Z ( I );
            Zt ( 2 , 2 ) = Z ( I + 1 );
            z01 = TZ0 ( I );
            z02 = TZ0 ( I + 1 );
            [Am,Bm,Qm,~,Pm,Am1,Am0]=feval ( @s_pre_t ,bc , mc , R , z01 , z02 , Qt , Zt );
            c1=2*R*(DT/DS)/Bm;                                         %R：系数θ
            e1=Z(I)+Z(I+1)+((1-R)/R)*c1*(Q(I)-Q(I+1));
            a2=2*R*(DT/DS)*((Qm/Am)^2*Bm-G*Am);
            c2=1-4*R*(DT/DS)*(Qm/Am);
            d2=1+4*R*(DT/DS)*(Qm/Am);
            e2=((1-R)/R)*a2*(Z(I+1)-Z(I))+(1-4*(1-R)*(DT/DS)*(Qm/Am))*Q(I+1)+(1+4*(1-R)*(DT/DS)*(Qm/Am))*Q(I);
            e2=e2+2*DT*(Qm/Am)^2*(Am1-Am0)/DS-2*DT*(G*nc^2*Qm^2*Pm^(4/3))/(Am^(7/3));
            Aeq(8*I-6,1)=2*I;Aeq(8*I-6,2)=2*I-1;Aeq(8*I-6,3)=1;
            Aeq(8*I-5,1)=2*I;Aeq(8*I-5,2)=2*I;  Aeq(8*I-5,3)=-1*c1;
            Aeq(8*I-4,1)=2*I;Aeq(8*I-4,2)=2*I+1;Aeq(8*I-4,3)=1;
            Aeq(8*I-3,1)=2*I;Aeq(8*I-3,2)=2*I+2;Aeq(8*I-3,3)=c1;
            Beq(2*I)=e1;
            Aeq(8*I-2,1)=2*I+1;Aeq(8*I-2,2)=2*I-1;Aeq(8*I-2,3)=a2;
            Aeq(8*I-1,1)=2*I+1;Aeq(8*I-1,2)=2*I;  Aeq(8*I-1,3)=c2;
            Aeq(8*I,1)=2*I+1;  Aeq(8*I,2)=2*I+1;  Aeq(8*I,3)=-1*a2;
            Aeq(8*I+1,1)=2*I+1;Aeq(8*I+1,2)=2*I+2;Aeq(8*I+1,3)=d2;
            Beq(2*I+1)=e2;
        elseif I == nl || I == nt + nl
            Aeq ( 8 * I - 6 , 1 : 3 ) = [ 2 * I , 2 * I - 1 , 1 ];
            Aeq ( 8 * I - 5 , 1 : 3 ) = [ 2 * I , 2 * I  , 0 ];
            Aeq ( 8 * I - 4 , 1 : 3 ) = [ 2 * I , 2 * I + 1 , -1 ];
            Aeq ( 8 * I - 3 , 1 : 3 ) = [ 2 * I , 2 * I + 2 , 0 ];
            Beq ( 2*I ) = hj;                 % 节点处水位降落边界条件
            Aeq ( 8 * I - 2 , 1 : 3 ) = [ 2 * I + 1 , 2 * I - 1 , 0 ];
            Aeq ( 8 * I - 1 , 1 : 3 ) = [ 2 * I + 1 , 2 * I  , 1 ];
            Aeq ( 8 * I , 1 : 3 ) = [ 2 * I + 1 , 2 * I + 1 , 0 ];
            Aeq ( 8 * I + 1 , 1 : 3 ) = [ 2 * I + 1 , 2 * I + 2 , -1 ];
            Beq ( 2*I+1 ) = 0;          % 节点处流量差值边界条件
        else
            DS = lt / ( nt - 1 );
            Zt ( 1 , 1 ) = ZQ1 ( 2 * I - 1 );
            Zt ( 1 , 2 ) = ZQ1 ( 2 * I + 1 );
            Zt ( 2 , 1 ) = Z ( I );
            Zt ( 2 , 2 ) = Z ( I + 1 );
            Qt ( 1 , 1 ) = ZQ1 ( 2 * I );
            Qt ( 1 , 2 ) = ZQ1 ( 2 * I + 2 );
            Qt ( 2 , 1 ) = Q ( I );
            Qt ( 2 , 2 ) = Q ( I + 1 );
            [Am,Bm,Qm,Pm,Am1,Am0]=feval ( @s_pret , R , Rt , Qt , Zt , Bsl );
            c1=2*R*(DT/DS)/Bm;                                         %R：系数θ
            e1=Z(I)+Z(I+1)+((1-R)/R)*c1*(Q(I)-Q(I+1));
            a2=2*R*(DT/DS)*((Qm/Am)^2*Bm-G*Am);
            c2=1-4*R*(DT/DS)*(Qm/Am);
            d2=1+4*R*(DT/DS)*(Qm/Am);
            e2=((1-R)/R)*a2*(Z(I+1)-Z(I))+(1-4*(1-R)*(DT/DS)*(Qm/Am))*Q(I+1)+(1+4*(1-R)*(DT/DS)*(Qm/Am))*Q(I);
            e2=e2+2*DT*(Qm/Am)^2*(Am1-Am0)/DS-2*DT*(G*n2^2*Qm^2*Pm^(4/3))/(Am^(7/3));
            Aeq(8*I-6,1)=2*I;Aeq(8*I-6,2)=2*I-1;Aeq(8*I-6,3)=1;
            Aeq(8*I-5,1)=2*I;Aeq(8*I-5,2)=2*I;  Aeq(8*I-5,3)=-1*c1;
            Aeq(8*I-4,1)=2*I;Aeq(8*I-4,2)=2*I+1;Aeq(8*I-4,3)=1;
            Aeq(8*I-3,1)=2*I;Aeq(8*I-3,2)=2*I+2;Aeq(8*I-3,3)=c1;
            Beq(2*I)=e1;
            Aeq(8*I-2,1)=2*I+1;Aeq(8*I-2,2)=2*I-1;Aeq(8*I-2,3)=a2;
            Aeq(8*I-1,1)=2*I+1;Aeq(8*I-1,2)=2*I;  Aeq(8*I-1,3)=c2;
            Aeq(8*I,1)=2*I+1;  Aeq(8*I,2)=2*I+1;  Aeq(8*I,3)=-1*a2;
            Aeq(8*I+1,1)=2*I+1;Aeq(8*I+1,2)=2*I+2;Aeq(8*I+1,3)=d2;
            Beq(2*I+1)=e2;
        end
    end
    SP=full(spconvert(Aeq));         % 将系数记录矩阵Aeq转换为稀疏存储矩阵，再转换为全矩阵
    
    Deq=full(spdiags((max(abs(SP),[],2).^(-1)),0,2*N,2*N));
    % max(abs(SP),[],2):SP绝对值矩阵每一行的最大值；
    % spdiags((max(abs(SP),[],2).^(-1)),0,2*N,2*N):建立带状存储矩阵;
    Aeq1=Deq*SP;
    Beq1=Deq*Beq;
    % Dep的作用是防止病态方程的出现
    ZQ2=Aeq1\Beq1;
    % ZQ2=SP\Beq;
end
for I=1:N
    Z1(I)=ZQ2(2*I-1);
    Q1(I)=ZQ2(2*I);
end