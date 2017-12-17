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
    ZQ ( 2 * I - 1 ) = Z ( I );                            % δ֪������
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
while norm ( ZQ1 - ZQ2 ) > Emax||temp <= 1                              %��ZQ1��ZQ2�㹻����ʱֹͣ����
    temp = temp + 1;
    ZQ1 = ZQ2;
    for I = 1 : N - 1                             %���㷽����
        if I < nl || I > nl + nt         %�������ܵ�������                           
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
            c1=2*R*(DT/DS)/Bm;                                         %R��ϵ����
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
            Beq ( 2*I ) = hj;                 % �ڵ㴦ˮλ����߽�����
            Aeq ( 8 * I - 2 , 1 : 3 ) = [ 2 * I + 1 , 2 * I - 1 , 0 ];
            Aeq ( 8 * I - 1 , 1 : 3 ) = [ 2 * I + 1 , 2 * I  , 1 ];
            Aeq ( 8 * I , 1 : 3 ) = [ 2 * I + 1 , 2 * I + 1 , 0 ];
            Aeq ( 8 * I + 1 , 1 : 3 ) = [ 2 * I + 1 , 2 * I + 2 , -1 ];
            Beq ( 2*I+1 ) = 0;          % �ڵ㴦������ֵ�߽�����
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
            c1=2*R*(DT/DS)/Bm;                                         %R��ϵ����
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
    SP=full(spconvert(Aeq));         % ��ϵ����¼����Aeqת��Ϊϡ��洢������ת��Ϊȫ����
    
    Deq=full(spdiags((max(abs(SP),[],2).^(-1)),0,2*N,2*N));
    % max(abs(SP),[],2):SP����ֵ����ÿһ�е����ֵ��
    % spdiags((max(abs(SP),[],2).^(-1)),0,2*N,2*N):������״�洢����;
    Aeq1=Deq*SP;
    Beq1=Deq*Beq;
    % Dep�������Ƿ�ֹ��̬���̵ĳ���
    ZQ2=Aeq1\Beq1;
    % ZQ2=SP\Beq;
end
for I=1:N
    Z1(I)=ZQ2(2*I-1);
    Q1(I)=ZQ2(2*I);
end