function H = sf_h_d ( DSL , i , n , b , m , T1 , Q , Hout )
N = n;
I = i;
B = b;
M = m;
D = 1;
G = 9.8;
E1 = 1e-8;
DS=DSL;
H=zeros(1,T1);R=zeros(1,T1);C=zeros(1,T1);U=zeros(1,T1);
I1=sqrt(1-I*I);
H(T1)=Hout;                        % 已知下游断面边界水深，向上游推求 %子渠段最后一个计算断面的参数
A=B*H(T1)+M*H(T1)*H(T1);             % A:断面面积
P=B+2*H(T1)*sqrt(1+M*M);           % P:湿周
R(T1)=A/P;                         % R:水力半径
C(T1)=exp(log(R(T1))/6)/N;          % C:谢才系数
U(T1)=Q/A;                         % U:流速
for T=(T1-1):-1:1
    H2=H(T+1);H1=H(T+1)-0.5;            % H1,H2:计算水深的初始边界
    A=B*H1+M*H1*H1;
    P=B+2*H1*sqrt(1+M*M);
    R1=A/P;
    C1=exp(log(R1)/6)/N;
    U1=Q/A;
    R3=(R1+R(T+1))/2;                   %取相邻两个计算断面水力参数的平均值计算沿程水头损失系数
    C3=(C1+C(T+1))/2;
    U3=(U1+U(T+1))/2;
    J=U3*U3/(C3*C3*R3);                   %Jf
    E2=(I-J)*DS;                          %断面比能差
    F6=I1*(H(T+1)-H1)+D*(U(T+1)^2-U1^2)/(2*G)-E2;     %可确定H1比实际的H(T)小，F6可能是正可能是负
    
    while abs(H1-H2)>E1                  %运用递代运算使H1趋近于断面实际水深
      H3=(H1+H2)/2;
      A=B*H3+M*H3*H3;
      P=B+2*H3*sqrt(1+M*M);
      R1=A/P;
      C1=exp(log(R1)/6)/N;
      U1=Q/A;
      R3=(R1+R(T+1))/2;
      C3=(C1+C(T+1))/2;
      U3=(U1+U(T+1))/2;
      J=U3*U3/(C3*C3*R3);
      E2=(I-J)*DS;
      F7=I1*(H(T+1)-H3)+D*(U(T+1)^2-U1^2)/(2*G)-E2;           % 与书中F(h)的符号相同
      if F6*F7<0                 %二分点换边  F6与F7反号，说明H3大于实际断面水深，向H减小的方向取二分点
        H2=H1;                         
      end
      H1=H3;
      F6=F7;
    end
    
H(T)=H3;
U(T)=U1;
C(T)=C1;
R(T)=R1;
end