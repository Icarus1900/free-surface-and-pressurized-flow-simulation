function Z=Ini_Z(CANAL,N1,Nout,Qout,TZ0,Zouth,Q0,Hd0, I )
% Calculate the water surface profile of a canal pool with constant flow as the initial condition of the simulation
% Z：计算结果（sum(N1)维数组）
% CANAL：渠段参数矩阵 例如CANAL(1:5,1)表示第一段子渠道的渠道参数(1糙率，2底坡，3底宽，4边坡，5渠段长度）; 
% N1：计算断面数目，N1是个维数为Nout+1的向量   
% Nout:渠段取水口数，或不同断面交界点数，标量
% Qout:取水口流量，是个维数为Nout的向量
% TZ0:渠底高程（sum(N1)维数组）;
% Zouth:取水口水位降落，是个维数为Nout的向量
% Q0:最下游流量;
% Hd:最下游水深
delta=0.7;                             %当量粗糙度，管道有抹灰面层并经抛光
mu=1.31*10^(-6);             %   水温10°时的粘滞系数
lambd=0.02;                           %设定初始值
lambd2=0.03;
N=sum(N1);                                     %原N=max（sum(N1))15-09-18改
Q=zeros(1,N);                                  % 声明变量
Z=zeros(1,N);                                  % 声明变量
J2=N;                                    % J1:赋值起点，J2:赋值终点
if I ==25
    for K=Nout+1:-1:1                     % K:一个渠段内的分段数
        J1=J2-N1(K)+1;
        if K==Nout+1
            Q(J1:J2)=Q0*ones(1,J2-J1+1);           %Q0：该渠段下游流量
            Hout=Hd0;                        % Hd0=Hd1二分法的下边界
            %Hout=Z(1,(N))+Zouth(K)-TZ0(J2);                        %Hout:子渠道最下游断面水深；
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);           %恒定非均匀流下各断面水深
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);           
        elseif K==Nout
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)是K子渠段末端的水头损失
            Hout=Z(J2+1)+Zouth(K)-TZ0(J2);
            Re=4*Q(J2)/2/(3.14259*mu*CANAL(3,K));
            while abs(lambd-lambd2)>0.0001
                lambd2=lambd;
                lambd1=(-2*log10(delta/(3.71*CANAL(3,K))+2.51/(Re*sqrt(lambd))))^(-2);
                lambd=lambd1;
            end
            hf=lambd*CANAL(5,K)*(Q(J2)/2)^2*8/(3.1315*3.1415*9.81*(CANAL(3,K)^5));
            Z(J1:J2)=(Hout+hf) : -hf / ( N1(K)-1 ) : Hout;
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        else
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)是K子渠段末端的水头损失
            Hout=Z(J2+1)+Zouth(K)-TZ0(J2);                        % Hout=下一个子渠道最上游计算断面的水位+水头损失-子渠道最下游计算断面渠底高程     ？  %%% 加上流速水头
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);         %Q(J2)：K子渠段流量
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        end
        J2=J1-1;
    end
else
    for K=Nout+1:-1:1                     % K:一个渠段内的分段数
        J1=J2-N1(K)+1;
        if K==Nout+1
            Q(J1:J2)=Q0*ones(1,J2-J1+1);           %Q0：该渠段下游流量
            Hout=Hd0;                        % Hd0=Hd1二分法的下边界
            %Hout=Z(1,(N))+Zouth(K)-TZ0(J2);                        %Hout:子渠道最下游断面水深；
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);           %恒定非均匀流下各断面水深
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        else
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)是K子渠段末端的水头损失
            Hout=Z(J2+1)+Zouth(K)-TZ0(J2);                        % Hout=下一个子渠道最上游计算断面的水位+水头损失-子渠道最下游计算断面渠底高程     ？  %%% 加上流速水头
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);         %Q(J2)：K子渠段流量
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        end
        J2=J1-1;                                                 %
    end
end
