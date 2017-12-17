function Z=Ini_Z(CANAL,N1,Nout,Qout,TZ0,Zouth,Q0,Hd0, I )
% Calculate the water surface profile of a canal pool with constant flow as the initial condition of the simulation
% Z����������sum(N1)ά���飩
% CANAL�����β������� ����CANAL(1:5,1)��ʾ��һ������������������(1���ʣ�2���£�3�׿�4���£�5���γ��ȣ�; 
% N1�����������Ŀ��N1�Ǹ�ά��ΪNout+1������   
% Nout:����ȡˮ��������ͬ���潻�����������
% Qout:ȡˮ���������Ǹ�ά��ΪNout������
% TZ0:���׸̣߳�sum(N1)ά���飩;
% Zouth:ȡˮ��ˮλ���䣬�Ǹ�ά��ΪNout������
% Q0:����������;
% Hd:������ˮ��
delta=0.7;                             %�����ֲڶȣ��ܵ���Ĩ����㲢���׹�
mu=1.31*10^(-6);             %   ˮ��10��ʱ��ճ��ϵ��
lambd=0.02;                           %�趨��ʼֵ
lambd2=0.03;
N=sum(N1);                                     %ԭN=max��sum(N1))15-09-18��
Q=zeros(1,N);                                  % ��������
Z=zeros(1,N);                                  % ��������
J2=N;                                    % J1:��ֵ��㣬J2:��ֵ�յ�
if I ==25
    for K=Nout+1:-1:1                     % K:һ�������ڵķֶ���
        J1=J2-N1(K)+1;
        if K==Nout+1
            Q(J1:J2)=Q0*ones(1,J2-J1+1);           %Q0����������������
            Hout=Hd0;                        % Hd0=Hd1���ַ����±߽�
            %Hout=Z(1,(N))+Zouth(K)-TZ0(J2);                        %Hout:�����������ζ���ˮ�
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);           %�㶨�Ǿ������¸�����ˮ��
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);           
        elseif K==Nout
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)��K������ĩ�˵�ˮͷ��ʧ
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
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)��K������ĩ�˵�ˮͷ��ʧ
            Hout=Z(J2+1)+Zouth(K)-TZ0(J2);                        % Hout=��һ�������������μ�������ˮλ+ˮͷ��ʧ-�����������μ���������׸߳�     ��  %%% ��������ˮͷ
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);         %Q(J2)��K����������
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        end
        J2=J1-1;
    end
else
    for K=Nout+1:-1:1                     % K:һ�������ڵķֶ���
        J1=J2-N1(K)+1;
        if K==Nout+1
            Q(J1:J2)=Q0*ones(1,J2-J1+1);           %Q0����������������
            Hout=Hd0;                        % Hd0=Hd1���ַ����±߽�
            %Hout=Z(1,(N))+Zouth(K)-TZ0(J2);                        %Hout:�����������ζ���ˮ�
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);           %�㶨�Ǿ������¸�����ˮ��
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        else
            Q(J1:J2)=(Q0+sum(Qout(K:Nout)))*ones(1,J2-J1+1);        %Zouth(K)��K������ĩ�˵�ˮͷ��ʧ
            Hout=Z(J2+1)+Zouth(K)-TZ0(J2);                        % Hout=��һ�������������μ�������ˮλ+ˮͷ��ʧ-�����������μ���������׸߳�     ��  %%% ��������ˮͷ
            Z(J1:J2)=sf_h_d(CANAL(:,K),N1(K),Q(J2),Hout);         %Q(J2)��K����������
            Z(J1:J2)=Z(J1:J2)+TZ0(J1:J2);
        end
        J2=J1-1;                                                 %
    end
end
