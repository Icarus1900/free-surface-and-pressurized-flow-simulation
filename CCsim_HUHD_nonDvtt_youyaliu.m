% Simulation and Control of multi canal systems
% Originally written by Guanghua Guan, Wei Cui and Jie Fan inspired by Changde Wang  @ 2006
% Upgraded By Xiong Yao, Zhiliang Ding, Mengkai Liu, Guoqiang Liu and Yibo Yan
clear;
clc;
tic

% Part I : Data input
% ����****************************************����ϵͳ�趨����****************************************����
NMinToHou = 60 ;                        % ʱ����ƣ�1Сʱ=60���ӣ�1����=60��
NHouToDay = 24 ;                        % 1��=24Сʱ
Limitlength = 200 ;                     % �������ض�������������ʱ��̳���
SPRaise = 0.1 ;                         %���ڿ���Ŀ��ˮλ̧��ֵ
ShowIceTime = 360;                      %���������ʱ����
CalcAcc = 1e-14;                         %���㾫�ȣ�ԭֵ 1e-15
% % ����**************************************����ϵͳ�趨����**************************************����
% disp ( 'ѡ��������ݶ��뷽ʽ 1 = txt; 2 = Excel; ����� = database' )
% Choice = input (' ����������ѡ��:','s');
% switch Choice
%     case '1'
%         [ Ncnl, Ncout, Hu1, Hd1, DT, SV_DT, Tmax, Wt, G, thta, Vmax, DBGate, DBWatLev, T_MDS1, T_MDS2, MDS1, MDS2, HeatCapa, LatentHeat, IceDen, WatDen, DTIce, UpTemp,......
%             IniRough, FinlRough, simtype, SilenceMode, GateAccuracy, GateDischargeAcc, CAAconst, IniIceThick, PID, Kf, tk, PIDS, CanalAccidentUp, CanalAccidentDn,......
%             FailTime, Nout, Nodtype, Kouth, Zoutd, N1, CANAL, z0_u, z0_d, Gb, CDGATE, Gmax, Qmax, Hd01, ZJU, ZJD, Temp, CloudAmount, WindVelo, radi, Tw, Ksn, Csn,......
%             Qoutdsgn, Qoutype, Tstart, Tend,Qstart, Qend,Qout0, Qdown_start, Qdown_end, Tdown_start, Tdown_end ] = dataTXT ( NMinToHou, NHouToDay, SPRaise ) ;
%     case '2'
% [ Ncnl, Ncout, Hu1, Hd1, DT, SV_DT, Tmax, Wt, G, thta, Vmax, DBGate, DBWatLev, T_MDS1, T_MDS2, MDS1, MDS2, HeatCapa, LatentHeat, IceDen, WatDen, DTIce, UpTemp,......
%     IniRough, FinlRough, simtype, SilenceMode, GateAccuracy, GateDischargeAcc, CAAconst, IniIceThick, PID, Kf, tk, PIDS, CanalAccidentUp, CanalAccidentDn,......
%     FailTime, Nout, Nodtype, Kouth, Zoutd, N1, CANAL, z0_u, z0_d, Gb, CDGATE, Gmax, Qmax, Hd01, ZJU, ZJD, Temp, CloudAmount, WindVelo, radi, Tw, Ksn, Csn,......
%     Qoutdsgn, Qoutype, Tstart, Tend,Qstart, Qend,Qout0, Qdown_start, Qdown_end, Tdown_start, Tdown_end ,IfTermate] = dataEXCEL ( NMinToHou, NHouToDay, SPRaise ) ;
%     otherwise
%         [ Ncnl, Ncout, Hu1, Hd1, DT, SV_DT, Tmax, Wt, G, thta, Vmax, DBGate, DBWatLev, T_MDS1, T_MDS2, MDS1, MDS2, HeatCapa, LatentHeat, IceDen, WatDen, DTIce, UpTemp,......
%             IniRough, FinlRough, simtype, SilenceMode, GateAccuracy, GateDischargeAcc, CAAconst, IniIceThick, PID, Kf, tk, PIDS, CanalAccidentUp, CanalAccidentDn,......
%             FailTime, Nout, Nodtype, Kouth, Zouth, Zoutd, N1, CANAL, z0_u, z0_d, Gb, CDGATE, Gmax, Qmax, Hd01, ZJU, ZJD, Temp, CloudAmount, WindVelo, radi, Tw, Ksn,......
%             Csn, Qoutdsgn, Qoutype, Tstart, Tend, Qstart, Qend, Qout0, Qdown_start, Qdown_end, Tdown_start, Tdown_end ] = dataDB ( NMinToHou, NHouToDay, SPRaise ) ;
% end
try load InputSIM
    disp ( '���������Ƿ��иĶ�?' )
    cprintf('red','���иĶ��밴��y�� ���س������������ݣ�');
    cprintf('Comments','ֱ�Ӱ��س�����ʼ���棻');
    Choice = input ('��������ѡ��:','s');
    switch Choice
        case 'y'
            disp ('�������ļ��������ݣ�');
%             delete InputSIM.mat
            dataEXCEL ( NMinToHou, NHouToDay, SPRaise ) ;
            load InputSIM
        otherwise
            disp ( '���õ�ǰĿ¼�е��������ݿ�ʼ���棡' ) ;
    end
catch
    disp ( '��ǰĿ¼��û���������ݣ����ļ����룡' ) ;
    dataEXCEL ( NMinToHou, NHouToDay, SPRaise ) ;
    load InputSIM
end
% ����***************************************�˻�������Ϣ*****************************************����

if IfTermate==1;
    return
end
% ����***************************************�˻�������Ϣ*****************************************����

% Part II : defination of variable and matrix
Zouthini = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;        %������ֲ�ˮͷ��ʧ����������
Zouthqtt = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;        %������ֲ�ˮͷ��ʧ��tʱ�̵�ˮͷ��ʧ��
Zouthunst = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;       %������ֲ�ˮͷ��ʧ����ʼʱ��ˮͷ��ʧ��
Len_Canal = zeros ( Ncnl , 1 ) ;                            %���س���
V_wave = zeros ( Ncnl , 1 ) ;                               %ˮ���ٶ�
T_wave = zeros ( Ncnl , 1 ) ;                               %ˮ������ʱ��
TZ0 = zeros ( Ncnl , max ( sum (  N1 ( : , 1 : Ncnl ) ) ) ) ;   %����������������׸߳�
Z0 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;     %�������������ˮλ
Qt_d = zeros ( Ncnl , 1 ) ;                                 % �����μƻ������������� (Q0)
Qsum_down = 0 ;                                             % ��¼������ˮ����
QP1 = zeros ( Ncnl , 2 ) ;                                  % ���������������߽�����,����������բ�Ź�բ����
QPn = zeros ( Ncnl , 1 ) ;                                  % ���������������߽�����������������բ�Ź�բ����
GA = zeros ( Ncnl , 2 ) ;                                   % բ�ſ���
N = zeros ( Ncnl , 1 ) ;                                    %һ�������ڵļ����������
Q = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;  %��������������
Qtt = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;%��������������
Z = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;  %����������ˮλ
Z1 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;
Z2 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;
Vini = zeros ( Ncnl , max ( Nout( : ) ) + 1 );              %�����������ʼƽ������%%%%%%%%%%%%%%%%%%%%%%�Լ��ӵ�
Vqtt = zeros ( Ncnl , max ( Nout( : ) ) + 1 ) ;             %�������������γ�ˮλ����ʱƽ������%%%%%%%%%%%%%%%%%%%%%%%�Լ��ӵ�
Vunsteady = zeros ( Ncnl , max ( Nout( : ) ) + 1 ) ;        %�Ǻ㶨�����������ٵļ���
Hd0 = zeros ( Ncnl , 1 ) ;                                  %���γ�ˮλ�������ζ˳�ʼ��ˮ�
V = zeros ( Ncnl , 1 ) ;                                    %������ˮ�����
Delay = zeros ( Ncnl , 1 ) ;                                %�����ͺ�ʱ��
Velo = zeros ( (max ( Nout)+ 1 ) , Ncnl ) ;                 %ÿ�������ڸ�����ƽ������
Velomax = zeros ( Ncnl , 1 ) ;                              %�������������
VeloENG = zeros (  max ( Nout ( : ) )+ 1 , Ncnl ) ;         %��������
% ����------------------------------------��¼������ˮ�����---------------------------------------����
Tmaxm = Tmax * NMinToHou ;                                  % ������ʱ��(����)
Volu = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                    %��ͬʱ��������ˮ�����
Delayl = zeros ( Ncnl , Tmaxm / DT + 1 ) ;
Qtal = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                    %��¼ǰ������
% ����------------------------------------��¼������ˮ�����---------------------------------------����
% ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ڳ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
for I = Ncnl : -1 : 1
    N ( I ) = sum ( N1 ( : , I ) ) ;                        %һ�������ڵļ��������
    % �����йر���������,���ڲ������迼�Ǳ��������Ƿ���Ҫ�������������ϴ�
    hi ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % ���Ǻ��
    %Tw ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % ����ƽ��ˮ��
    Tws ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;   	% �����ϱ����¶�
    Ca ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % �������ܶ�
    Ci ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % ����Ũ��
    TTT ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % �ⶳʱ�̣���Ҫ����Ĵ������
    Ciaver ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;    % �ռ䲽���ϵ�ƽ������Ũ��
    CAA ( I ) = zeros ( 1 ) ;                               % �������ܶȷⶳ��׼
    d0 ( I ) = zeros ( 1 ) ;                                % �趨����������
    Vi ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % ����
    u1max ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;     % u����ƽ������
    % �����ʼֵ
    % !!�����ʼ���������⣬����¶ȵ������ʱ���Ǻ�ȡ������ܶȵȲ�����һ��ȫ���㡣
    %Tw ( I , 1 : N ( I ) ) = 1 ;                           % ˮ�³�ʼ����
    hi ( I , 1 : N ( I ) ) = 0 ;                            % ���Ǻ��
    Ci ( I , 1 : N ( I ) ) = 0 ;                            % ����Ũ��
    Ca ( I , 1 : N ( I ) ) = 0 ;                            % �������ܶ�
    Vii1 ( I , 1 : N ( I ) ) = 0 ;                          % ���ǿ�϶���
end
% ����------------------------------------������Ƿⶳ�ʱ�����ʼ��---------------------------------------����
PercLength = zeros ( Ncnl , Tmaxm / DT + 1 );
LengthIce = zeros ( Ncnl , Tmaxm / DT + 1);
% ����------------------------------------������Ƿⶳ�ʱ�����ʼ��---------------------------------------����
% ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ڳ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
Qt_u = zeros ( Ncnl , 1 ) ;                                 % �����μƻ������������� (Q1)
YT = zeros ( Ncnl , 1 ) ;                                   % ���Ƶ�㶨��Ŀ��ˮλ
Qsum_up = zeros ( Ncnl + 1 , 1 ) ;                          % ���������������ۼ�ˮ��(Qsum)
Huall = zeros ( Tmaxm / DT , Ncnl ) ;                       % ��ͬʱ����ϵ�բǰˮ��
Hdall = zeros ( Tmaxm / DT , Ncnl ) ;                       % ��ͬʱ����ϵ�բ��ˮ��
Hdiff = zeros ( Tmaxm / DT , Ncnl ) ;                       % ��ͬʱ����ϵ�բǰբ��ˮ��֮��
Q0_d_0= zeros ( Ncnl , 1 ) ;
V0_u = zeros ( Ncnl , 1 ) ;
Tx = zeros ( Tmaxm / DT + 1 , 1 ) ;                           % ʱ�����ݺ�����
Q_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % �����������������
Qoutsum = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                   % ����������ȡˮ���������
Z_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % ����������ˮλ���
Z_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % ����������ˮλ���
Z_M = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % �������е����ˮλ���
GA1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % ����������բ�ſ������
V1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                        % ������ˮ��������
YT1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % ���Ƶ�㶨��Ŀ��ˮλ���
YF1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % ���Ƶ�ʵ��ˮλ���
YTorig_U = zeros ( Ncnl , 1 ) ;                               % �㶨������Ŀ��ˮλ
YTorig_D = zeros ( Ncnl , 1 ) ;                               % �㶨������Ŀ��ˮλ
YTorig_UP = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                 % �㶨������Ŀ��ˮλ���
YTorig_DOWN = zeros ( Ncnl , Tmaxm / DT + 1 ) ;               % �㶨������Ŀ��ˮλ���
Velocity = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                  % ������������������ı仯����
 Velocity1N = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                  % ���������׹�բ�������ٴ�С�仯���
Vttl = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                      % ��������ǰ��Ŀ����������ֵ���
E = zeros ( Ncnl , 2 ) ;                            % ���Ƶ�ˮλ���
Ef = zeros ( Ncnl , 2 ) ;                           % �˲���Ŀ��Ƶ�ˮλ���
EC = zeros ( Ncnl , 2 ) ;                           % ���Ƶ�ˮλ���仯��
ECf = zeros ( Ncnl , 2 ) ;                          % �˲���Ŀ��Ƶ�ˮλ���仯��
Jug1 = zeros ( Ncnl , 1 ) ;                         % ������ˮλ���ָ�꣨ISE ��������ָ�꣩�м����
IAE1 = zeros ( Ncnl , 1 ) ;                         % ���������ˮλָ���м����
IAQ1 = zeros ( Ncnl , 1 ) ;                         % �������������ָ���м����
IAW1 = zeros ( Ncnl , 1 ) ;                         % ���������բ��ָ���м����
IAQ = zeros ( Ncnl , 1 ) ;                          % �������������ָ��
IAW = zeros ( Ncnl , 1 ) ;                          % ���������բ��ָ��
DQGA_sum = zeros ( Ncnl , 1 ) ;                     % ���������е����ι�բ�����仯�ۻ��������ж�ϵͳ�Ƿ�ﵽ�ȶ���
QdownT = zeros ( Tmaxm / DT + 1 , 1 ) ;             % Tʱ�̵�ϵͳ�����ζ�ȡˮ����
YF = zeros ( Ncnl , 1 ) ;                             % ���Ƶ�Ǻ㶨��Ŀ��ˮλ
Vtt = zeros ( Ncnl , 1 ) ;                            % ��������ǰ��Ŀ����������
DQGA = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % ������ΪDT������SV_DTΪ������������������������޸�
Ddz_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ����������Уʵ��
Ddz_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ����������Уʵ��
Ddz_M = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ����������Уʵ��
EH_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % ����������ˮλ���
EH_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % ����������ˮλ���
EH = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                % ����������ˮλ���
DV_U1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������1h����ˮλ��
DV_D1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������1h����ˮλ��
DV_M1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������1h����ˮλ��
DV_U2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������24h����ˮλ��
DV_D2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������24h����ˮλ��
DV_M2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % ������24h����ˮλ��
DQsum_up = zeros ( Ncnl + 1 , 1 ) ;                   % �����������ι�ˮ��֮��
LengthPool = zeros ( Ncnl , 1 ) ;
len_1d = zeros ( Ncnl , 1 ) ;
len_3d = zeros ( Ncnl , 1 ) ;
len_5d = zeros ( Ncnl , 1 ) ;
DQff= zeros ( Ncnl , Tmaxm / DT + 1 ) ;
DQfb= zeros ( Ncnl , Tmaxm / DT + 1 ) ;
Dvttl= zeros(Ncnl,Tmaxm / DT);
% Part III : Initialize of variables
SV_N = DT/SV_DT ;           %ÿ�ο��ƶ���ʱ���ڵķǺ㶨������ʱ�䲽��
% ����***************************************������������*****************************************����
% ����--------------------------�������������ˮ�������δ�������������ʱ��--------------------------����
for I = 1 : Ncnl
    if I==51
        I=I;
    end
    Len_Canal ( I ) =  sum( CANAL ( 5 , 1 : Nout ( I ) + 1 , I ) ) ;          % ���س���
    V_wave ( I ) = ( G * Hd01 ( I ) ) ^ 0.5 * NMinToHou ;                        % ��λ��m/min  %"Hd01"���γ�ˮλ�������ζ˳�ʼ��ˮ�
    T_wave ( I ) = Len_Canal ( I ) / V_wave ( I ) ;   %������ʱ��
    T_wave ( I ) = round ( T_wave ( I ) / DT ) * DT ;                    % ȡ��
    % ����--------------------------�������������ˮ�������δ�������������ʱ��--------------------------����
    % ����--------------------------------------�������׸߳�TZ0---------------------------------------����
    zt1 = z0_u ( I ) ;                                      % %�������������߳�
    Jst = 1 ;                                                % J1:��ֵ��㣬J2:��ֵ�յ�
    for K = 1 : Nout ( I ) + 1                                    % K:һ�������ڵķֶ��� Nout:�����νڵ����
        Jend = Jst + N1 ( K , I ) - 1 ;                                    % N1�� ���ڵ�֮��������ĸ���
        zt2 = zt1 - CANAL ( 2 , K , I ) * CANAL ( 5 , K , I ) ;
        TZ0 ( I , Jst : Jend ) = zt1 : ( zt2 - zt1 ) / ( Jend - Jst ) : zt2  ; %TZ0:����������������׸̣߳�������������ȷֲ�
        Z0 ( I , Jst : Jend ) = TZ0 ( I , Jst : Jend ) + Hd01 ( I ) ;
        if  K~= Nout ( I ) + 1
            zt1 = zt2 - Zoutd ( K , I ) ; %zt1:�������׸߳�ʱ���������߳�;zt2:�������׸߳�ʱ�������յ�߳�  % Zoutd����ȡˮ��(������)��ʼ���׽���
            Jst = Jend + 1 ;
        end
    end
    % ����--------------------------------------�������׸߳�TZ0---------------------------------------����
end
% ����***************************************������������*****************************************����
% ����--------------------------ȡˮ�ڱ仯ʼ��ʱ�估������ʼ��--------------------------����
if Qdown_start == Qdown_end                     %���ϵͳĩ�������ޱ仯�������仯ʱ�������ȡˮ�ڱ仯 ����ʼʱ�̷�ˮ����=ĩ��ʱ�̷�ˮ����
    TS = min ( Tstart );
    TE = max ( Tend );
else                                            %���ϵͳĩ�������б仯�������仯ʱ��Ӧ����ȡˮ�ں�ĩ��ʱ��ĵ���
    TS = min ( min ( Tstart ) , Tdown_start );  %�����仯��ʼʱ��
    TE = max ( max ( Tend ) , Tdown_end ) ;     %�����仯����ʱ��
end
for I = Ncnl : -1 : 1                           %�������γ�ʼ����
    if I == Ncnl
        Qt_d ( I ) = Qdown_start ;                %Qdown_start�����γ�ʼ����     Qt_d ( I )����I�����ε����γ�ʼ����
    else
        Qt_d ( I ) = Qt_d ( I + 1 ) + sum ( Qout0 ( : , I + 1 ) );         %Qout0�� ��ȡˮ�ڳ�ʼ����
    end
end
Qout1 = Qout0 ;                         % ��¼��һ�ε�Qout0
Q01_d = Qt_d ( Ncnl ) ;                 % ��¼��һʱ��������ˮ�����仯
% ����--------------------------ȡˮ�ڱ仯ʼ��ʱ�估������ʼ��--------------------------����
% ����------------------------------���빤���ļ��������ʼբ������----------------------------------����
NodtypeGate ( 1 : Ncnl ) = 5 ;    % Nodtype2��ʲô��˼��
Nodtype2 = vertcat ( Nodtype , NodtypeGate );
% for I = 1 : Ncnl
%     Nodtype2 ( Nout ( I ) + 1 , I ) = 5 ;         % բ�Žڵ㣬��Nout(I)+1�����������ټ������Ҫ
% end

% Part IV : Initial state calculation of canal system
% ����---------------------------------------������ˮλ��ʼ��������--------------------------------����
for I = 1 : Ncnl
    % ����---------------------------------------������ʼ����-------------------------------------����
    Q ( I , 1 : N ( I ) ) =feval(@ Ini_Q , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , Qt_d ( I ) ) ; %�����������������=����������+����������ȡˮ������
    % ����---------------------------------------������ʼ����-------------------------------------����
    % ����---------------------------------------ˮλ��ʼ����-------------------------------------����
    % ����-------------------------------------�ֲ�ˮͷ��ʧ����-----------------------------------����
    Vini ( I , 1 : Nout ( I ) + 1 ) = feval(@VSegmentavg , CANAL , Q , Z0 , TZ0 , I , N1 , Nout( I ) ) ;                           %Z0��������ˮ����
    Zouthini (  1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vini , Nout ( I ) , I );
    % ����-------------------------------------�ֲ�ˮͷ��ʧ����-----------------------------------����
    % ����--------------------------------���γ�ˮλ��������ˮ�����-------------------------------����
    if Wt == 0
        Hd0 ( I ) = Hd01 ( I ) ;
    end
    % ����--------------------------------���γ�ˮλ��������ˮ�����-------------------------------����
    % ����---------------��������У�ͨ�����ˮ����������������������γ�ʼ����ˮ��----------------����
    if Wt == 0.5
%         for J = 1 : N ( I )
%             Z ( I , J ) = Hd01 ( I ) + TZ0 ( I , J ) ;
%         end
        Z ( I , : ) = Hd01 ( I ) *ones (  1 , N ( I ) ) + TZ0 ( I , : );
        [ V( I ) , Delay( I )] = feval(@Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I )+1, I ) , N1 ( 1 : Nout ( I )+1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) ,Q(I,1:N(I))) ;      % ���ڽ�����ȵ�Ӱ�죬��ʵ�ʵ���ʵ��������
        Hd0 ( I ) =feval(@ V2Hd , V ( I ) , CANAL( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthini ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , [ 2 8 ] ) ;   % ��������и��������ζ˳�ʼ��ˮ�
    end
    % ����---------------��������У�ͨ�����ˮ����������������������γ�ʼ����ˮ��----------------����
    % ����---------------�ݴ������ʼʱ��ˮ����---------------------------------------------------����
    Z1 ( I , 1 : N ( I ) )=feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthini ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I ) , I ) ;
    Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z1 , TZ0 , I , N1 , Nout ( I ) ) ;               % ���������������������
    Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    Z2 ( I , 1 : N ( I ) ) =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I )  , I ) ;     %15-9-26��
    Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z2 , TZ0 , I , N1 , Nout ( I ) ) ;               % ���������������������
    Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    if max ( Z1 - Z2 ) > CalcAcc;
        Z1 = Z2 ;
        Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z1 , TZ0 , I , N1 , Nout ( I ) ) ;               % ���������������������
        Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
        Z2 ( I , 1 : N ( I ) ) =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I ) , I  ) ;     %15-9-26��
        Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z2 , TZ0 , I , N1 , Nout ( I ) ) ;               % ���������������������
        Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    else
        Z = Z2;
    end
    Z11 = Z ( 1 , 1 : N ( 1 ) );
    % ����---------------------------------------ˮλ��ʼ����-------------------------------------����
    % ����---------------------------------------�������ټ���-------------------------------------����
    for J = 1 : Nout ( I ) + 1
        Velo ( J , I ) = feval ( @Vcmax , CANAL , Q , Z , TZ0 , I , J , N1 ) ;               % ������������������
    end
    Velomax ( I ) = feval ( @F_velomax , VeloENG (1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Limitlength , Nodtype2 ( 1 : Nout ( I ) + 1 , I ) , Velo ( 1 : Nout ( I ) + 1 , I ) , CANAL ( 5 ,1 : Nout ( I ) + 1 , I ) );
    % ����---------------------------------------�������ټ���-------------------------------------����
end
% ����---------------------------------------������ˮλ��ʼ��������--------------------------------����

% ����--------------------------------------------��ʼբ�ſ���------------------------------------����
for I = 1 : Ncnl
    if I==30
        I=I;
    end
    if I == 1
        Q0_d_0 ( I , 1 ) = Qt_d ( I ) + sum ( Qout0 ( : , I ) );
        V0_u ( I , 1 ) =Q0_d_0 ( I , 1 )/ ( CANAL ( 3 ,1 ,I ) + CANAL ( 4 ,1 ,I ) * ( Z ( I ,1 ) - TZ0 ( I , 1 ) ) )/( Z ( I , 1 ) - TZ0 ( I , 1 ) );
%         Hu = Hu1 - Kouth ( 1 , 1 ) * V0_u( I , 1 )^ 2 /2/G ;%��һ���Σ����أ�����ˮ��            %Kouth�� ���ڵ�ˮͷ��ʧϵ��
        Hu = Hu1;
        Hd = Z ( I , 1 ) - z0_u ( I ) ;%+ Kouth ( Nout ( I -1) + 1 , I ) *  V0_u ( I ,1 )^ 2/2/G ;%����ˮ��   Hu1��������ˮ��ˮ��
    else
        Q0_d_0 ( I , 1 ) = Qt_d ( I - 1 ) ;
        V0_u ( I , 1 ) =feval(@ VPoolavg , CANAL , Q0_d_0 , Z , TZ0 , I-1 , N1 , Nout ( I-1 ) ) ;               %15-9-22��
        Hu = Z ( I - 1 , N ( I - 1 ) ) - z0_d ( I - 1 ) ;
        Hd = Z ( I , 1 ) - z0_d ( I - 1 ) + Kouth ( Nout ( I -1) + 1 , I-1 ) *  V0_u ( I ,1 )^ 2/2/G ;          %z0_d ( I - 1 ) ��һ�������յ�����׸߳�  Z ( I , 1 )����һ����������ˮ��߳�
    end
    Huall ( 1 , I ) = Hu ;
    Hdall ( 1 , I ) = Hd ;
    Hdiff ( 1 , I ) = Hu - Hd ;
    
    
    if I==1
        Nodtypelogical=( Nodtype( :, I )==1);
        Qt_u ( I ) = Qt_d ( I ) +sum( Qout0 ( : , I ).*  Nodtypelogical);        %�����ط�ˮ���żƻ���ˮ��
    else
        
        Qt_u ( I ) = Qt_d ( I-1 );
    end
    if I == 1                                            % GateDischarge����բ�ſ���
        Au = ( CANAL ( 3 , 1 , I ) + CANAL ( 4 , 1 , I ) * Hu ) * Hu ; %Au����������բǰ��ˮ���
    else
        Au = ( CANAL ( 3 , Nout ( I - 1 ) + 1 , I - 1 ) + CANAL ( 4 , Nout ( I - 1 ) + 1 , I - 1 ) * Hu ) * Hu ;
    end
    GA ( I , : ) =  feval(@Q2Greal , I  , 0 , Hu , Hd , Qt_u ( I ) , 0 , Gb ( I ) , Au , Qt_u ( I ) , Gmax ( I ) , G , CDGATE ( I ) , DBGate , 1000 , repmat([1000 1000],[Ncnl 1]) , DBWatLev , SilenceMode ,  CalcAcc ) ;
    GA ( I , : ) =  round ( GA ( I , : ) * GateAccuracy  ) / GateAccuracy  ; %GA����������բ�ſ���
    QP1 ( I , : ) = Qt_u ( I ) ;                               % �����γ�ʼ�ϱ߽�
    QPn ( I , : ) = Qt_d( I ) ;                               % �����γ�ʼ�±߽�
    YT ( I ) = Hd01 ( I ) + z0_d ( I ) ; %YT���Ƶ�㶨��Ŀ��ˮλ
    [ V( I ) , Delay( I )] = feval(@Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 : N ( I ) ) ) ;% ������ˮ�����
    Volu ( I , 1 ) = V ( I ) ;
    Delayl ( I , 1 ) = Delay ( I ) ;
    Qtal ( I , 1 ) = Qt_u ( I ) ;
    Qsum_up ( I ) = 0 ;                                   % ���������������ۼ�ˮ��
end
% ����--------------------------------------------��ʼբ�ſ���------------------------------------����

% ����-----------------------------------------��¼���ݺͲ�������---------------------------------����

for I = 1 : Ncnl
    Tx ( I , 1 ) = 0 ;                                        % ʱ�����ݺ�����
    Q_U ( I , 1 ) = QP1 ( I , 1 ) ;                           % �����������������
    Qoutsum ( I , 1 ) = sum ( Qout0 ( : , I ) ) ;             % ����������ȡˮ���������
    Z_U ( I , 1 ) = Z ( I , 1 ) ;                             % ����������ˮλ���
    Z_D ( I , 1 ) = Z ( I , N ( I ) ) ;                       % ����������ˮλ���
    Z_M ( I , 1 ) =0.5 * ( Z_U ( I , 1 )+ Z_D ( I , 1 ));
    GA1 ( I , 1 ) = GA( I , 1 ) ;                             % ����������բ�ſ������
    V1 ( I , 1 ) = V ( I ) ;                                  % ������ˮ��������
    Vttl ( I , 1 ) = V ( I ) ;                                % ��������ǰ��Ŀ����������ֵ���
    YT1 ( I , 1 ) = YT ( I ) ;                                % ���Ƶ�㶨��Ŀ��ˮλ���
    YTorig_U ( I ) = Z ( I , 1 ) ;                            % �㶨������Ŀ��ˮλ���
    YTorig_D ( I ) = Z ( I , N ( I ) ) ;                      % �㶨������Ŀ��ˮλ���
    YTorig_UP ( I , 1 ) = YTorig_U ( I ) ;                    % �㶨������Ŀ��ˮλ���
    YTorig_DOWN ( I , 1 ) = YTorig_D ( I ) ;                  % �㶨������Ŀ��ˮλ���
    Velocity ( I , 1 ) = Velomax ( I ) ;
end
Zall ( 1 , : , : ) = Z ;                                      %��ͬʱ����ϸ����ε�ˮλ��Zall(5,3,11)��ʾ��5�����ڵ��������ε�11���ڵ��ϵ�ˮλ
Qall ( 1 , : , : ) = Q ;                                      %��ͬʱ����ϸ����ε�ˮλ
% ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ڳ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
Twall ( 1 , : , : ) = Tw ;                                    %��ͬʱ����϶���ƽ��ˮ��
hiall ( 1 , : , : ) = hi ;                                    %�洢��ϵ��������Ǻ�ȵı仯����
Ciall ( 1 , : , : ) = Ci ;                                    %��ͬʱ����ϵ�����Ũ��
% ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ڳ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����

Qta = Qt_u ;                                        % Ŀ������������ǰ���ã�
% ����-----------------------------------------��¼���ݺͲ�������---------------------------------����
% ����*********************************************��ʼ����***************************************����

% Part V : Unsteady flow with canal control simulation
% ����*********************************************���㲿��***************************************����
DB = Vmax * DT ;                                      %һ������������բ�����ֵ����
T = 0 ;                                               %ʱ�䣨���ӣ�
T1 = 0 ;                                              % ������ʼ��
if Tmax <= 72
    Tshow = 60 ;                             % pops time every 60 minutes
else if Tmax <= 720
        Tshow = 360 ;                        % pops time every 6 hours
    else
        Tshow = 1440 ;                       % Pops time every day
    end
end
while T < Tmaxm     %Tmaxm������ʱ�䣨���ӣ�
    T1 = T1 + 1 ;
    T = T + DT ;    %DT����ʱ�䲽����min
  if T==1440
      T=T;
  end
    for I = Ncnl : -1 : 1           %������ϵͳ���������仯ʱ�ڹر�ǰ������
        %         if T < TE * NMinToHou
        %             PID ( I , 1 : 3 ) = 0 ;
        %         else
        PID ( I , 1 : 3 ) = PIDS ( I , 1 : 3 ) ;
        %         end
    end
    if Tshow == 60
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['Ŀǰ����ʱ����е���' , num2str( T ) , '���ӣ�������ʱ��Ϊ' , num2str( Tmax ) , 'Сʱ�� =' , num2str( Tmaxm ) , '���ӣ� ��������С�' ] ) ;
        end
    elseif Tshow == 360
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['Ŀǰ����ʱ����е���' , num2str( T / NMinToHou ) , 'Сʱ��������ʱ��Ϊ' , num2str( Tmax ) , 'Сʱ����������С�' ] ) ;
        end
    else
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['Ŀǰ����ʱ����е���' , num2str( T / ( NMinToHou * NHouToDay ) ) , '�죬������ʱ��Ϊ' , num2str( Tmax / NHouToDay ) ,'�죬���������===================================' ] ) ;
        elseif floor ( T / 180 ) == T / 180 ;
            disp ( ['Ŀǰ����ʱ����е���' , num2str( floor ( T / ( NMinToHou * NHouToDay ) ) ) ,'��' , num2str( T / NMinToHou - floor ( T / ( NMinToHou * NHouToDay ) ) * NHouToDay ) , 'Сʱ��������ʱ��Ϊ' , num2str( Tmax / NHouToDay ) , '�죬��������С�' ] ) ;
        end
    end
    if T < FailTime
        Qout0 = feval(@Qout , T, Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0 ) ;                           % ����Qout
    else
        Qout0 = AccidentResponse ( Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, T , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout );
    end                          % ����Qout
    Qsum_down = Qsum_down +feval(@ Qdown , T , Tdown_start , Qdown_start , Tdown_end , Qdown_end ) * DT * NMinToHou ;           % ��������ˮ��   ��������������⣬Ҫ��ˮ�������ۻ�ֵ��ʲô��
    % ����---------------------------------�����̬�����仯�������¼���Ŀ������---------------------------����
    for I = Ncnl : -1 : 1
        if I == Ncnl
            Qt_d ( I ) = feval(@Qdown , T , Tdown_start , Qdown_start , Tdown_end , Qdown_end ) ;             % Qt_d�������������� (Qt: Ŀ������)
        else
            Qt_d ( I ) = Qt_d ( I + 1 ) + sum ( Qout0 ( : , I + 1 ) );
        end
    end
    for I = Ncnl : -1 : 1
        if I == 1
            Qt_u ( I ) = Qt_d ( I ) + sum ( Qout0 ( : , I ) );
        else
            Qt_u ( I ) = Qt_d ( I - 1 ); % Qt_u �����μƻ������������� (Q1)
        end     %Qttѭ����ʼ��ĺ㶨������
    end
    for I = Ncnl : -1 : 1
        Qtt ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % �����������������鶨��
        Qtt ( I , 1 : N ( I ) ) =feval(@ Ini_Q , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , Qt_d ( I ) ) ;
        Vqtt ( I , 1 : ( Nout ( I ) + 1 ) ) = zeros ( 1 , Nout ( I ) + 1 ) ;       % ����������ƽ���������鶨��
        Vqtt ( I , 1 : ( Nout ( I ) + 1 ) )  = feval( @VSegmentavg , CANAL , Qtt , Z , TZ0 , I , N1 , Nout ( I ) ) ;               % �����������ƽ������
        Zouthqtt ( 1 : ( Nout ( I ) + 1 ) , I ) = zeros ( Nout ( I ) + 1 , 1 ) ;       % ����������ƽ���������鶨��
        Zouthqtt (  1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vqtt , Nout ( I ) , I );
        %for K = 1 : Nout ( I )
        %    HeadLoss ( K , I ) = ZOUTH ( Nodtype ( K , I ) , Kouth ( K , I ) , Qtt ( sum ( N1 ( 1 : K , I ) ) ) ) ;    %Qttѭ����ʼ��ĺ㶨������
        %end
        if Wt == 0                                    % ���γ�ˮλ����
            Hdtt ( I ) = Hd01 ( I ) ; %Hdttѭ����ʼ��㶨��ˮ���߼��������ˮ��
        else
            Hdtt = feval(@V2Hd , V1 ( I , 1 ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthqtt ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , [ 2 8 ] ) ;  % ���������
        end
        Ztt =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I  ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthqtt ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hdtt ( I ) , I  ) ;
        Zttl (T1+1, I , 1 : N ( I ) ) = Ztt ; %Zttѭ����ʼ��㶨��ˮ���߼���
        YTorig_U ( I ) = Ztt ( 1 ) ;
        YTorig_D ( I ) = Ztt ( N ( I ) ) ;
        YT ( I ) = Hdtt ( I ) + TZ0 ( I , N ( I ) ) ;
        Vtt ( I ) =feval(@ Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Ztt ( 1 : N ( I ) ) , Q ( I , 1 :N ( I ) ) ) ;
        Vttl ( I , T1 + 1 ) = Vtt ( I ) ;
        if Delay ( I ) > DT
            Delay ( I ) = DT ;
        end
        if T < FailTime
%             sumQout = sum ( Qout ( ( T + Delay ( I ) ), Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0)  ) ;                           % sumQout ���ؼ�����ȡˮ������֮��
              sumQout ( 1 : Ncnl , T1 ) = sum ( Qout ( ( T ), Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0)  ) ;
        else
            sumQout ( 1 : Ncnl , T1 ) = sum ( AccidentResponse ( Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, T + Delay ( I ) , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout ) );
        end
        Dvttl ( I , T1 ) = Vttl ( I , T1 + 1 ) - Vttl ( I , T1 );
        if I == Ncnl
            Qta ( I ) = Qt_d ( I ) + sumQout ( I , T1 ); %+  Dvttl ( I , T1 ) / ( DT * NMinToHou * tk ( I ) ) ; %�����׶���������β������+����ȡˮ+�������仯�����ˮ��/ʱ�䲽��
        else
            Qta ( I ) = Qta ( I + 1 ) + sumQout ( I , T1 ) ;%+  Dvttl ( I , T1 ) / ( DT * NMinToHou * tk ( I ) ) ;
        end
        Qtal ( I , T1 + 1 ) = Qta ( I ) ;
        if   T > NMinToHou * TS && T < NMinToHou * ( TE + 3 )   %TS�����仯��ʼʱ��
            Qta ( I ) = Qta ( I ) ;
        else
            Qta ( I ) = Qt_u ( I ) ;
        end
    end
    % ����---------------------------------�����̬�����仯�������¼���Ŀ������---------------------------����
    % ����------------------------------------------�����������������-----------------------------------����
    for I = Ncnl : -1 : 1
        YF ( I ) = Wt * Z(  I, 1 ) + ( 1 - Wt ) * Z(  I, N ( I ) ) ;
        E ( I , 2 ) = YT ( I ) - YF ( I ) ;                        % E(I,2):I���δ�ʱ�̵���E(I,1):I������һʱ�̵���
        % ����---------------------------------------------�˲�-----------------------------------------����
        Ef ( I , 2 ) = ( 1 - Kf ( I ) ) * E( I , 2 ) + Kf ( I ) * E ( I , 1 );    % �˲���ˮλƫ��
        ECf ( I , 2 )=( Ef ( I , 2 ) - Ef ( I , 1 ) ) / DT;            % �˲�|��|ˮλƫ����
        % ����---------------------------------------------�˲�-----------------------------------------����
        EC ( I , 2 ) = ( ( E ( I , 2 ) - E ( I , 1 ) ) ) / DT ;              % �˲�|ǰ|ˮλƫ����
        Jug1 ( I ) = Jug1 ( I ) + ( E ( I , 2 ) ^ 2 ) / ( YT ( I ) - 0.5 * ( z0_u ( I ) + z0_d ( I ) ) ) ;
        IAE1 ( I ) = IAE1 ( I ) + abs ( E ( I , 2 ) ) / ( YT ( I ) - 0.5 * ( z0_u ( I ) + z0_d ( I ) ) ) ;
        [ QP1( I , 2 ) , DQff(I,T/DT) , DQfb(I,T/DT) ]= feval(@FFwithFB , Ef ( I , 2 ) , ECf ( I , 2 ) , ECf ( I , 1 ) , PID ( I , 1 : 3 ) , Qta ( I ) , QP1 ( I , 1 ) , Qmax ( I ) ) ; % PID������������˲�
        
   
        % ����---------------�������仯�ڼ����������ʱ�䲽����ͨ��������ʱ�������仯�������--------------����
        if   T > NMinToHou * TS && T < NMinToHou * ( TE + 3 )
            if fix ( T_wave ( I ) / DT ) == T_wave ( I ) / DT
                QP1 ( I , 2 ) = QP1 ( I , 2 ) ;
            else
                QP1 ( I , 2 ) = QP1 ( I , 1 ) ;
            end
        end
        
       
        
        % ����---------------�������仯�ڼ����������ʱ�䲽����ͨ��������ʱ�������仯�������--------------����
    end
    % ����------------------------------------------�����������������------------------------------------����
    
    % ����---------------------------------------------բ�ſ��ȼ���---------------------------------------����
    Nodtype2 = Nodtype ;
    for I = 1 : Ncnl
        Nodtype2 ( Nout ( I ) + 1 , I ) = 5 ;                % բ�Žڵ㣬��Nout(I)+1�����������ټ������Ҫ
    end
    for I = Ncnl : -1 : 1
        if I==55
            I=I;
        end
        if  I == 1
            %Hu = Hu1 - Kouth ( Nout ( 1 ) + 1 , 1 ) * V0_d ( I , 1 ) ^ 2 /2/G ;
            Hu = Hu1;                     %15-9-22��
            Hd = Z ( I , 1 ) - z0_u ( I ) ;
        else
            Hu = Z ( I - 1 , N ( I - 1 ) ) - z0_d ( I - 1 ) ;
            Hd = Z ( I , 1 ) - z0_d ( I - 1 ) + Kouth ( Nout ( I-1 ) + 1 , I-1 ) * V0_u ( I , 1 ) ^ 2 /2/G  ;            
%             Hu = Z ( I - 1 , N ( I - 1 ) ) - z0_d ( I - 1 )  ;
%             Hd = Z ( I , 1 ) - z0_d ( I - 1 ) + Kouth ( Nout ( I ) + 1 , I ) * V0_u ( I , 1 ) ^ 2 /2/G ;
        end
        Huall ( T1 , I ) = Hu ;
        Hdall ( T1 , I ) = Hd ;
        Hdiff ( T1 , I ) = Hu - Hd ;
        if I == 1
            Au = ( CANAL ( 3 , 1 , I ) + CANAL ( 4 , 1 , I ) * Hu ) * Hu ;
        else
            Au = ( CANAL ( 3 , Nout ( I - 1 ) + 1 , I - 1 ) + CANAL ( 4 , Nout ( I - 1 ) + 1 , I - 1 ) * Hu )* Hu ;
        end
        if Hu < Hd
            disp ( [ ' T= ' , num2str( T ) , '��I= ' , num2str( I ) , 'ʱ��Hu-Hd=' , num2str( Hu - Hd )  , 'բ��ˮλ����բǰˮλ�������쳣��'] ) ;
            pause ( 3 ) ; % ͣ��3��
        end
        GA ( I , 2 ) = feval(@Q2Greal , I , T , Hu , Hd , QP1 ( I , 2 ) , GA ( I , 1 ) , Gb ( I ) , Au , QP1 ( I , 2 ) , Gmax ( I ) , G , CDGATE ( I ) , DBGate , DB , E , DBWatLev , SilenceMode , CalcAcc ) ;
        GA ( I , 2 ) = round ( GA ( I , 2 ) * GateAccuracy ) / GateAccuracy ;
        if Hu < Hd
            QP1 ( I , 2 ) = QP1 ( I , 1 ) ;
        else
            QP1 ( I , 2 ) = feval(@GateDischarge , Hu , Hd , GA ( I , 2 ) , Gb ( I ) , Au , QP1 ( I , 2 ) , G , CDGATE ( I ) ) ; %cancel comment by GGH, 15.09.29  ����������ָ��ͻ���ֳ�ʼ״̬���ȣ���Ϊ����բ�ž��ȷ����Ŀ���������ʼ��������ƫ�ʹ�ĵ�����ʽ�������������ˮλ�������
        end
        DQGA ( I , T1 ) = QP1 ( I , 2 ) - QP1 ( I , 1 ) ;
        if DQGA ( I , T1 ) < GateDischargeAcc
            DQGA (I , T1 ) = 0 ;
        end
        % ����---------------------------------------------բ�ſ��ȼ���---------------------------------------����
        Qsum_up ( I ) = Qsum_up ( I ) + QP1 ( I , 2 ) * DT * NMinToHou ;          % �Ƚ�����������ֲο���
        % ����---------------------------------------------�Ǻ㶨������---------------------------------------����
        if I == Ncnl
            QPn ( I ) = Qt_d ( I ) ;
        else
            QPn ( I ) = QP1 ( I + 1 , 2 ) ;
        end

        SV_DQP1 = ( QP1 ( I , 2 ) - Q ( I , 1 ) ) / SV_N ;
        SV_DQPn = ( QPn ( I ) - Q( I , N ( I ) ) ) / SV_N ;
        for SV_I = 1 : SV_N  %��ϵͳ��һ������ʱ�䲽�������ηǺ㶨�����ɹ���
            if T < FailTime
                SV_Qout =feval(@ Qout , (T - DT + SV_I * SV_DT) , Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0) ;                           % ����Qout
            else
                SV_Qout =feval(@ AccidentResponse , Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, ( T - DT + SV_I * SV_DT) , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout );
            end
            %for K = 1 : Nout ( I )
            %    Zouth ( K , I ) = ZOUTH ( Nodtype ( K , I ) , Kouth ( K , I ) , Q ( I , sum ( N1 ( 1 : K , I ) ) ) ) ;
            %end
            if simtype == 1     %������
                if hiall ( T1 , I , N ( I ) ) > 0   %�ڱ��Ǻ�ȴ�����ʱ����ˮ��ѧ����
                    [Q( I , 1 : N ( I ) ) , Z( I  , 1 :N ( I ) )] = feval(@s_Preice , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , hiall ( T1 , I , 1 : N ( I ) ) , T , TTT ( I , : ) , IceDen , WatDen , IniRough , FinlRough ) ;
                else                                %�ڱ��Ǻ�ȵ�����ʱ����������Ǻ㶨������
                    [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) )] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
                end
            else                %��������
                [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) ] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) ,Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q( I , N(I) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
          
              %  [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) ] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) ,Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn, 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
            end
        end
        Hall = Z - TZ0 ;
        % ����---------------------------------------------�Ǻ㶨������---------------------------------------����
        % ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ļ��㲿��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����
        if simtype == 1
            % ����-------------------------------------------̫������-----------------------------------------����
            T_radi = floor ( T / ( NMinToHou * NHouToDay ) ) + 1 ;
            T_sunrise = ( T_radi - 1 ) * NHouToDay * NMinToHou * NMinToHou + 6 * NMinToHou * NMinToHou ;    % �ճ�ʱ��6��
            T_sunset = ( T_radi - 1 ) * NHouToDay * NMinToHou * NMinToHou + 18 * NMinToHou * NMinToHou ;    % ����ʱ��18��
            if T * NMinToHou > T_sunrise && T * NMinToHou <= T_sunset                                       % �ռ�12Сʱ������̫������ֲ�
                Radi = radi ( T_radi ) * ( 0.5 * pi / ( 12 * NMinToHou * NMinToHou ) ) * sin ( ( T * NMinToHou - T_sunrise ) * pi / ( T_sunset - T_sunrise ) ) ;
            else
                Radi = 0 ;
            end
            % ����-------------------------------------------̫������------------------------------------------����
            for j = 1 : N ( I ) %��ѭ������ӷ�ף�Ӧ�ð�������-������-�������ʽ��һ��ѭ������Ӽ��
                for II = 1 : Nout ( I ) + 1
                    if j - sum ( N1 ( 1 : II , I ) ) <= 0
                        a = II ;        % aΪj��ǰ�������ص��Ӷκţ����ڶ�ȡ�����������
                        break ;
                    end
                end
                D ( I , j ) = Zall ( T1 , I , j ) - TZ0 ( I , j ) - IceDen * hi ( I , j ) / WatDen ;        % ��Чˮ��
                B ( I , j ) = CANAL ( 3 , a , I ) + 2 * CANAL ( 4 , a , I ) * D ( I , j ) ;                 % ˮ����
                A ( I , j ) = ( CANAL ( 3 , a , I ) + CANAL ( 4 , a , I ) * D ( I , j ) ) * D ( I , j ) ;   % ��ˮ�������
                u ( I , j ) = Qall ( T1 , I , j ) / A ( I , j ) ;                                           % ����ƽ������
                DD ( I , j ) = A ( I , j ) / B ( I , j ) ;                                                  % ˮ�
                if u ( I , j ) > u1max ( I , j )
                    u1max ( I , j ) = u ( I , j ) ;
                end
                D1 ( I , j ) = Z ( I , j ) - TZ0 ( I , j ) - IceDen * hi ( I , j ) / WatDen ;               %D1��D��ʲô����Zall����һʱ�̵�ˮλ��
                B1 ( I , j ) = CANAL ( 3 , a , I ) + 2 * CANAL ( 4 , a , I ) * D1 ( I , j ) ;
                A1 ( I , j ) = ( CANAL ( 3 , a , I ) + CANAL ( 4 , a , I ) * D1 ( I , j ) ) * D1 ( I , j ) ;
                u1 ( I , j ) = Q ( I , j ) / A1 ( I , j ) ;
                Fr1 ( I , j ) = u1 ( I , j ) / sqrt ( G * D1 ( I , j ) ) ;   % �������
                DD1 ( I , j ) = A1 ( I , j ) / B1 ( I , j ) ;
                Ta ( T1 , I ) = Temp ( I,T1 * DT * NMinToHou ) ;              % ��ȡ����ֵ
                C( T1 , I ) = CloudAmount ( I,T1 * DT * NMinToHou ) ;         % ��ȡ����
                Va( T1 , I ) = WindVelo ( I,T1 * DT * NMinToHou ) ;           % ��ȡ����
                CAA ( I , : ) = CAAconst ;                                  % �趨�������ܶȷⶳ��׼
                d0 ( I , : ) = IniIceThick ;                                % �趨��������
                Ds1 = CANAL ( 5 , a , I ) / ( N1 ( a , I ) - 1 ) ;          % ���㲽��
                hiw = 1622 * u ( I , j ) ^ 0.8 * DD ( I , j ) ^ 0.2 ;       % ����ˮ���Ƚ���ϵ��
                %��������������������������仯
                fai ( I , j ) = HeatLossCalc ( Ta ( T1 , I ) , Radi , C( T1 , I ) , Va( T1 , I ) , hiall ( T1 , I , j ) , Twall ( T1 , I , j ) , Ca ( I , j ) , Ciall ( T1 , I , j ) ) ;  % ����ˮ�徻������ʧ
                if Twall ( T1 , I , j ) < 0
                    %ˮ�µ������ʱ��ˮ��ɱ��ų�������Ǳ�ȣ�ʹˮ�±�Ϊ���
                    Ci ( I , j ) = - ( WatDen * HeatCapa * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat ) + Ciall ( T1 , I , j ) ;
                    Tw ( I , j ) = 0 ;
                    if hiall ( T1 , I , j ) == 0        % �����������ܶ����ж��Ƿ��γɱ���
                        Ca ( I , j ) = Ci ( I , j ) * DD ( I , j ) / d0 ( I ) ;   %������Ũ�Ȼ�����������ܶ�    ������Ǻ�Ȳ��������أ�����Ҫ�����ʱ�������ܶ�ô��
                    end
                else
                    if j == 1
                        if I == 1
                            Tw ( I , j ) = UpTemp ;     % ��һ���ڵ��ˮ�µ���������ˮ�¶�
                        else                            %ÿ�����ص�һ���ڵ���¶ȼ��㣨��ȥ��һ�����أ�
                            if hiall ( T1 , I , j ) > 0 % �б���ʱ��ˮ�²��ܴ�����������Ӱ�죬������Ҫ��������
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - hiw * Twall ( T1 , I , j ) / ( WatDen * HeatCapa * DD ( I , j ) ) ;
                            elseif Ciall ( T1 , I , j ) > 0
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - DT * NMinToHou * ( fai ( I , j ) + hiw * Twall ( T1 , I , j ) ) / ( WatDen * HeatCapa * DD ( I , j ) ) ;
                            else
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - DT * NMinToHou * fai ( I , j ) / ( WatDen * HeatCapa * DD ( I , j ) ) ; %�¶ȱ仯=�ϸ��ڵ��¶�-����������ʧ
                            end
                        end
                    else
                        Tw ( I , j ) = Tww ( Tw ( I , j - 1 ) , Twall ( T1 , I , j - 1 ) , Twall ( T1 , I , j ) , fai ( I , j ) , u ( I , j - 1 ) , u ( I , j ) , DD ( I , j ) , hiall ( T1 , I , j ) , Ciall ( T1 , I , j ) , HeatCapa , hiw , DT , Ds1 , sum ( N1 ( 1 : a - 1 , I ) ) , j , DTIce ) ;
                    end
                    if hiall ( T1 , I , j ) == 0
                        if Ciall ( T1 , I , j ) < 0
                            Tw ( I , j ) = - Ciall ( T1 , I , j ) * IceDen * LatentHeat / ( WatDen * HeatCapa ) + Tw ( I , j ) ;
                            hi ( I , j ) = 0 ;
                            Ci ( I , j ) = 0 ;
                            Ca ( I , j ) = 0 ;
                        elseif Ciall ( T1 , I , j ) > 0
                            if j ~= sum ( N1 ( 1 : a , I ) )
                                if hiall ( T1 , I , j + 1 ) > 0    % �����ڱ���ǰԵ���ѻ�
                                    Ti = TTT ( I , j + 1 ) + 1 ;
                                    if j == 1
                                        if Ciaver ( I , j ) == 0
                                            Ciaver ( I , j ) = ( Ciall ( Ti , I , j ) * A ( I , j ) * Ds1 + ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                        else
                                            Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                        end
                                    elseif j == sum ( N1 ( 1 : a - 1 , I ) ) + 1
                                        if Nodtype ( a - 1 , I ) == 3
                                            if a == 2
                                                if Nodtype ( Nout ( I - 1 ) , I - 1 ) == 3
                                                    if  Ciaver ( I , j ) == 0
                                                        Ciaver ( I , j ) = ( Ciall ( Ti , I , j ) * A ( I , j ) * Ds1 + ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                                    else
                                                        Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                                    end
                                                else
                                                    if  Ciaver ( I , j ) == 0
                                                        Ciaver ( I , j ) = ( Ciall ( Ti , I , j ) * A ( I , j ) * Ds1 + ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                                    else
                                                        Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                                    end
                                                end
                                            else
                                                if  Ciaver ( I , j ) == 0
                                                    Ciaver ( I , j ) = ( Ciall ( Ti , I , j ) * A ( I , j ) * Ds1 + ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                                else
                                                    Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) ) / ( A ( I , j ) * Ds1 ) ;
                                                end
                                            end
                                        else
                                            if  Ciaver ( I , j ) == 0
                                                Ciaver ( I , j ) = ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                            else
                                                Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                            end
                                        end
                                    else
                                        if  Ciaver ( I , j ) == 0
                                            Ciaver ( I , j ) = ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                        else
                                            Ciaver ( I , j ) = Ciaver ( I , j ) + ( ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) * Ds1 * B ( I , j ) * DT * NMinToHou / ( IceDen * LatentHeat ) + Qall ( T1 , I , j ) * Ciall ( T1 , I , j ) * DT * NMinToHou ) / ( A ( I , j ) * Ds1 ) ;
                                        end
                                    end
                                else
                                    Ciaver ( I , j ) = Ciall ( T1 , I , j ) ;
                                end
                            else
                                Ciaver ( I , j ) = Ciall ( T1 , I , j ) ;
                            end
                            Ca ( I , j ) = Ciaver ( I , j ) * DD ( I , j ) / d0 ( I ) ;
                            if Ca ( I , j ) >= CAA ( I )               % �ﵽ�ⶳ��׼��ⶳ
                                if Ca ( I , j ) >= 1
                                    hi ( I , j ) = Ca( I , j ) * d0 ( I ) ;   % Vi(I,j)/(B1(I,j)*Ds1);
                                    Vii1 ( I , j ) = 0 ;               % ��϶��
                                else
                                    hi ( I , j ) = d0 ( I ) ;
                                    Vii1 ( I , j ) = B1 ( I , j ) * d0 ( I ) * ( 1 - Ca ( I , j ) ) ;  % ��϶�����
                                end
                                Ci ( I , j ) = 0 ;
                                TTT ( I , j ) = T1 ;
                            else                               % û�ﵽ�ⶳ��׼�������������
                                if a == 1
                                    if j == 1
                                        Ci ( I , j ) = DT * NMinToHou * ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat * DD ( I , j ) ) ;
                                    else
                                        Ciall1 = Ciall ( T1 , I , j ) ;
                                        u0 = u ( I , j ) / ( 1 + ( DTIce * NMinToHou / Ds1 ) * ( u ( I , j ) - u ( I , j - 1 ) ) ) ;
                                        for CiI = 1 : DT / DTIce
                                            Ci ( I , j ) = DTIce * NMinToHou * ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat * DD ( I , j ) ) + ( Ds1 - u0 * DTIce * NMinToHou ) * ( Ciall1 - Ciall ( T1 , I , j-1 ) ) / Ds1 + Ciall ( T1 , I , j - 1 ) ;
                                            Ciall1 = Ci ( I , j ) ;
                                        end
                                    end
                                else
                                    if I == 1
                                        if j == sum ( N1 ( 1 : a - 1 , I ) ) + 1    % �ڵ���һ������
                                            if Nodtype ( a - 1 , I ) == 3    % �ڵ�Ϊ������
                                                Ci ( I , j ) = DT * NMinToHou * ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat * DD ( I , j ) ) ;
                                            else
                                                Ci ( I , j ) = Ci ( I , j - 1 ) ;
                                            end
                                        else
                                            Ciall11 = Ciall ( T1 , I , j ) ;
                                            u0 = u ( I , j ) / ( 1 + ( DTIce * NMinToHou / Ds1 ) * ( u ( I , j ) - u ( I , j - 1 ) ) ) ;
                                            for CiI = 1 : DT / DTIce
                                                Ci ( I , j ) = DTIce * NMinToHou * ( fai ( I , j ) - hiw * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat * DD ( I , j ) ) + ( Ds1 - u0 * DTIce * NMinToHou ) * ( Ciall11 - Ciall ( T1 , I , j - 1 ) ) / Ds1 + Ciall ( T1 , I , j - 1 ) ;
                                                Ciall11=Ci(I,j);
                                            end
                                        end
                                    else
                                        Ci ( I , j ) = Cii ( fai ( I , j ) , DD ( I , j ) , u ( I , j - 1 ) , u ( I , j ) , Ci ( I , j - 1 ) , Ciall ( T1 , I , j - 1 ) , Ciall ( T1 , I , j ) , Twall ( T1 , I , j ) , Ds1 , a , sum ( N1 ( 1 : a - 1 , I ) ) , IceDen , LatentHeat , DT , I , j , hiw ,  Nodtype ( Nout ( I - 1 ) , I - 1 ) , Nodtype ( a - 1 , I ) , DTIce ) ;
                                    end
                                end
                                if j == sum ( N1 ( 1 : a , I ) )
                                    if a == Nout ( I ) + 1
                                        Ci ( I , j ) = Ci ( I , j ) + Ciall ( T1 , I , j ) ; % ����������ĩ���������ѻ�
                                    elseif Nodtype ( a , I ) == 3             % �����ڵ�����ǰ���������ѻ�
                                        if a == 1
                                            if I ~= 1
                                                if Nodtype ( Nout ( I - 1 ) , I - 1 ) ~= 3
                                                    Ci ( I , j ) = Ci ( I , j ) + Ciall ( T1 , I , j ) ;
                                                else
                                                    Ci ( I , j ) = Ci ( I , j ) ;
                                                end
                                            else
                                                Ci ( I , j ) = Ci ( I , j ) + Ciall ( T1 , I , j ) ;
                                            end
                                        else
                                            Ci ( I , j ) = Ci ( I , j ) + Ciall ( T1 , I , j ) ;
                                        end
                                    end
                                end
                            end
                        else
                            if j ~= 1
                                if j ~= sum ( N1 ( 1 : a-1 , I ) ) + 1
                                    if Ciall ( T1 , I , j - 1 ) > 0
                                        Ciall11 = Ciall ( T1 , I , j ) ;
                                        u0 = u ( I , j ) / ( 1 + ( DTIce * NMinToHou / Ds1 ) * ( u ( I , j ) - u ( I , j - 1 ) ) ) ;
                                        for CiI = 1 : DT / DTIce
                                            Ci ( I , j ) = ( Ds1 - u0 * DTIce * NMinToHou ) * ( Ciall11 - Ciall ( T1 , I , j-1 ) ) / Ds1 + Ciall ( T1 , I , j - 1 ) ;
                                            Ciall11 = Ci ( I , j ) ;
                                        end
                                    end
                                end
                            end
                        end
                    elseif hiall ( T1 , I , j ) > 0
                        [hi( I , j ), Vii1( I , j )] = hii ( fai ( I , j ) , Twall ( T1 , I , j ) , hiw , B1 ( I , j ) , IceDen , LatentHeat , Vii1 ( I , j ) , hiall ( T1 , I , j ) , Ca ( I , j ) , DT ) ;
                        if Ciall ( T1 , I , j ) > 0
                            ab = Ciall ( T1 , I , j ) * IceDen * LatentHeat / ( HeatCapa * 1000 ) ;
                            if Tw ( I , j ) >= ab
                                Tw ( I , j ) = Tw ( I , j ) - ab ;
                                Ci ( I , j ) = 0 ;
                            else
                                Tw ( I , j ) = 0 ;
                                Ci ( I , j ) = 0 ;
                                hi ( I , j ) = hi ( I , j ) + abs ( Tw ( I , j ) - ab ) * HeatCapa * 1000 * A ( I ,j ) / ( B ( I , j ) * IceDen * LatentHeat ) ;
                            end
                        elseif Ciall ( T1 , I , j ) < 0
                            Tw ( I , j ) = Tw ( I , j ) - Ciall ( T1 , I , j ) * IceDen * LatentHeat / ( 1000 * HeatCapa ) ;
                            Ci ( I , j ) = 0 ;
                        end
                    else
                        Tw ( I , j ) = abs ( hiall ( T1 , I , j ) ) * IceDen * LatentHeat / ( DD1 ( I , j ) * 1000 * HeatCapa ) + Tw ( I , j ) ;
                        hi ( I , j ) = 0 ;
                        Ci ( I , j ) = 0 ;
                        Ca ( I , j ) = 0 ;
                    end
                end
            end
        end
        for J = 1 : Nout ( I ) + 1
            Velo ( J , I ) = Vcmax ( CANAL , Q , Z , TZ0 , I , J , N1 ) ;   % ������������������
        end
        Velomax ( I ) = F_velomax( VeloENG (1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Limitlength , Nodtype2 ( 1 : Nout ( I ) + 1 , I ) , Velo ( 1 : Nout ( I ) + 1 , I ) , CANAL ( 5 ,1 : Nout ( I ) + 1 , I ) );
        [V( I ) , Delay( I )] = Vol ( Nout( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 : N ( I ) ) ) ;% ������ˮ�����
        % ����-------------------------------------��¼��һʱ�̵ĸ�����-----------------------------------------����
        IAQ1 ( I ) = IAQ1 ( I ) + abs ( QP1 ( I , 2 ) - QP1 ( I , 1 ) ) ;
        IAW1 ( I ) = IAW1 ( I ) + abs ( GA ( I , 2 ) - GA ( I , 1 ) ) ;
        E ( I , 1 ) = E ( I , 2 ) ;
        EC ( I , 1 ) = EC ( I , 2 ) ;
        Ef ( I , 1 ) = Ef ( I , 2 ) ;
        ECf ( I , 1 ) = ECf ( I , 2 ) ;
        GA ( I , 1 ) = GA ( I , 2 ) ;
        QP1 ( I , 1 ) = QP1 ( I , 2 ) ;
        % ����-------------------------------------��¼��һʱ�̵ĸ�����-----------------------------------------����
        T2 = T1 / 1 + 1 ;
        if T2 == fix ( T2 )    % ����ʱ�̼�¼����
            Tx ( T2 ) = T / NMinToHou ;
            Q_U( I , T2 ) = QP1 ( I , 2 ) ;
            Qoutsum ( I , T2 ) = sum ( SV_Qout ( : , I ) ) ;
            QdownT(1)=Qdown ( 0 ,Tdown_start,Qdown_start,Tdown_end,Qdown_end ) ;
            QdownT ( T2 ) = Qdown ( T ,Tdown_start,Qdown_start,Tdown_end,Qdown_end ) ;
            GA1 ( I , T2 ) = GA ( I , 2 ) ;
            YTorig_UP ( I , T2 ) = YTorig_U ( I ) ;
            YTorig_DOWN ( I , T2 ) = YTorig_D ( I ) ;
            Velocity ( I , T2 ) = Velomax ( I ) ;
            Velocity1N( I , T2 ) =V0_u ( I ,1 );
            Z_U ( I , T2 ) = Z ( I , 1 ) ;                                      %���������ζ���ˮλ
            Z_D ( I , T2 ) = Z ( I , N ( I ) ) ;                                %���������ζ���ˮλ
            Z_M ( I , T2 ) = 0.5 * Z ( I , 1 ) + 0.5 * Z ( I , N ( I ) ) ;
            Volu ( I , T2 ) = V ( I ) ;
            Delayl ( I , T2 ) = Delay ( I ) ;
            Ddz_U ( I , T2 ) = ZJU ( I ) - Z ( I , 1 ) ;                         % ����У��ˮλ��ʵ��ˮλ��
            Ddz_D ( I , T2 ) = ZJD ( I ) - Z ( I , N ( I ) ) ;                   % ����У��ˮλ��ʵ��ˮλ��
            Ddz_M ( I , T2 ) = Ddz_U ( I , T2 ) * 0.5 + Ddz_D ( I , T2 ) * 0.5;  % ����У��ˮλ��ʵ��ˮλ��
            EH_U ( I , T2 ) = Z ( I , 1 ) - YTorig_U ( I ) ;                     % ����������ˮλ���
            EH_D ( I , T2 ) = Z ( I, N ( I ) ) - YTorig_D ( I ) ;                % ����������ˮλ���
            EH ( I , T2 ) = 0.5 * EH_U ( I , T2 ) + EH_D ( I , T2 ) * 0.5 ;      % ����������ˮλ���
            V1 ( I , T2) = V ( I ) ;
            YT1 ( I , T2 ) = YT ( I ) ;
            YF1 ( I , T2 ) = YF ( I ) ;
            if Tx ( T2 ) > T_MDS1
                DT2 = floor ( T_MDS1 * NMinToHou / DT ) ;
                T3 = T2 - DT2 ;
                DV_U1 ( I , T2 ) = Z_U ( I , T2 ) - Z_U ( I , T3 ) ;
                DV_D1 ( I , T2 ) = Z_D ( I , T2 ) - Z_D ( I , T3 ) ;
                DV_M1 ( I , T2 ) = Z_M ( I , T2 ) - Z_M ( I, T3 ) ;
                if SilenceMode == 0
                    if abs ( DV_U1 ( I , T2 ) ) > MDS1
                        disp ( ['����',num2str( I ),'����1hˮλ������ޣ�'] )
                    elseif abs ( DV_D1 ( I , T2 ) ) > MDS1
                        disp ( ['����',num2str( I ),'����1hˮλ������ޣ�'] )
                    elseif abs ( DV_M1 ( I , T2 ) ) > MDS1
                        disp ( ['����',num2str( I ),'����1hˮλ������ޣ�'] )
                    end
                end
            else
                DV_U1 ( I , T2 ) = Z_U ( I , T2 ) - Z_U ( I , 1 ) ;
                DV_D1 ( I , T2 ) = Z_D ( I , T2 ) - Z_D ( I , 1 ) ;
                DV_M1 ( I , T2 ) = Z_M ( I , T2 ) - Z_M ( I , 1 ) ;
            end
            if Tx ( T2 ) > T_MDS2
                DT3 = floor ( T_MDS2 * NMinToHou / DT ) ;
                T3 = T2 - DT3 ;
                DV_U2 ( I , T2 ) = Z_U ( I , T2 ) - Z_U ( I , T3 ) ;
                DV_D2 ( I , T2 ) = Z_D ( I , T2 ) - Z_D ( I , T3 ) ;
                DV_M2 ( I , T2 ) = Z_M ( I , T2 ) - Z_M ( I , T3 ) ;
                if SilenceMode == 0
                    if abs ( DV_U2 ( I , T2 ) ) > MDS2
                        disp ( ['����' , num2str( I ) , '����24hˮλ������ޣ�'] )
                    elseif abs ( DV_D2 ( I , T2 ) ) > MDS2
                        disp ( ['����',num2str( I ),'����24hˮλ������ޣ�'] )
                    elseif abs ( DV_M2 ( I , T2 ) ) > MDS2
                        disp ( ['����' , num2str( I ) , '����24hˮλ������ޣ�'] )
                    end
                end
            else
                DV_U2 ( I , T2 ) = Z_U ( I , T2 ) - Z_U ( I , 1 ) ;
                DV_D2 ( I , T2 ) = Z_D ( I , T2 ) - Z_D ( I , 1 ) ;
                DV_M2 ( I , T2 ) = Z_M ( I , T2 ) - Z_M ( I , 1 ) ;
            end
        end
        for K = Ncnl : -1 : 1                         % ���������������ۼƹ�ˮ��֮��
            if I == Ncnl                          % 2010-10-20��ѭ��������I��ΪK
                DQsum_up ( K ) = Qsum_up ( K ) - Qsum_down ;
            else
                DQsum_up ( K ) = Qsum_up ( K ) - Qsum_up ( K + 1 ) ;
            end
        end
    end
    Qout1 = Qout0 ;                                % ��¼��һ�ε�Qout0
    Q01_d = Qt_d ( Ncnl ) ;                           % ��¼��һʱ��������ˮ�����仯
    Zall ( T2 , : , : ) = Z ;                             % ��¼��ͬʱ����ϸ����ε�ˮλ������Zall(5,3,11)��ʾ��5�����ڵ��������ε�11���ڵ��ϵ�ˮλ
    Qall ( T2 , : , : ) = Q ;
    %����
    Twall ( T2 , : , : ) = Tw ;
    hiall ( T2 , : , : ) = hi ;
    Ciall ( T2 , : , : ) = Ci ;
    %����
    CurrentTimeSec = mod ( T , 360 );
%    [ LengthIce(1:Ncnl,T/DT+1) , PercLength(1:Ncnl,T/DT+1) ] = ShowIceCoverPerc ( LengthIce(:,:) , PercLength(:,:) , hiall (: , : , : ) , N1(:,:) , CANAL(:,:,:) , DT , T , Len_Canal ( : ) , Ncnl , SilenceMode , N(:) , CurrentTimeSec , ShowIceTime );
end

% Part VI : Data recording and output
for I = Ncnl : -1 : 1
    HiceD1 ( I , 1 : N ( I ) ) = zeros ( 1 ,  N ( I ) ) ;
    HiceD3 ( I , 1 : N ( I ) ) = zeros ( 1 ,  N ( I ) ) ;
    HiceD5 ( I , 1 : N ( I ) ) = zeros ( 1 ,  N ( I ) ) ;
    Jice  = 0;
    LengthPool ( I ) = sum ( CANAL ( 5 , : , I ) ) ;
    for K = 1 : N ( I )
        if Tmax <= 24
            HiceD1( I , K  ) = hiall ( 1 * NHouToDay * NMinToHou / DT , I , K ) ;
        elseif Tmax <= 72
            HiceD1( I , K  ) = hiall ( 1 * NHouToDay * NMinToHou / DT , I , K ) ;
            HiceD3( I , K  ) = hiall ( 3 * NHouToDay * NMinToHou / DT , I , K ) ;
        else
            HiceD1( I , K  ) = hiall ( 1 * NHouToDay * NMinToHou / DT , I , K ) ;
            HiceD3( I , K  ) = hiall ( 3 * NHouToDay * NMinToHou / DT , I , K ) ;
            HiceD5( I , K  ) = hiall ( 5 * NHouToDay * NMinToHou / DT , I , K ) ;
        end
        for Jice  = 1 : Nout(I) + 1 ;
            if HiceD1 ( I , K ) > 0 % ������1��ʱ�̣�01��00:00:00�������Χ
                if K > sum ( N1 ( 1 : Jice ,I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_1d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_1d ( I ) = 0;
            end
            if HiceD3 ( I , K ) > 0 % ������1��ʱ�̣�03��00:00:00�������Χ
                if K > sum ( N1 ( 1 : Jice , I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_3d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_3d ( I ) = 0;
            end
            if HiceD5 ( I , K ) > 0 % ������1��ʱ�̣�05��00:00:00�������Χ
                if K > sum ( N1 ( 1 : Jice  , I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_5d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_5d ( I ) = 0;
            end
        end
    end
end
Tw_u = Twall ( : , : , 1 ) ;    % ����բ��ˮ��
hi_u = hiall ( : , : , 1 ) ;    % ����բ��ˮ��
[ XTw_u, YTw_u ] = size ( Tw_u );
Tw_d = zeros( XTw_u, YTw_u );
for I = Ncnl : -1 : 1
    Tw_d = Twall ( : , : , N(I) ) ; % ��������ˮ��
    hi_d = hiall ( : , : , N(I) ) ; % �������α���
    IAQ ( I ) = IAQ1 ( I ) - abs ( Q_U ( I , 1 ) - Q_U ( I , T2 ) ) ;
    IAW ( I ) = IAW1 ( I ) - abs ( GA1 ( I , 1 ) - GA1 ( I , T2 ) ) ;
end
IAE = IAE1 * DT / NMinToHou / Tmaxm ;
Jug = Jug1 * DT / Tmaxm ;
Jug_sum = sum ( Jug1 ) * DT / Tmaxm ;             % Ѱ�Ų���
disp ( ['��ǰʱ��Ϊ______' , datestr( now, 31 ) ] ) ;
disp ('=================�������===================') ;
save SIM_RES
Tsim_S = toc;
if SilenceMode == 0
    disp ('============������������Ŀ��Ŀ¼�µ�SIM_RES.mat�ļ���============') ;
    disp ( ['�����ʱ' , num2str( Tsim_S / NMinToHou ) , '����'] );
end
Output;
if SilenceMode == 0
    disp ('===բ�ſ��ȡ���բ���������Ƶ�ˮλ��������١����Ǻ�ȡ�����ˮ�¡����Ƿ�Χ�ֱ���txt��ʽ������Ŀ��Ŀ¼��===') ;
    disp ('===բ�ſ��ȡ���բ���������Ƶ�ˮλ��������١����Ǻ�ȡ�����ˮ�¡����Ƿ�Χ�ֱ���xls��ʽ������Ŀ��Ŀ¼��===') ;
    pause ( 5 ) ; % ͣ��5��
end
% ����****************************************��������ͼ****************************************����
if SilenceMode == 0
    disp ( '�Ƿ���Ҫ�Է������ݽ���ѡ���Ի�ͼ 1 = �� ����� = ��' )
    Choice = input (' ����������ѡ��:','s');
    switch Choice
        case '1'
            disp ('======��������SIM_RES.mat�ļ��е����ݿ�ѡ���Խ��л�ͼ========');
            PPNi
        otherwise
            disp ('=============�����ͼ��������������������===============');
    end
end
%����****************************************��������ͼ****************************************����
