% Simulation and Control of multi canal systems
% Originally written by Guanghua Guan, Wei Cui and Jie Fan inspired by Changde Wang  @ 2006
% Upgraded By Xiong Yao, Zhiliang Ding, Mengkai Liu, Guoqiang Liu and Yibo Yan
clear;
clc;
tic

% Part I : Data input
% 以上****************************************读入系统设定常量****************************************以上
NMinToHou = 60 ;                        % 时间进制，1小时=60分钟，1分钟=60秒
NHouToDay = 24 ;                        % 1天=24小时
Limitlength = 200 ;                     % 定义渠池断面计算冰期流速时最短长度
SPRaise = 0.1 ;                         %冰期控制目标水位抬高值
ShowIceTime = 360;                      %冰情输出的时间间隔
CalcAcc = 1e-14;                         %计算精度，原值 1e-15
% % 以上**************************************读入系统设定常量**************************************以上
% disp ( '选择仿真数据读入方式 1 = txt; 2 = Excel; 任意键 = database' )
% Choice = input (' 请作出您的选择:','s');
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
    disp ( '输入数据是否有改动?' )
    cprintf('red','如有改动请按“y” 并回车重新载入数据；');
    cprintf('Comments','直接按回车键开始仿真；');
    Choice = input ('请您作出选择:','s');
    switch Choice
        case 'y'
            disp ('重新由文件读入数据！');
%             delete InputSIM.mat
            dataEXCEL ( NMinToHou, NHouToDay, SPRaise ) ;
            load InputSIM
        otherwise
            disp ( '采用当前目录中的输入数据开始仿真！' ) ;
    end
catch
    disp ( '当前目录中没有输入数据，从文件读入！' ) ;
    dataEXCEL ( NMinToHou, NHouToDay, SPRaise ) ;
    load InputSIM
end
% 以下***************************************人机交互信息*****************************************以下

if IfTermate==1;
    return
end
% 以上***************************************人机交互信息*****************************************以上

% Part II : defination of variable and matrix
Zouthini = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;        %建筑物局部水头损失（均匀流）
Zouthqtt = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;        %建筑物局部水头损失（t时刻的水头损失）
Zouthunst = zeros ( max ( Nout ( : ) )+ 1  , Ncnl ) ;       %建筑物局部水头损失（初始时刻水头损失）
Len_Canal = zeros ( Ncnl , 1 ) ;                            %渠池长度
V_wave = zeros ( Ncnl , 1 ) ;                               %水波速度
T_wave = zeros ( Ncnl , 1 ) ;                               %水波传播时间
TZ0 = zeros ( Ncnl , max ( sum (  N1 ( : , 1 : Ncnl ) ) ) ) ;   %渠道各计算断面渠底高程
Z0 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;     %渠道各断面设计水位
Qt_d = zeros ( Ncnl , 1 ) ;                                 % 各渠段计划下游流量过程 (Q0)
Qsum_down = 0 ;                                             % 记录下游需水总量
QP1 = zeros ( Ncnl , 2 ) ;                                  % 各渠段上游流量边界条件,即渠段上游闸门过闸流量
QPn = zeros ( Ncnl , 1 ) ;                                  % 各渠段下游流量边界条件，即渠段下游闸门过闸流量
GA = zeros ( Ncnl , 2 ) ;                                   % 闸门开度
N = zeros ( Ncnl , 1 ) ;                                    %一个渠段内的计算断面总数
Q = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;  %渠道各断面流量
Qtt = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;%渠道各断面流量
Z = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;  %渠道各断面水位
Z1 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;
Z2 = zeros ( Ncnl , max ( sum ( N1 ( : , 1 : Ncnl ) ) ) ) ;
Vini = zeros ( Ncnl , max ( Nout( : ) ) + 1 );              %渠道各断面初始平均流速%%%%%%%%%%%%%%%%%%%%%%自己加的
Vqtt = zeros ( Ncnl , max ( Nout( : ) ) + 1 ) ;             %渠道各断面下游常水位运行时平均流速%%%%%%%%%%%%%%%%%%%%%%%自己加的
Vunsteady = zeros ( Ncnl , max ( Nout( : ) ) + 1 ) ;        %非恒定流隔断面流速的计算
Hd0 = zeros ( Ncnl , 1 ) ;                                  %下游常水位运行下游端初始“水深”
V = zeros ( Ncnl , 1 ) ;                                    %各渠段水体体积
Delay = zeros ( Ncnl , 1 ) ;                                %渠池滞后时间
Velo = zeros ( (max ( Nout)+ 1 ) , Ncnl ) ;                 %每个渠池内各断面平均流速
Velomax = zeros ( Ncnl , 1 ) ;                              %各渠池最大流速
VeloENG = zeros (  max ( Nout ( : ) )+ 1 , Ncnl ) ;         %工程流速
% 以下------------------------------------记录各渠段水体体积---------------------------------------以下
Tmaxm = Tmax * NMinToHou ;                                  % 仿真总时间(分钟)
Volu = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                    %不同时间点各渠池水体体积
Delayl = zeros ( Ncnl , Tmaxm / DT + 1 ) ;
Qtal = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                    %记录前馈流量
% 以上------------------------------------记录各渠段水体体积---------------------------------------以上
% 以下%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%冰期程序所增部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以下
for I = Ncnl : -1 : 1
    N ( I ) = sum ( N1 ( : , I ) ) ;                        %一个渠段内的计算断面数
    % 冰期有关变量的声明,冰期参数，需考虑变量名称是否需要调整，工作量较大。
    hi ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % 冰盖厚度
    %Tw ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % 断面平均水温
    Tws ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;   	% 冰盖上表面温度
    Ca ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % 面流冰密度
    Ci ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % 流冰浓度
    TTT ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % 封冻时刻，需要计算的次数表达
    Ciaver ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;    % 空间步长上的平均流冰浓度
    CAA ( I ) = zeros ( 1 ) ;                               % 面流冰密度封冻标准
    d0 ( I ) = zeros ( 1 ) ;                                % 设定的流冰层厚度
    Vi ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;        % 冰量
    u1max ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;     % u断面平均流速
    % 冰情初始值
    % !!冰情初始条件有问题，如果温度低于零度时冰盖厚度、流冰密度等参数不一定全是零。
    %Tw ( I , 1 : N ( I ) ) = 1 ;                           % 水温初始条件
    hi ( I , 1 : N ( I ) ) = 0 ;                            % 冰盖厚度
    Ci ( I , 1 : N ( I ) ) = 0 ;                            % 流冰浓度
    Ca ( I , 1 : N ( I ) ) = 0 ;                            % 面流冰密度
    Vii1 ( I , 1 : N ( I ) ) = 0 ;                          % 冰盖空隙体积
end
% 以下------------------------------------计算冰盖封冻率变量初始化---------------------------------------以下
PercLength = zeros ( Ncnl , Tmaxm / DT + 1 );
LengthIce = zeros ( Ncnl , Tmaxm / DT + 1);
% 以上------------------------------------计算冰盖封冻率变量初始化---------------------------------------以上
% 以上%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%冰期程序所增部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上
Qt_u = zeros ( Ncnl , 1 ) ;                                 % 各渠段计划上游流量过程 (Q1)
YT = zeros ( Ncnl , 1 ) ;                                   % 控制点恒定流目标水位
Qsum_up = zeros ( Ncnl + 1 , 1 ) ;                          % 各渠段渠首流量累计水量(Qsum)
Huall = zeros ( Tmaxm / DT , Ncnl ) ;                       % 不同时间点上的闸前水深
Hdall = zeros ( Tmaxm / DT , Ncnl ) ;                       % 不同时间点上的闸后水深
Hdiff = zeros ( Tmaxm / DT , Ncnl ) ;                       % 不同时间点上的闸前闸后水深之差
Q0_d_0= zeros ( Ncnl , 1 ) ;
V0_u = zeros ( Ncnl , 1 ) ;
Tx = zeros ( Tmaxm / DT + 1 , 1 ) ;                           % 时间数据横坐标
Q_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 各渠段上游流量输出
Qoutsum = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                   % 各渠段区间取水口流量输出
Z_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 各渠段上游水位输出
Z_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 各渠段下游水位输出
Z_M = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 各渠段中点计算水位输出
GA1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 各渠段上游闸门开度输出
V1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                        % 各渠段水体体积输出
YT1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 控制点恒定流目标水位输出
YF1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                       % 控制点实际水位输出
YTorig_U = zeros ( Ncnl , 1 ) ;                               % 恒定流上游目标水位
YTorig_D = zeros ( Ncnl , 1 ) ;                               % 恒定流下游目标水位
YTorig_UP = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                 % 恒定流上游目标水位输出
YTorig_DOWN = zeros ( Ncnl , Tmaxm / DT + 1 ) ;               % 恒定流下游目标水位输出
Velocity = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                  % 各渠池最大流速随仿真的变化过程
 Velocity1N = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                  % 各渠池渠首过闸处的流速大小变化情况
Vttl = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                      % 蓄量补偿前馈目标蓄量计算值输出
E = zeros ( Ncnl , 2 ) ;                            % 控制点水位误差
Ef = zeros ( Ncnl , 2 ) ;                           % 滤波后的控制点水位误差
EC = zeros ( Ncnl , 2 ) ;                           % 控制点水位误差变化率
ECf = zeros ( Ncnl , 2 ) ;                          % 滤波后的控制点水位误差变化率
Jug1 = zeros ( Ncnl , 1 ) ;                         % 各渠段水位误差指标（ISE 二次性能指标）中间变量
IAE1 = zeros ( Ncnl , 1 ) ;                         % 计算各渠段水位指标中间变量
IAQ1 = zeros ( Ncnl , 1 ) ;                         % 计算各渠段流量指标中间变量
IAW1 = zeros ( Ncnl , 1 ) ;                         % 计算各渠段闸门指标中间变量
IAQ = zeros ( Ncnl , 1 ) ;                          % 计算各渠段流量指标
IAW = zeros ( Ncnl , 1 ) ;                          % 计算各渠段闸门指标
DQGA_sum = zeros ( Ncnl , 1 ) ;                     % 各个渠池中的上游过闸流量变化累积（用于判断系统是否达到稳定）
QdownT = zeros ( Tmaxm / DT + 1 , 1 ) ;             % T时刻的系统最下游端取水流量
YF = zeros ( Ncnl , 1 ) ;                             % 控制点非恒定流目标水位
Vtt = zeros ( Ncnl , 1 ) ;                            % 蓄量补偿前馈目标蓄量计算
DQGA = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % 除数因为DT，采用SV_DT为了增大容量，方便可能做得修改
Ddz_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段上游校实差
Ddz_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段下游校实差
Ddz_M = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段中游校实差
EH_U = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % 各渠段上游水位误差
EH_D = zeros ( Ncnl , Tmaxm / DT + 1 ) ;              % 各渠段下游水位误差
EH = zeros ( Ncnl , Tmaxm / DT + 1 ) ;                % 各渠段中游水位误差
DV_U1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段1h上游水位降
DV_D1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段1h下游水位降
DV_M1 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段1h中游水位降
DV_U2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段24h上游水位降
DV_D2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段24h下游水位降
DV_M2 = zeros ( Ncnl , Tmaxm / DT + 1 ) ;             % 各渠段24h中游水位降
DQsum_up = zeros ( Ncnl + 1 , 1 ) ;                   % 各渠段上下游过水量之差
LengthPool = zeros ( Ncnl , 1 ) ;
len_1d = zeros ( Ncnl , 1 ) ;
len_3d = zeros ( Ncnl , 1 ) ;
len_5d = zeros ( Ncnl , 1 ) ;
DQff= zeros ( Ncnl , Tmaxm / DT + 1 ) ;
DQfb= zeros ( Ncnl , Tmaxm / DT + 1 ) ;
Dvttl= zeros(Ncnl,Tmaxm / DT);
% Part III : Initialize of variables
SV_N = DT/SV_DT ;           %每次控制动作时间内的非恒定流计算时间步数
% 以下***************************************渠道参数计算*****************************************以下
% 以下--------------------------计算各个渠池中水波由下游传播到上游所需时间--------------------------以下
for I = 1 : Ncnl
    if I==51
        I=I;
    end
    Len_Canal ( I ) =  sum( CANAL ( 5 , 1 : Nout ( I ) + 1 , I ) ) ;          % 渠池长度
    V_wave ( I ) = ( G * Hd01 ( I ) ) ^ 0.5 * NMinToHou ;                        % 单位：m/min  %"Hd01"下游常水位运行下游端初始“水深”
    T_wave ( I ) = Len_Canal ( I ) / V_wave ( I ) ;   %波传播时间
    T_wave ( I ) = round ( T_wave ( I ) / DT ) * DT ;                    % 取整
    % 以上--------------------------计算各个渠池中水波由下游传播到上游所需时间--------------------------以上
    % 以下--------------------------------------计算渠底高程TZ0---------------------------------------以下
    zt1 = z0_u ( I ) ;                                      % %各渠段渠底起点高程
    Jst = 1 ;                                                % J1:赋值起点，J2:赋值终点
    for K = 1 : Nout ( I ) + 1                                    % K:一个渠段内的分段数 Nout:各渠段节点个数
        Jend = Jst + N1 ( K , I ) - 1 ;                                    % N1： 各节点之间计算断面的个数
        zt2 = zt1 - CANAL ( 2 , K , I ) * CANAL ( 5 , K , I ) ;
        TZ0 ( I , Jst : Jend ) = zt1 : ( zt2 - zt1 ) / ( Jend - Jst ) : zt2  ; %TZ0:渠道各计算断面渠底高程；假设计算断面均匀分布
        Z0 ( I , Jst : Jend ) = TZ0 ( I , Jst : Jend ) + Hd01 ( I ) ;
        if  K~= Nout ( I ) + 1
            zt1 = zt2 - Zoutd ( K , I ) ; %zt1:计算渠底高程时子渠段起点高程;zt2:计算渠底高程时子渠段终点高程  % Zoutd：各取水口(建筑物)初始渠底降落
            Jst = Jend + 1 ;
        end
    end
    % 以上--------------------------------------计算渠底高程TZ0---------------------------------------以上
end
% 以上***************************************渠道参数计算*****************************************以上
% 以下--------------------------取水口变化始终时间及流量初始化--------------------------以下
if Qdown_start == Qdown_end                     %如果系统末端流量无变化，流量变化时间仅考虑取水口变化 ；初始时刻分水流量=末了时刻分水流量
    TS = min ( Tstart );
    TE = max ( Tend );
else                                            %如果系统末端流量有变化，流量变化时间应考虑取水口和末端时间的迭加
    TS = min ( min ( Tstart ) , Tdown_start );  %流量变化开始时刻
    TE = max ( max ( Tend ) , Tdown_end ) ;     %流量变化结束时刻
end
for I = Ncnl : -1 : 1                           %各段下游初始流量
    if I == Ncnl
        Qt_d ( I ) = Qdown_start ;                %Qdown_start：下游初始流量     Qt_d ( I )：第I个渠段的下游初始流量
    else
        Qt_d ( I ) = Qt_d ( I + 1 ) + sum ( Qout0 ( : , I + 1 ) );         %Qout0： 各取水口初始流量
    end
end
Qout1 = Qout0 ;                         % 记录上一次的Qout0
Q01_d = Qt_d ( Ncnl ) ;                 % 记录上一时刻下游需水流量变化
% 以上--------------------------取水口变化始终时间及流量初始化--------------------------以上
% 以上------------------------------读入工况文件，计算初始闸门流量----------------------------------以上
NodtypeGate ( 1 : Ncnl ) = 5 ;    % Nodtype2是什么意思？
Nodtype2 = vertcat ( Nodtype , NodtypeGate );
% for I = 1 : Ncnl
%     Nodtype2 ( Nout ( I ) + 1 , I ) = 5 ;         % 闸门节点，第Nout(I)+1个子渠段流速计算的需要
% end

% Part IV : Initial state calculation of canal system
% 以下---------------------------------------流量、水位初始条件计算--------------------------------以下
for I = 1 : Ncnl
    % 以下---------------------------------------流量初始条件-------------------------------------以下
    Q ( I , 1 : N ( I ) ) =feval(@ Ini_Q , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , Qt_d ( I ) ) ; %渠道各断面计算流量=最下游流量+所有流经的取水口流量
    % 以上---------------------------------------流量初始条件-------------------------------------以上
    % 以下---------------------------------------水位初始条件-------------------------------------以下
    % 以下-------------------------------------局部水头损失计算-----------------------------------以下
    Vini ( I , 1 : Nout ( I ) + 1 ) = feval(@VSegmentavg , CANAL , Q , Z0 , TZ0 , I , N1 , Nout( I ) ) ;                           %Z0：均匀流水面线
    Zouthini (  1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vini , Nout ( I ) , I );
    % 以上-------------------------------------局部水头损失计算-----------------------------------以上
    % 以下--------------------------------下游常水位运行下游水深计算-------------------------------以下
    if Wt == 0
        Hd0 ( I ) = Hd01 ( I ) ;
    end
    % 以上--------------------------------下游常水位运行下游水深计算-------------------------------以上
    % 以下---------------等体积运行，通过设计水深求各渠段体积，反算各渠段初始下游水深----------------以下
    if Wt == 0.5
%         for J = 1 : N ( I )
%             Z ( I , J ) = Hd01 ( I ) + TZ0 ( I , J ) ;
%         end
        Z ( I , : ) = Hd01 ( I ) *ones (  1 , N ( I ) ) + TZ0 ( I , : );
        [ V( I ) , Delay( I )] = feval(@Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I )+1, I ) , N1 ( 1 : Nout ( I )+1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) ,Q(I,1:N(I))) ;      % 由于建筑物等的影响，和实际的真实体积有误差
        Hd0 ( I ) =feval(@ V2Hd , V ( I ) , CANAL( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthini ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , [ 2 8 ] ) ;   % 等体积运行各渠段下游端初始“水深”
    end
    % 以上---------------等体积运行，通过设计水深求各渠段体积，反算各渠段初始下游水深----------------以上
    % 以下---------------递代法求初始时刻水面线---------------------------------------------------以下
    Z1 ( I , 1 : N ( I ) )=feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthini ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I ) , I ) ;
    Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z1 , TZ0 , I , N1 , Nout ( I ) ) ;               % 计算各子渠段最下游流速
    Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    Z2 ( I , 1 : N ( I ) ) =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I )  , I ) ;     %15-9-26改
    Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z2 , TZ0 , I , N1 , Nout ( I ) ) ;               % 计算各子渠段最下游流速
    Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    if max ( Z1 - Z2 ) > CalcAcc;
        Z1 = Z2 ;
        Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z1 , TZ0 , I , N1 , Nout ( I ) ) ;               % 计算各子渠段最下游流速
        Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
        Z2 ( I , 1 : N ( I ) ) =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hd0 ( I ) , I  ) ;     %15-9-26改
        Vunsteady ( I , 1 : ( Nout ( I ) + 1 ) )  =feval(@ VSegmentavg , CANAL , Q , Z2 , TZ0 , I , N1 , Nout ( I ) ) ;               % 计算各子渠段最下游流速
        Zouthunst ( 1 : ( Nout ( I ) + 1 ) , I ) = feval(@HeadLoss , Nodtype2 , Kouth  ,Vunsteady , Nout ( I ) , I );
    else
        Z = Z2;
    end
    Z11 = Z ( 1 , 1 : N ( 1 ) );
    % 以上---------------------------------------水位初始条件-------------------------------------以上
    % 以下---------------------------------------渠池流速计算-------------------------------------以下
    for J = 1 : Nout ( I ) + 1
        Velo ( J , I ) = feval ( @Vcmax , CANAL , Q , Z , TZ0 , I , J , N1 ) ;               % 计算各子渠段最大流速
    end
    Velomax ( I ) = feval ( @F_velomax , VeloENG (1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Limitlength , Nodtype2 ( 1 : Nout ( I ) + 1 , I ) , Velo ( 1 : Nout ( I ) + 1 , I ) , CANAL ( 5 ,1 : Nout ( I ) + 1 , I ) );
    % 以上---------------------------------------渠池流速计算-------------------------------------以上
end
% 以上---------------------------------------流量、水位初始条件计算--------------------------------以上

% 以下--------------------------------------------初始闸门开度------------------------------------以下
for I = 1 : Ncnl
    if I==30
        I=I;
    end
    if I == 1
        Q0_d_0 ( I , 1 ) = Qt_d ( I ) + sum ( Qout0 ( : , I ) );
        V0_u ( I , 1 ) =Q0_d_0 ( I , 1 )/ ( CANAL ( 3 ,1 ,I ) + CANAL ( 4 ,1 ,I ) * ( Z ( I ,1 ) - TZ0 ( I , 1 ) ) )/( Z ( I , 1 ) - TZ0 ( I , 1 ) );
%         Hu = Hu1 - Kouth ( 1 , 1 ) * V0_u( I , 1 )^ 2 /2/G ;%第一渠段（渠池）上游水深            %Kouth： 各节点水头损失系数
        Hu = Hu1;
        Hd = Z ( I , 1 ) - z0_u ( I ) ;%+ Kouth ( Nout ( I -1) + 1 , I ) *  V0_u ( I ,1 )^ 2/2/G ;%下游水深   Hu1：最上游水库水深
    else
        Q0_d_0 ( I , 1 ) = Qt_d ( I - 1 ) ;
        V0_u ( I , 1 ) =feval(@ VPoolavg , CANAL , Q0_d_0 , Z , TZ0 , I-1 , N1 , Nout ( I-1 ) ) ;               %15-9-22改
        Hu = Z ( I - 1 , N ( I - 1 ) ) - z0_d ( I - 1 ) ;
        Hd = Z ( I , 1 ) - z0_d ( I - 1 ) + Kouth ( Nout ( I -1) + 1 , I-1 ) *  V0_u ( I ,1 )^ 2/2/G ;          %z0_d ( I - 1 ) 上一个渠段终点的渠底高程  Z ( I , 1 )：第一个计算断面的水面高程
    end
    Huall ( 1 , I ) = Hu ;
    Hdall ( 1 , I ) = Hd ;
    Hdiff ( 1 , I ) = Hu - Hd ;
    
    
    if I==1
        Nodtypelogical=( Nodtype( :, I )==1);
        Qt_u ( I ) = Qt_d ( I ) +sum( Qout0 ( : , I ).*  Nodtypelogical);        %本渠池分水口门计划分水量
    else
        
        Qt_u ( I ) = Qt_d ( I-1 );
    end
    if I == 1                                            % GateDischarge计算闸门开度
        Au = ( CANAL ( 3 , 1 , I ) + CANAL ( 4 , 1 , I ) * Hu ) * Hu ; %Au各渠段上游闸前过水面积
    else
        Au = ( CANAL ( 3 , Nout ( I - 1 ) + 1 , I - 1 ) + CANAL ( 4 , Nout ( I - 1 ) + 1 , I - 1 ) * Hu ) * Hu ;
    end
    GA ( I , : ) =  feval(@Q2Greal , I  , 0 , Hu , Hd , Qt_u ( I ) , 0 , Gb ( I ) , Au , Qt_u ( I ) , Gmax ( I ) , G , CDGATE ( I ) , DBGate , 1000 , repmat([1000 1000],[Ncnl 1]) , DBWatLev , SilenceMode ,  CalcAcc ) ;
    GA ( I , : ) =  round ( GA ( I , : ) * GateAccuracy  ) / GateAccuracy  ; %GA各渠段上游闸门开度
    QP1 ( I , : ) = Qt_u ( I ) ;                               % 各渠段初始上边界
    QPn ( I , : ) = Qt_d( I ) ;                               % 各渠段初始下边界
    YT ( I ) = Hd01 ( I ) + z0_d ( I ) ; %YT控制点恒定流目标水位
    [ V( I ) , Delay( I )] = feval(@Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 : N ( I ) ) ) ;% 各渠段水体体积
    Volu ( I , 1 ) = V ( I ) ;
    Delayl ( I , 1 ) = Delay ( I ) ;
    Qtal ( I , 1 ) = Qt_u ( I ) ;
    Qsum_up ( I ) = 0 ;                                   % 各渠段渠首流量累计水量
end
% 以上--------------------------------------------初始闸门开度------------------------------------以上

% 以下-----------------------------------------记录数据和参数定义---------------------------------以下

for I = 1 : Ncnl
    Tx ( I , 1 ) = 0 ;                                        % 时间数据横坐标
    Q_U ( I , 1 ) = QP1 ( I , 1 ) ;                           % 各渠段上游流量输出
    Qoutsum ( I , 1 ) = sum ( Qout0 ( : , I ) ) ;             % 各渠段区间取水口流量输出
    Z_U ( I , 1 ) = Z ( I , 1 ) ;                             % 各渠段上游水位输出
    Z_D ( I , 1 ) = Z ( I , N ( I ) ) ;                       % 各渠段下游水位输出
    Z_M ( I , 1 ) =0.5 * ( Z_U ( I , 1 )+ Z_D ( I , 1 ));
    GA1 ( I , 1 ) = GA( I , 1 ) ;                             % 各渠段上游闸门开度输出
    V1 ( I , 1 ) = V ( I ) ;                                  % 各渠段水体体积输出
    Vttl ( I , 1 ) = V ( I ) ;                                % 蓄量补偿前馈目标蓄量计算值输出
    YT1 ( I , 1 ) = YT ( I ) ;                                % 控制点恒定流目标水位输出
    YTorig_U ( I ) = Z ( I , 1 ) ;                            % 恒定流上游目标水位输出
    YTorig_D ( I ) = Z ( I , N ( I ) ) ;                      % 恒定流下游目标水位输出
    YTorig_UP ( I , 1 ) = YTorig_U ( I ) ;                    % 恒定流上游目标水位输出
    YTorig_DOWN ( I , 1 ) = YTorig_D ( I ) ;                  % 恒定流上游目标水位输出
    Velocity ( I , 1 ) = Velomax ( I ) ;
end
Zall ( 1 , : , : ) = Z ;                                      %不同时间点上各渠段的水位，Zall(5,3,11)表示第5个周期第三个渠段第11个节点上的水位
Qall ( 1 , : , : ) = Q ;                                      %不同时间点上各渠段的水位
% 以下%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%冰期程序所增部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以下
Twall ( 1 , : , : ) = Tw ;                                    %不同时间点上断面平均水温
hiall ( 1 , : , : ) = hi ;                                    %存储渠系各断面冰盖厚度的变化过程
Ciall ( 1 , : , : ) = Ci ;                                    %不同时间点上的流冰浓度
% 以上%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%冰期程序所增部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以上

Qta = Qt_u ;                                        % 目标流量（流量前馈用）
% 以上-----------------------------------------记录数据和参数定义---------------------------------以上
% 以上*********************************************初始条件***************************************以上

% Part V : Unsteady flow with canal control simulation
% 以下*********************************************计算部分***************************************以下
DB = Vmax * DT ;                                      %一个采样周期内闸门最大值增量
T = 0 ;                                               %时间（分钟）
T1 = 0 ;                                              % 变量初始化
if Tmax <= 72
    Tshow = 60 ;                             % pops time every 60 minutes
else if Tmax <= 720
        Tshow = 360 ;                        % pops time every 6 hours
    else
        Tshow = 1440 ;                       % Pops time every day
    end
end
while T < Tmaxm     %Tmaxm仿真总时间（分钟）
    T1 = T1 + 1 ;
    T = T + DT ;    %DT仿真时间步长，min
  if T==1440
      T=T;
  end
    for I = Ncnl : -1 : 1           %用于在系统流量产生变化时期关闭前馈环节
        %         if T < TE * NMinToHou
        %             PID ( I , 1 : 3 ) = 0 ;
        %         else
        PID ( I , 1 : 3 ) = PIDS ( I , 1 : 3 ) ;
        %         end
    end
    if Tshow == 60
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['目前仿真时间进行到第' , num2str( T ) , '分钟，仿真总时长为' , num2str( Tmax ) , '小时（ =' , num2str( Tmaxm ) , '分钟） 仿真进行中。' ] ) ;
        end
    elseif Tshow == 360
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['目前仿真时间进行到第' , num2str( T / NMinToHou ) , '小时，仿真总时长为' , num2str( Tmax ) , '小时，仿真进行中。' ] ) ;
        end
    else
        if floor ( T / Tshow ) == T / Tshow ;
            disp ( ['目前仿真时间进行到第' , num2str( T / ( NMinToHou * NHouToDay ) ) , '天，仿真总时长为' , num2str( Tmax / NHouToDay ) ,'天，仿真进行中===================================' ] ) ;
        elseif floor ( T / 180 ) == T / 180 ;
            disp ( ['目前仿真时间进行到第' , num2str( floor ( T / ( NMinToHou * NHouToDay ) ) ) ,'天' , num2str( T / NMinToHou - floor ( T / ( NMinToHou * NHouToDay ) ) * NHouToDay ) , '小时，仿真总时长为' , num2str( Tmax / NHouToDay ) , '天，仿真进行中。' ] ) ;
        end
    end
    if T < FailTime
        Qout0 = feval(@Qout , T, Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0 ) ;                           % 引用Qout
    else
        Qout0 = AccidentResponse ( Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, T , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout );
    end                          % 引用Qout
    Qsum_down = Qsum_down +feval(@ Qdown , T , Tdown_start , Qdown_start , Tdown_end , Qdown_end ) * DT * NMinToHou ;           % 下游总需水量   此语句疑似有问题，要需水流量的累积值做什么？
    % 以下---------------------------------如果流态发生变化，需重新计算目标流量---------------------------以下
    for I = Ncnl : -1 : 1
        if I == Ncnl
            Qt_d ( I ) = feval(@Qdown , T , Tdown_start , Qdown_start , Tdown_end , Qdown_end ) ;             % Qt_d各渠段下游流量 (Qt: 目标流量)
        else
            Qt_d ( I ) = Qt_d ( I + 1 ) + sum ( Qout0 ( : , I + 1 ) );
        end
    end
    for I = Ncnl : -1 : 1
        if I == 1
            Qt_u ( I ) = Qt_d ( I ) + sum ( Qout0 ( : , I ) );
        else
            Qt_u ( I ) = Qt_d ( I - 1 ); % Qt_u 各渠段计划上游流量过程 (Q1)
        end     %Qtt循环开始后的恒定流流量
    end
    for I = Ncnl : -1 : 1
        Qtt ( I , 1 : N ( I ) ) = zeros ( 1 , N ( I ) ) ;       % 渠道各断面流量数组定义
        Qtt ( I , 1 : N ( I ) ) =feval(@ Ini_Q , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , Qt_d ( I ) ) ;
        Vqtt ( I , 1 : ( Nout ( I ) + 1 ) ) = zeros ( 1 , Nout ( I ) + 1 ) ;       % 渠道各断面平均流速数组定义
        Vqtt ( I , 1 : ( Nout ( I ) + 1 ) )  = feval( @VSegmentavg , CANAL , Qtt , Z , TZ0 , I , N1 , Nout ( I ) ) ;               % 计算各子渠段平均流速
        Zouthqtt ( 1 : ( Nout ( I ) + 1 ) , I ) = zeros ( Nout ( I ) + 1 , 1 ) ;       % 渠道各断面平均流速数组定义
        Zouthqtt (  1 : ( Nout ( I ) + 1 ) , I ) =feval(@ HeadLoss , Nodtype2 , Kouth  ,Vqtt , Nout ( I ) , I );
        %for K = 1 : Nout ( I )
        %    HeadLoss ( K , I ) = ZOUTH ( Nodtype ( K , I ) , Kouth ( K , I ) , Qtt ( sum ( N1 ( 1 : K , I ) ) ) ) ;    %Qtt循环开始后的恒定流流量
        %end
        if Wt == 0                                    % 下游常水位运行
            Hdtt ( I ) = Hd01 ( I ) ; %Hdtt循环开始后恒定流水面线计算的下游水深
        else
            Hdtt = feval(@V2Hd , V1 ( I , 1 ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthqtt ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , [ 2 8 ] ) ;  % 等体积运行
        end
        Ztt =feval(@ Ini_Z , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I  ) + 1 , I ) , Nout ( I ) , Qout0 ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthqtt ( 1 : Nout ( I ) , I ) , Qt_d ( I ) , Hdtt ( I ) , I  ) ;
        Zttl (T1+1, I , 1 : N ( I ) ) = Ztt ; %Ztt循环开始后恒定流水面线计算
        YTorig_U ( I ) = Ztt ( 1 ) ;
        YTorig_D ( I ) = Ztt ( N ( I ) ) ;
        YT ( I ) = Hdtt ( I ) + TZ0 ( I , N ( I ) ) ;
        Vtt ( I ) =feval(@ Vol , Nout ( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Ztt ( 1 : N ( I ) ) , Q ( I , 1 :N ( I ) ) ) ;
        Vttl ( I , T1 + 1 ) = Vtt ( I ) ;
        if Delay ( I ) > DT
            Delay ( I ) = DT ;
        end
        if T < FailTime
%             sumQout = sum ( Qout ( ( T + Delay ( I ) ), Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0)  ) ;                           % sumQout 渠池间所有取水口流量之和
              sumQout ( 1 : Ncnl , T1 ) = sum ( Qout ( ( T ), Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0)  ) ;
        else
            sumQout ( 1 : Ncnl , T1 ) = sum ( AccidentResponse ( Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, T + Delay ( I ) , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout ) );
        end
        Dvttl ( I , T1 ) = Vttl ( I , T1 + 1 ) - Vttl ( I , T1 );
        if I == Ncnl
            Qta ( I ) = Qt_d ( I ) + sumQout ( I , T1 ); %+  Dvttl ( I , T1 ) / ( DT * NMinToHou * tk ( I ) ) ; %渠池首端流量等于尾端流量+区间取水+由蓄量变化引起的水量/时间步长
        else
            Qta ( I ) = Qta ( I + 1 ) + sumQout ( I , T1 ) ;%+  Dvttl ( I , T1 ) / ( DT * NMinToHou * tk ( I ) ) ;
        end
        Qtal ( I , T1 + 1 ) = Qta ( I ) ;
        if   T > NMinToHou * TS && T < NMinToHou * ( TE + 3 )   %TS流量变化开始时刻
            Qta ( I ) = Qta ( I ) ;
        else
            Qta ( I ) = Qt_u ( I ) ;
        end
    end
    % 以上---------------------------------如果流态发生变化，需重新计算目标流量---------------------------以上
    % 以下------------------------------------------控制器流量输出计算-----------------------------------以下
    for I = Ncnl : -1 : 1
        YF ( I ) = Wt * Z(  I, 1 ) + ( 1 - Wt ) * Z(  I, N ( I ) ) ;
        E ( I , 2 ) = YT ( I ) - YF ( I ) ;                        % E(I,2):I渠段此时刻的误差；E(I,1):I渠段上一时刻的误差；
        % 以下---------------------------------------------滤波-----------------------------------------以下
        Ef ( I , 2 ) = ( 1 - Kf ( I ) ) * E( I , 2 ) + Kf ( I ) * E ( I , 1 );    % 滤波后水位偏差
        ECf ( I , 2 )=( Ef ( I , 2 ) - Ef ( I , 1 ) ) / DT;            % 滤波|后|水位偏差率
        % 以上---------------------------------------------滤波-----------------------------------------以上
        EC ( I , 2 ) = ( ( E ( I , 2 ) - E ( I , 1 ) ) ) / DT ;              % 滤波|前|水位偏差率
        Jug1 ( I ) = Jug1 ( I ) + ( E ( I , 2 ) ^ 2 ) / ( YT ( I ) - 0.5 * ( z0_u ( I ) + z0_d ( I ) ) ) ;
        IAE1 ( I ) = IAE1 ( I ) + abs ( E ( I , 2 ) ) / ( YT ( I ) - 0.5 * ( z0_u ( I ) + z0_d ( I ) ) ) ;
        [ QP1( I , 2 ) , DQff(I,T/DT) , DQfb(I,T/DT) ]= feval(@FFwithFB , Ef ( I , 2 ) , ECf ( I , 2 ) , ECf ( I , 1 ) , PID ( I , 1 : 3 ) , Qta ( I ) , QP1 ( I , 1 ) , Qmax ( I ) ) ; % PID流量输出，带滤波
        
   
        % 以下---------------在流量变化期间变控制器输出时间步长，通过将非整时刻流量变化输出置零--------------以下
        if   T > NMinToHou * TS && T < NMinToHou * ( TE + 3 )
            if fix ( T_wave ( I ) / DT ) == T_wave ( I ) / DT
                QP1 ( I , 2 ) = QP1 ( I , 2 ) ;
            else
                QP1 ( I , 2 ) = QP1 ( I , 1 ) ;
            end
        end
        
       
        
        % 以上---------------在流量变化期间变控制器输出时间步长，通过将非整时刻流量变化输出置零--------------以上
    end
    % 以上------------------------------------------控制器流量输出计算------------------------------------以上
    
    % 以下---------------------------------------------闸门开度计算---------------------------------------以下
    Nodtype2 = Nodtype ;
    for I = 1 : Ncnl
        Nodtype2 ( Nout ( I ) + 1 , I ) = 5 ;                % 闸门节点，第Nout(I)+1个子渠段流速计算的需要
    end
    for I = Ncnl : -1 : 1
        if I==55
            I=I;
        end
        if  I == 1
            %Hu = Hu1 - Kouth ( Nout ( 1 ) + 1 , 1 ) * V0_d ( I , 1 ) ^ 2 /2/G ;
            Hu = Hu1;                     %15-9-22改
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
            disp ( [ ' T= ' , num2str( T ) , '，I= ' , num2str( I ) , '时，Hu-Hd=' , num2str( Hu - Hd )  , '闸后水位高于闸前水位，控制异常。'] ) ;
            pause ( 3 ) ; % 停顿3秒
        end
        GA ( I , 2 ) = feval(@Q2Greal , I , T , Hu , Hd , QP1 ( I , 2 ) , GA ( I , 1 ) , Gb ( I ) , Au , QP1 ( I , 2 ) , Gmax ( I ) , G , CDGATE ( I ) , DBGate , DB , E , DBWatLev , SilenceMode , CalcAcc ) ;
        GA ( I , 2 ) = round ( GA ( I , 2 ) * GateAccuracy ) / GateAccuracy ;
        if Hu < Hd
            QP1 ( I , 2 ) = QP1 ( I , 1 ) ;
        else
            QP1 ( I , 2 ) = feval(@GateDischarge , Hu , Hd , GA ( I , 2 ) , Gb ( I ) , Au , QP1 ( I , 2 ) , G , CDGATE ( I ) ) ; %cancel comment by GGH, 15.09.29  如果开启这句指令就会出现初始状态不稳，因为根据闸门精度反算的目标流量与初始流量出现偏差，使四点隐格式计算出的流量、水位出现误差
        end
        DQGA ( I , T1 ) = QP1 ( I , 2 ) - QP1 ( I , 1 ) ;
        if DQGA ( I , T1 ) < GateDischargeAcc
            DQGA (I , T1 ) = 0 ;
        end
        % 以上---------------------------------------------闸门开度计算---------------------------------------以上
        Qsum_up ( I ) = Qsum_up ( I ) + QP1 ( I , 2 ) * DT * NMinToHou ;          % 比较流量面积积分参考用
        % 以下---------------------------------------------非恒定流计算---------------------------------------以下
        if I == Ncnl
            QPn ( I ) = Qt_d ( I ) ;
        else
            QPn ( I ) = QP1 ( I + 1 , 2 ) ;
        end

        SV_DQP1 = ( QP1 ( I , 2 ) - Q ( I , 1 ) ) / SV_N ;
        SV_DQPn = ( QPn ( I ) - Q( I , N ( I ) ) ) / SV_N ;
        for SV_I = 1 : SV_N  %在系统的一个动作时间步长内算多次非恒定流过渡过程
            if T < FailTime
                SV_Qout =feval(@ Qout , (T - DT + SV_I * SV_DT) , Ksn, Csn, Qstart, Qend, Tstart, Tend, Ncout, Qout0) ;                           % 引用Qout
            else
                SV_Qout =feval(@ AccidentResponse , Ksn, Csn, Qstart, Qend, Tstart, Tend, Qoutdsgn, Qoutype, CanalAccidentUp, CanalAccidentDn, ( T - DT + SV_I * SV_DT) , FailTime,  Z( : , : ), ZJD ( : ), N, Qout0, Ncout );
            end
            %for K = 1 : Nout ( I )
            %    Zouth ( K , I ) = ZOUTH ( Nodtype ( K , I ) , Kouth ( K , I ) , Q ( I , sum ( N1 ( 1 : K , I ) ) ) ) ;
            %end
            if simtype == 1     %计算结冰
                if hiall ( T1 , I , N ( I ) ) > 0   %在冰盖厚度大于零时做冰水力学计算
                    [Q( I , 1 : N ( I ) ) , Z( I  , 1 :N ( I ) )] = feval(@s_Preice , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) , hiall ( T1 , I , 1 : N ( I ) ) , T , TTT ( I , : ) , IceDen , WatDen , IniRough , FinlRough ) ;
                else                                %在冰盖厚度等于零时还是做常规非恒定流计算
                    [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) )] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
                end
            else                %不计算结冰
                [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) ] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) ,Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q( I , N(I) ) + SV_DQPn , 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
          
              %  [Q( I , 1 : N ( I ) ) , Z( I , 1 : N ( I ) ) ] =feval(@ s_Pre , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , SV_DT * NMinToHou , N1 ( 1 : Nout ( I ) + 1 , I ) , Q ( I , 1 : N ( I ) ) ,Z ( I , 1 : N ( I ) ) , Q ( I , 1 ) + SV_DQP1 , Q ( I , N ( I ) ) + SV_DQPn, 1 , Nout ( I ) , SV_Qout ( 1 : Nout ( I ) , I ) , TZ0 ( I , 1 : N ( I ) ) , Zouthunst ( 1 : Nout ( I ) , I ) ) ;
            end
        end
        Hall = Z - TZ0 ;
        % 以上---------------------------------------------非恒定流计算---------------------------------------以上
        % 以下%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%冰的计算部分%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%以下
        if simtype == 1
            % 以下-------------------------------------------太阳辐射-----------------------------------------以下
            T_radi = floor ( T / ( NMinToHou * NHouToDay ) ) + 1 ;
            T_sunrise = ( T_radi - 1 ) * NHouToDay * NMinToHou * NMinToHou + 6 * NMinToHou * NMinToHou ;    % 日出时刻6点
            T_sunset = ( T_radi - 1 ) * NHouToDay * NMinToHou * NMinToHou + 18 * NMinToHou * NMinToHou ;    % 日落时刻18点
            if T * NMinToHou > T_sunrise && T * NMinToHou <= T_sunset                                       % 日间12小时半正弦太阳辐射分布
                Radi = radi ( T_radi ) * ( 0.5 * pi / ( 12 * NMinToHou * NMinToHou ) ) * sin ( ( T * NMinToHou - T_sunrise ) * pi / ( T_sunset - T_sunrise ) ) ;
            else
                Radi = 0 ;
            end
            % 以上-------------------------------------------太阳辐射------------------------------------------以上
            for j = 1 : N ( I ) %此循环过于臃肿，应该按照渠段-子渠段-断面的形式多一层循环会更加简洁
                for II = 1 : Nout ( I ) + 1
                    if j - sum ( N1 ( 1 : II , I ) ) <= 0
                        a = II ;        % a为j当前所在渠池的子段号，用于读取渠道断面参数
                        break ;
                    end
                end
                D ( I , j ) = Zall ( T1 , I , j ) - TZ0 ( I , j ) - IceDen * hi ( I , j ) / WatDen ;        % 等效水深
                B ( I , j ) = CANAL ( 3 , a , I ) + 2 * CANAL ( 4 , a , I ) * D ( I , j ) ;                 % 水面宽度
                A ( I , j ) = ( CANAL ( 3 , a , I ) + CANAL ( 4 , a , I ) * D ( I , j ) ) * D ( I , j ) ;   % 过水断面面积
                u ( I , j ) = Qall ( T1 , I , j ) / A ( I , j ) ;                                           % 断面平均流速
                DD ( I , j ) = A ( I , j ) / B ( I , j ) ;                                                  % 水深？
                if u ( I , j ) > u1max ( I , j )
                    u1max ( I , j ) = u ( I , j ) ;
                end
                D1 ( I , j ) = Z ( I , j ) - TZ0 ( I , j ) - IceDen * hi ( I , j ) / WatDen ;               %D1与D有什么区别。Zall是上一时刻的水位？
                B1 ( I , j ) = CANAL ( 3 , a , I ) + 2 * CANAL ( 4 , a , I ) * D1 ( I , j ) ;
                A1 ( I , j ) = ( CANAL ( 3 , a , I ) + CANAL ( 4 , a , I ) * D1 ( I , j ) ) * D1 ( I , j ) ;
                u1 ( I , j ) = Q ( I , j ) / A1 ( I , j ) ;
                Fr1 ( I , j ) = u1 ( I , j ) / sqrt ( G * D1 ( I , j ) ) ;   % 佛汝德数
                DD1 ( I , j ) = A1 ( I , j ) / B1 ( I , j ) ;
                Ta ( T1 , I ) = Temp ( I,T1 * DT * NMinToHou ) ;              % 读取气温值
                C( T1 , I ) = CloudAmount ( I,T1 * DT * NMinToHou ) ;         % 读取云量
                Va( T1 , I ) = WindVelo ( I,T1 * DT * NMinToHou ) ;           % 读取风速
                CAA ( I , : ) = CAAconst ;                                  % 设定面流冰密度封冻标准
                d0 ( I , : ) = IniIceThick ;                                % 设定流冰层厚度
                Ds1 = CANAL ( 5 , a , I ) / ( N1 ( a , I ) - 1 ) ;          % 计算步长
                hiw = 1622 * u ( I , j ) ^ 0.8 * DD ( I , j ) ^ 0.2 ;       % 冰与水的热交换系数
                %计算相邻两个计算断面间的热量变化
                fai ( I , j ) = HeatLossCalc ( Ta ( T1 , I ) , Radi , C( T1 , I ) , Va( T1 , I ) , hiall ( T1 , I , j ) , Twall ( T1 , I , j ) , Ca ( I , j ) , Ciall ( T1 , I , j ) ) ;  % 计算水体净热量损失
                if Twall ( T1 , I , j ) < 0
                    %水温低于零度时，水结成冰放出热量（潜热）使水温变为零度
                    Ci ( I , j ) = - ( WatDen * HeatCapa * Twall ( T1 , I , j ) ) / ( IceDen * LatentHeat ) + Ciall ( T1 , I , j ) ;
                    Tw ( I , j ) = 0 ;
                    if hiall ( T1 , I , j ) == 0        % 计算面流冰密度以判断是否形成冰盖
                        Ca ( I , j ) = Ci ( I , j ) * DD ( I , j ) / d0 ( I ) ;   %将流冰浓度换算成面流冰密度    如果冰盖厚度不等于零呢？不需要计算此时面流冰密度么？
                    end
                else
                    if j == 1
                        if I == 1
                            Tw ( I , j ) = UpTemp ;     % 第一个节点的水温等于上游来水温度
                        else                            %每个渠池第一个节点的温度计算（除去第一个渠池）
                            if hiall ( T1 , I , j ) > 0 % 有冰盖时，水温不受大气热量交换影响，热量主要用于生冰
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - hiw * Twall ( T1 , I , j ) / ( WatDen * HeatCapa * DD ( I , j ) ) ;
                            elseif Ciall ( T1 , I , j ) > 0
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - DT * NMinToHou * ( fai ( I , j ) + hiw * Twall ( T1 , I , j ) ) / ( WatDen * HeatCapa * DD ( I , j ) ) ;
                            else
                                Tw ( I , j ) = Twall ( T1 , I - 1 , N ( I - 1 ) ) - DT * NMinToHou * fai ( I , j ) / ( WatDen * HeatCapa * DD ( I , j ) ) ; %温度变化=上个节点温度-区间热量损失
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
                                if hiall ( T1 , I , j + 1 ) > 0    % 冰花在冰盖前缘处堆积
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
                            if Ca ( I , j ) >= CAA ( I )               % 达到封冻标准则封冻
                                if Ca ( I , j ) >= 1
                                    hi ( I , j ) = Ca( I , j ) * d0 ( I ) ;   % Vi(I,j)/(B1(I,j)*Ds1);
                                    Vii1 ( I , j ) = 0 ;               % 空隙率
                                else
                                    hi ( I , j ) = d0 ( I ) ;
                                    Vii1 ( I , j ) = B1 ( I , j ) * d0 ( I ) * ( 1 - Ca ( I , j ) ) ;  % 空隙的体积
                                end
                                Ci ( I , j ) = 0 ;
                                TTT ( I , j ) = T1 ;
                            else                               % 没达到封冻标准则进行流冰迭代
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
                                        if j == sum ( N1 ( 1 : a - 1 , I ) ) + 1    % 节点后第一个断面
                                            if Nodtype ( a - 1 , I ) == 3    % 节点为倒虹吸
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
                                        Ci ( I , j ) = Ci ( I , j ) + Ciall ( T1 , I , j ) ; % 冰花在渠段末拦冰锁处堆积
                                    elseif Nodtype ( a , I ) == 3             % 冰花在倒虹吸前拦冰锁处堆积
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
            Velo ( J , I ) = Vcmax ( CANAL , Q , Z , TZ0 , I , J , N1 ) ;   % 计算各子渠段最大流速
        end
        Velomax ( I ) = F_velomax( VeloENG (1 : Nout ( I ) + 1 , I ) , Nout ( I ) , Limitlength , Nodtype2 ( 1 : Nout ( I ) + 1 , I ) , Velo ( 1 : Nout ( I ) + 1 , I ) , CANAL ( 5 ,1 : Nout ( I ) + 1 , I ) );
        [V( I ) , Delay( I )] = Vol ( Nout( I ) , CANAL ( : , 1 : Nout ( I ) + 1 , I ) , N1 ( 1 : Nout ( I ) + 1 , I ) , TZ0 ( I , 1 : N ( I ) ) , Z ( I , 1 : N ( I ) ) , Q ( I , 1 : N ( I ) ) ) ;% 各渠段水体体积
        % 以下-------------------------------------记录上一时刻的各变量-----------------------------------------以下
        IAQ1 ( I ) = IAQ1 ( I ) + abs ( QP1 ( I , 2 ) - QP1 ( I , 1 ) ) ;
        IAW1 ( I ) = IAW1 ( I ) + abs ( GA ( I , 2 ) - GA ( I , 1 ) ) ;
        E ( I , 1 ) = E ( I , 2 ) ;
        EC ( I , 1 ) = EC ( I , 2 ) ;
        Ef ( I , 1 ) = Ef ( I , 2 ) ;
        ECf ( I , 1 ) = ECf ( I , 2 ) ;
        GA ( I , 1 ) = GA ( I , 2 ) ;
        QP1 ( I , 1 ) = QP1 ( I , 2 ) ;
        % 以上-------------------------------------记录上一时刻的各变量-----------------------------------------以上
        T2 = T1 / 1 + 1 ;
        if T2 == fix ( T2 )    % 整数时刻记录数据
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
            Z_U ( I , T2 ) = Z ( I , 1 ) ;                                      %渠道最上游断面水位
            Z_D ( I , T2 ) = Z ( I , N ( I ) ) ;                                %渠道最下游断面水位
            Z_M ( I , T2 ) = 0.5 * Z ( I , 1 ) + 0.5 * Z ( I , N ( I ) ) ;
            Volu ( I , T2 ) = V ( I ) ;
            Delayl ( I , T2 ) = Delay ( I ) ;
            Ddz_U ( I , T2 ) = ZJU ( I ) - Z ( I , 1 ) ;                         % 上游校核水位和实际水位差
            Ddz_D ( I , T2 ) = ZJD ( I ) - Z ( I , N ( I ) ) ;                   % 下游校核水位和实际水位差
            Ddz_M ( I , T2 ) = Ddz_U ( I , T2 ) * 0.5 + Ddz_D ( I , T2 ) * 0.5;  % 中游校核水位和实际水位差
            EH_U ( I , T2 ) = Z ( I , 1 ) - YTorig_U ( I ) ;                     % 各渠段上游水位误差
            EH_D ( I , T2 ) = Z ( I, N ( I ) ) - YTorig_D ( I ) ;                % 各渠段下游水位误差
            EH ( I , T2 ) = 0.5 * EH_U ( I , T2 ) + EH_D ( I , T2 ) * 0.5 ;      % 各渠段中游水位误差
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
                        disp ( ['渠池',num2str( I ),'上游1h水位变幅超限！'] )
                    elseif abs ( DV_D1 ( I , T2 ) ) > MDS1
                        disp ( ['渠池',num2str( I ),'下游1h水位变幅超限！'] )
                    elseif abs ( DV_M1 ( I , T2 ) ) > MDS1
                        disp ( ['渠池',num2str( I ),'中游1h水位变幅超限！'] )
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
                        disp ( ['渠池' , num2str( I ) , '上游24h水位变幅超限！'] )
                    elseif abs ( DV_D2 ( I , T2 ) ) > MDS2
                        disp ( ['渠池',num2str( I ),'下游24h水位变幅超限！'] )
                    elseif abs ( DV_M2 ( I , T2 ) ) > MDS2
                        disp ( ['渠池' , num2str( I ) , '中游24h水位变幅超限！'] )
                    end
                end
            else
                DV_U2 ( I , T2 ) = Z_U ( I , T2 ) - Z_U ( I , 1 ) ;
                DV_D2 ( I , T2 ) = Z_D ( I , T2 ) - Z_D ( I , 1 ) ;
                DV_M2 ( I , T2 ) = Z_M ( I , T2 ) - Z_M ( I , 1 ) ;
            end
        end
        for K = Ncnl : -1 : 1                         % 计算渠段上下游累计过水量之差
            if I == Ncnl                          % 2010-10-20将循环变量由I改为K
                DQsum_up ( K ) = Qsum_up ( K ) - Qsum_down ;
            else
                DQsum_up ( K ) = Qsum_up ( K ) - Qsum_up ( K + 1 ) ;
            end
        end
    end
    Qout1 = Qout0 ;                                % 记录上一次的Qout0
    Q01_d = Qt_d ( Ncnl ) ;                           % 记录上一时刻下游需水流量变化
    Zall ( T2 , : , : ) = Z ;                             % 记录不同时间点上各渠段的水位参数，Zall(5,3,11)表示第5个周期第三个渠段第11个节点上的水位
    Qall ( T2 , : , : ) = Q ;
    %冰期
    Twall ( T2 , : , : ) = Tw ;
    hiall ( T2 , : , : ) = hi ;
    Ciall ( T2 , : , : ) = Ci ;
    %冰期
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
            if HiceD1 ( I , K ) > 0 % 截至第1天时刻（01天00:00:00）结冰范围
                if K > sum ( N1 ( 1 : Jice ,I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_1d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_1d ( I ) = 0;
            end
            if HiceD3 ( I , K ) > 0 % 截至第1天时刻（03天00:00:00）结冰范围
                if K > sum ( N1 ( 1 : Jice , I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_3d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_3d ( I ) = 0;
            end
            if HiceD5 ( I , K ) > 0 % 截至第1天时刻（05天00:00:00）结冰范围
                if K > sum ( N1 ( 1 : Jice  , I ) ) && K < sum ( N1 ( 1 : Jice + 1 , I ) )
                    len_5d ( I ) = max ( sum ( CANAL ( 5 , Jice(I) : Nout ( I ) + 1 , I ) ) );
                end
            else
                len_5d ( I ) = 0;
            end
        end
    end
end
Tw_u = Twall ( : , : , 1 ) ;    % 渠池闸后水温
hi_u = hiall ( : , : , 1 ) ;    % 渠池闸后水温
[ XTw_u, YTw_u ] = size ( Tw_u );
Tw_d = zeros( XTw_u, YTw_u );
for I = Ncnl : -1 : 1
    Tw_d = Twall ( : , : , N(I) ) ; % 渠池下游水温
    hi_d = hiall ( : , : , N(I) ) ; % 渠池下游冰厚
    IAQ ( I ) = IAQ1 ( I ) - abs ( Q_U ( I , 1 ) - Q_U ( I , T2 ) ) ;
    IAW ( I ) = IAW1 ( I ) - abs ( GA1 ( I , 1 ) - GA1 ( I , T2 ) ) ;
end
IAE = IAE1 * DT / NMinToHou / Tmaxm ;
Jug = Jug1 * DT / Tmaxm ;
Jug_sum = sum ( Jug1 ) * DT / Tmaxm ;             % 寻优参数
disp ( ['当前时间为______' , datestr( now, 31 ) ] ) ;
disp ('=================仿真结束===================') ;
save SIM_RES
Tsim_S = toc;
if SilenceMode == 0
    disp ('============仿真结果保存在目标目录下的SIM_RES.mat文件中============') ;
    disp ( ['仿真耗时' , num2str( Tsim_S / NMinToHou ) , '分钟'] );
end
Output;
if SilenceMode == 0
    disp ('===闸门开度、过闸流量、控制点水位、最大流速、冰盖厚度、渠池水温、冰盖范围分别以txt格式保存在目标目录下===') ;
    disp ('===闸门开度、过闸流量、控制点水位、最大流速、冰盖厚度、渠池水温、冰盖范围分别以xls格式保存在目标目录下===') ;
    pause ( 5 ) ; % 停顿5秒
end
% 以下****************************************仿真结果绘图****************************************以下
if SilenceMode == 0
    disp ( '是否需要对仿真数据进行选择性绘图 1 = 是 任意键 = 否' )
    Choice = input (' 请作出您的选择:','s');
    switch Choice
        case '1'
            disp ('======将仿真结果SIM_RES.mat文件中的数据可选择性进行绘图========');
            PPNi
        otherwise
            disp ('=============不需绘图，将数据输出，仿真结束===============');
    end
end
%以上****************************************仿真结果绘图****************************************以上
