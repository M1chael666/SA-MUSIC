%%动平台MUSIC算法

function [X,Y,Z] = move_music_1(Pos_signal,Pos_receive,X_search,Y_search,lambda,snr,fig_mark)
K = size(Pos_signal,1);                        %信号源个数                                           
M = size(Pos_receive,1);                       %阵元个数                                        
P = M;                                         %节点个数
Z_search = Pos_signal(1,3);                    %信号源Z坐标


%% search arrival direction                                     
V = 7000;                                      %卫星速度
N_sample = 1000;                                %采样点数 
t = 1/(100e6);      %！这个值会影响定位性能！      %快时域时间间隔，采样频率100MHz 
Num = 1;                                      %集中采样次数
T = 0.5;                                       %慢时域时间间隔     
None_sample = T/t;
[DisX,DisY,DisZ] = DisT(V,T);                  %慢时域位移
signal = [];                                   %LFM信号                                                             
A_MxK_move = [];                               %导向矢量
yin = zeros(M,Num*N_sample);
Pos_receive_1 = Pos_receive;
for kk= 1:K
    A_MxK_move = [];
    Pos_receive = Pos_receive_1;
    signal(kk,:) = randn(1,Num*N_sample);
    %signal(kk,:) = sig_generation(1,Num*N_sample);
 for ii = 1:Num
    %signal =[signal,sig_generation((ii-1)*None_sample+1,N_sample)];       %采样时间未超出脉宽200ms 
    %signal(kk,:) = [signal(kk,:),sig_generation(1,N_sample)];             %采样时间超出脉宽
    [receive_x,receive_y,receive_z] = P_move(Pos_receive,N_sample,t,V);
    for jj = 1:N_sample
        A_MxK_move(:,(ii-1)*N_sample+jj) = asteer_far(Pos_signal(kk,:),[receive_x(:,jj),receive_y(:,jj),receive_z(:,jj)],lambda);
    end
    Pos_receive = Pos_receive + [DisX,DisY,DisZ];
 end
 yin = yin + A_MxK_move*diag(signal(kk,:));
end
%%调试
%A_MxK_move_com = Compensate_1(Pos_signal,A_MxK_move,N_sample,Num,lambda,t,T,V); 
c = 3e8;
f0 = c/lambda;
       

%yin = A_MxK_move*diag(signal);
yin = awgn(yin,snr,'measured');
%yin_com  = Compensate_1([X_search(ii),Y_search(jj),Z_search],yin,N_sample,Num,lambda,t,T,V);   %只补偿Num次
%[Es,En] = Eigen_decom(yin_com,K,M);
 
compensate_mark = 1;        %对比实验 1为只补偿Num次 0为不补偿 2为全补偿
%%
%%未完成
if compensate_mark == 2
P_MUSIC = zeros(length(X_search),length(Y_search));
 for ii = 1:length(X_search)
     for jj = 1:length(Y_search)
         yin_com  = Compensate([X_search(ii),Y_search(jj),Z_search],yin,N_sample,lambda,t,V);        %全补偿，共Num*N_sample次
         [Es,En] = Eigen_decom(yin_com,K,M);
         asteer = asteer_far([X_search(ii),Y_search(jj),Z_search],[receive_x(:,N_sample*Num*10),receive_y(:,N_sample*Num*10),receive_z(:,N_sample*Num*10)],lambda);
         P_MUSIC(ii,jj) = M-(asteer'*En*En'*asteer);          
     end
 end
end
%%
if compensate_mark == 1
P_MUSIC = zeros(length(X_search),length(Y_search));
 for ii = 1:length(X_search)
     for jj = 1:length(Y_search)
         yin_com  = Compensate_1([X_search(ii),Y_search(jj),Z_search],yin,N_sample,Num,lambda,t,T,V);   %只补偿Num次
         [Es,En] = Eigen_decom(yin_com,K,M);
         asteer = asteer_far([X_search(ii),Y_search(jj),Z_search],[receive_x(:,N_sample),receive_y(:,N_sample),receive_z(:,N_sample)],lambda);
         P_MUSIC(ii,jj) = (asteer'*En*En'*asteer)^-1;          
     end
 end
end
%%
if compensate_mark == 0
P_MUSIC = zeros(length(X_search),length(Y_search));
 for ii = 1:length(X_search)
     for jj = 1:length(Y_search)
         [Es,En] = Eigen_decom(yin,K,M);                                                             %不补偿 
         asteer = asteer_far([X_search(ii),Y_search(jj),Z_search],[receive_x(:,N_sample*Num*10),receive_y(:,N_sample*Num*10),receive_z(:,N_sample*Num*10)],lambda);
         P_MUSIC(ii,jj) = M-(asteer'*En*En'*asteer);          
     end
 end
end
%%
Pmax = max(max(abs(P_MUSIC)));
 %归一化
if(fig_mark)
figure
mesh(Y_search,X_search,abs(P_MUSIC));
xlabel('y/(m)','FontSize',13);
ylabel('x/(m)','FontSize',13);
zlabel('谱功率/dB','FontSize',13);
else
end


%%%%%%%%%%% 搜索谱峰
[max1,locs1] = max(abs(P_MUSIC));
[~,locs2] = max(max1);
Y = Y_search(locs2);
X = X_search(locs1(locs2));
Z = Z_search;
end