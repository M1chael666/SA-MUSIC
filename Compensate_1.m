%%输入：
%导向矢量：

%%输出：
%补偿后导向矢量：

function out = Compensate_1(Pos_signal,A_MxK_move,N_sample,Num,lambda,t,T,V)
[x_dis,y_dis,z_dis] = DisT(V,T);
R_dis = [x_dis,y_dis,z_dis];                        %集中采样间的阵元位移 
u = sig_u(Pos_signal);                              %单位向量   列向量
com = exp(1i*2*pi/lambda*R_dis*u);                  %相位补偿因子   列向量

for ii = 1:Num
    A_MxK_com(:,(ii-1)*N_sample+1:ii*N_sample) = A_MxK_move(:,(ii-1)*N_sample+1:ii*N_sample).*(com.^(Num-ii));    %迭代补偿
end

out = A_MxK_com;
end