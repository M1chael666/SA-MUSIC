%%自相关函数特征值分解
%输入： X1  接收信号
%      K    信号源个数 
%      M    阵元个数
%输出：En    噪声子空间
function [Es,En] = Eigen_decom(X1,K,M)
Rxx = X1*X1'/M;         %计算协方差矩阵
[Ev,D] = eig(Rxx);      %特征值分解
EVA = diag(D)';         %将特征值矩阵对角线拍成一行
[EVA,I] = sort(EVA);    %将特征值从小到大排序
Ev = fliplr(Ev(:,I));   %对应特征矢量排序
Es = Ev(:,1:K);
En = Ev(:,K+1:M);

end