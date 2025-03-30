function [DisX,DisY,DisZ] = DisT(V,T)
    dread = pi/180;  %角度转弧度制
    V_alpha=[145;138;133;142;-309];
    %V_alpha =[196;145;208;133;138;-309;142;214;297;122;-148];   %速度的方向角
    DisX = cos(V_alpha.*dread)*V*T;
    DisY = sin(V_alpha.*dread)*V*T;
    DisZ = zeros(length(DisX),1);
end