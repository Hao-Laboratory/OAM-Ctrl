function [Sx,Sy,Sz] = PoyntingVector(Ex,Ey,Ez,Hx,Hy,Hz)
%POYNTINGVECTOR calculates the poynting vector of 
% electromagnetic field
%   
% coded by Xin Liu
% email: liuxin2018@zju.edu.cn
% Apr.14, 2020

Hx_conj = conj(Hx);
Hy_conj = conj(Hy);
Hz_conj = conj(Hz);
Sx = real(Ey.*Hz_conj-Ez.*Hy_conj);
Sy = real(Ez.*Hx_conj-Ex.*Hz_conj);
Sz = real(Ex.*Hy_conj-Ey.*Hx_conj);
end