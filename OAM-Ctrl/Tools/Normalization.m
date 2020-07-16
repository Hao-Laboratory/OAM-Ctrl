function [OutputMatrix] = Normalization(InputMatrix)
%NORMALIZATION normalizes all data in of input matrix
%
% -------------------------------------------------------------------------
% coded by Liu,Xin
% Email: liuxin2018@zju.edu.cn
% Jan.4, 2018
% -------------------------------------------------------------------------

OutputMatrix = InputMatrix./max(InputMatrix(:));

end