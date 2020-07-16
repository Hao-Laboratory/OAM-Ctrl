function pupilshow(pupil_arg)
%PUPILSHOW �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

M = length(pupil_arg);
x = linspace(-1,1,M);
y = x;
[xx,yy] = meshgrid(x,y);
r = xx.^2 + yy.^2;

pupil = ones(M);
pupil(r>1) = nan;
h = imagesc(pupil_arg);axis image xy off;
set(h,'AlphaData',~isnan(pupil));
end