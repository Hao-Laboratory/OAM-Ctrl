function struct=AL(theta)
% ��ɫ��͸����AL���ṹ����,���ô���Ϊ
% lensStructureOption��1
% first version compiled by Hao,Xiang 2010-6-2
% last modified and/or upgraded by Hao,Xiang 2010-6-2
% 
% �������������
% ------------------------------------------------------------------------
% theta:   ������XZƽ��ͶӰ��Z�ᣨ����)����ļнǣ�0<=theta<90
% phi:     ������XZƽ��ͶӰ��X��������ļнǣ�0<=phi<=2��
% �����޸���ʷ
% ------------------------------------------------------------------------
% 2010-6-2              v0.1����

struct=sqrt(cos(theta));