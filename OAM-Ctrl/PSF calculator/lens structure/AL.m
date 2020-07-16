function struct=AL(theta)
% 消色差透镜（AL）结构参数,调用代码为
% lensStructureOption＝1
% first version compiled by Hao,Xiang 2010-6-2
% last modified and/or upgraded by Hao,Xiang 2010-6-2
% 
% 输入变量声明：
% ------------------------------------------------------------------------
% theta:   光线在XZ平面投影与Z轴（光轴)方向的夹角，0<=theta<90
% phi:     光线在XZ平面投影与X轴正方向的夹角，0<=phi<=2π
% 编译修改历史
% ------------------------------------------------------------------------
% 2010-6-2              v0.1建立

struct=sqrt(cos(theta));