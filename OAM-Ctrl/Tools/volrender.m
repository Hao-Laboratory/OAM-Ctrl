function h = volrender(xs,ys,zs,v,cmap,amap)
h = vol3d('xdata',xs,'ydata',ys,'zdata',zs,'cdata',v);

box on;
colormap(cmap);
alphamap(amap);
view(45,30);
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');