% SDP_CONE plots the VSDP logo.

hold on;
axis equal;
[X1, X2] = meshgrid (0:.08:2);
X3 = sqrt(X1.*X2);

axis ([0,2,-2,2,0,2]);
mesh (X2,  X3, X1,'edgecolor','none');
mesh (X2, -X3, X1,'edgecolor','none');
shading interp;

view([-40,30]);
print -dpng sdp_cone.png
