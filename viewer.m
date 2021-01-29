clear variables
close all;
clc;

set(groot,'defaulttextinterpreter','latex');  
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');

global maxL

%% Main Code

u = load('output/displacement.txt');
nXYZ = load('output/coordinates.txt');

x = unique(nXYZ(:,1));
y = unique(nXYZ(:,2));
z = unique(nXYZ(:,3));
maxL = max([x;y;z]);

xLen = size(x,1);
yLen = size(y,1);
zLen = size(z,1);
nLen = size(u,1)/3;

[Y,X,Z] = meshgrid(y,x,z);
Ux = zeros(xLen,yLen,zLen);
Uy = zeros(xLen,yLen,zLen);
Uz = zeros(xLen,yLen,zLen);

for i = 1:xLen
    for j = 1:yLen
        for k = 1:zLen
               
            idx = (i-1)*yLen*zLen+(j-1)*zLen+k;
            Ux(i,j,k) = u(idx+0*nLen);
            Uy(i,j,k) = u(idx+1*nLen);
            Uz(i,j,k) = u(idx+2*nLen);
        end
    end
end

xs = x(end);
ys = y(end);
zs = z(1);

%% Figure x

fig = figure(1);
slice(Y,X,Z,Ux,ys,xs,zs);

param(fig,'$u_x(x,y,z)$ [-]')
save(fig,'output/ux')

%% Figure y

fig = figure(2);
slice(Y,X,Z,Uy,ys,xs,zs);

param(fig,'$u_y(x,y,z)$ [-]')
save(fig,'output/uy')

%% Figure z

fig = figure(3);
slice(Y,X,Z,Uz,ys,xs,zs);

param(fig,'$u_z(x,y,z)$ [-]')
save(fig,'output/uz')

%% Function

function param(fig,name)
    
    set(fig,'visible','off'); 
    set(gca,'FontSize',12)
    
    xlabel('$y$ [-]')
    ylabel('$x$ [-]')
    zlabel('$z$ [-]')
    
    global maxL
    xlim([0,maxL])
    ylim([0,maxL])
    zlim([0,maxL])
    
    cb = colorbar;
    cb.Label.String = name;
    cb.Label.Interpreter = 'latex';
end

function save(fig,name)

    set(gcf, 'color', 'none');    
    set(gca, 'color', 'none');

    set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
    print(fig,name,'-painters','-dpdf','-r500')
end
