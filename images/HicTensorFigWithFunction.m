function [ h ] = HicTensorFigWithFunction( a,c,outline )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    outline = 0;
end

scale = .01;%1/size(a,3);

b = a;

for i = 1:size(a,3)
    aa(:,:,i) = [a(:,1,i),a(:,:,i)];
    aaa(:,:,i) = [aa(1,:,i);aa(:,:,i)];
    aaa(:,:,i) = rot90(aaa(:,:,i)',2);
end
a = aaa;

a = permute(a,[1,3,2]);
zslice = [];
yslice = [];
xslice = [1:size(a,2)];

h = figure;
[x,y,z] = meshgrid(1:1:size(a,2),1:1:size(a,1),-size(a,3)+1:1:0);

slice(x,y,z,a,xslice,yslice,zslice)
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
xlabel('x');ylabel('y');zlabel('z');
ylim([-1 size(b,1)+1])
% zlim([-1 size(a,1)+1])
daspect([scale 1 1])
view(-71,12)
axis off
hold on

if outline == 1
    for i = 1:size(b,3)
        plot3([i,i],[1,1],[0,-size(b,1)],'k-')
    plot3([i,i],[1,size(b,1)+1],[-size(b,1),-size(b,1)],'k-')
    plot3([i,i],[size(b,1)+1,size(b,1)+1],[-size(b,1),0],'k-')
    plot3([i,i],[size(b,1)+1,1],[0,0],'k-')
    end
end
colormap hot
colorbar

bar3_sr(flipud(c'))
end

