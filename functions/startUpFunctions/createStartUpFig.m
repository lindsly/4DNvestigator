tempH = load('HiC1M_FIB_raw.mat', 'H', 'chr');
tempR = load('RNAseq1M.mat', 'RNAseq_av');

chrSelect = 14;
tempChr = tempH.chr;
a = log(tempH.H(tempChr(chrSelect):tempChr(chrSelect+1)-1,...
    tempChr(chrSelect):tempChr(chrSelect+1)-1,1:8));
c = log(tempR.RNAseq_av(tempChr(chrSelect):tempChr(chrSelect+1)-1,:)+1)';

% trim Chr
a(~isfinite(a)) = 0;
a(isnan(a)) = 0;
remLocs = any(squeeze(sum(a,1) == 0),2);

a(remLocs,:,:) = [];
a(:,remLocs,:) = [];
c(:,remLocs) = [];

% begin plot
outline = 0;
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

[x,y,z] = meshgrid(1:1:size(a,2),1:1:size(a,1),-size(a,3)+1:1:0);

h = slice(app.UIAxes,x,y,z,a,xslice,yslice,zslice);
set(h,'EdgeColor','none')
xlabel(app.UIAxes,'Time');ylabel(app.UIAxes,'y');zlabel(app.UIAxes,'z');
ylim(app.UIAxes,[-1 size(b,1)+1])
daspect(app.UIAxes,[scale 1 1])
view(app.UIAxes,-61,22)
axis(app.UIAxes,'off')
hold(app.UIAxes,'on')
app.UIAxes.XLabel.Visible = 'on';
title(app.UIAxes,sprintf('Proliferating Fibroblast 4DN data, Chr%i',chrSelect))

if outline == 1
    for i = 1:size(b,3)
        plot3(app.UIAxes,[i,i],[1,1],[0,-size(b,1)],'k-')
        plot3(app.UIAxes,[i,i],[1,size(b,1)+1],[-size(b,1),-size(b,1)],'k-')
        plot3(app.UIAxes,[i,i],[size(b,1)+1,size(b,1)+1],[-size(b,1),0],'k-')
        plot3(app.UIAxes,[i,i],[size(b,1)+1,1],[0,0],'k-')
    end
end
colormap(app.UIAxes,'hot')

bar3_sr(app.UIAxes,flipud(c'))