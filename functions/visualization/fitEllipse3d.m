function [Xhigh, xArea] = fitEllipse3d(a,N)
% fitEllipse3d fits an ellipse to 3 pts in 3d space

if nargin < 2;N = 20;end

% center data and project to low dimension space
aMean = mean(a,1);
aCenter = a-aMean;
[u, s, v] = svd(aCenter);
aLowD = u(:,1:2)*s(1:2,1:2);

% fit ellipse
[A , C] = MinVolEllipse(aLowD', .01);

% get the major and minor axes
[U, D, V] = svd(A);
majAxis = 1/sqrt(D(1,1));
minAxis = 1/sqrt(D(2,2));
xArea = pi*majAxis*minAxis;
theta = [0:1/N:2*pi+1/N];

% Parametric equation of the ellipse
state(1,:) = majAxis*cos(theta);
state(2,:) = minAxis*sin(theta);

% Coordinate transform
X = V * state;
X(1,:) = X(1,:) + C(1);
X(2,:) = X(2,:) + C(2);

% project to higher dimension
XhighCenter = X'*v(:,1:2)';
Xhigh = XhighCenter + aMean;

% plot check
if 1==0
    figure, plot3(a(:,1),a(:,2),a(:,3),'b.'), axis equal, hold on
    plot3(Xhigh(:,1),Xhigh(:,2),Xhigh(:,3),'r-')
end

end