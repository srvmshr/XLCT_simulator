function [O] = gen_phan(name,R)

%% 1 voxel = 0.01 cm
vcm = 0.01;
v=[500 500 10];

if (nargin==1)
 R = 225;
end

%% Water phantom

[X,Y]=meshgrid( linspace(-(v(1)-1)/2,(v(1)-1)/2,v(1)), linspace(-(v(2)-1)/2,(v(2)-1)/2,v(2)));
Z = linspace(-(v(3)-1)/2,(v(3)-1)/2,v(3));

VV=(X.^2 + Y.^2 <= R^2);
h = fspecial('gaussian',5,1);
VV = imfilter(VV,h,'same');                 %%smooth out the jagged edges
W = repmat(VV,[1 1 v(3)]);

if (name=='W')
O = double(W);
end



%% Sphere Phantom
if (name=='S')
th = [0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
Si = zeros(5,6);
for i=1:6                                   %% spheres centers
  Si(1:2,i) = [cos(th(i)) sin(th(i)); -sin(th(i)) cos(th(i))] * [0;125];
end

% sphere z coordinate
Si(3,:) = 0;

% sphere radii
Si(4,:) = 5*[0.25, 0.5, 1, 2, 4, 8];

O = zeros(v);

for i=1:6                                   %%Create actual spheres
    for j=1:v(3)
        O(:,:,j) = O(:,:,j) + ((X-Si(1,i)).^2 + (Y-Si(2,i)).^2 + (Z(j)-Si(3,i)).^2 <= Si(4,i)^2);
    end
end

O = O + 0.1*W;
end

%% Contrast resolution
if (name=='C')
th = [0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
Si = zeros(4,6);
for i=1:6
  Si(1:2,i) = [cos(th(i)) sin(th(i)); -sin(th(i)) cos(th(i))] * [0;125];
end

% sphere z coordinate
Si(3,:) = 0;

% sphere radii
Si(4,:) = 5*[2, 2, 2, 2, 2, 2];

% sphere activity              %%Just change activity by another attribute
Si(5,:) = [2, 1, 0.75, 0.5, 0.25, 0.1];


O = zeros(v);

for i=1:6
    for j=1:v(3)
        O(:,:,j) = O(:,:,j) + Si(5,i)*((X-Si(1,i)).^2 + (Y-Si(2,i)).^2 + (Z(j)-Si(3,i)).^2 <= Si(4,i)^2);
    end
end

%  no background?
%O = O + 0.1*W;
end

%% calculate sensitivity map
if (name=='A')
load mc321.out
mc321(:,3) = mc321(:,3)/mc321(1,3);
%mc321(:,3) = 0.5*mc321(:,3)/mc321(1,3); %% 0.5W appears source power
mua=0.17;
musp=10.84;
D=1/(3*(mua+musp));
mueff = sqrt(mua/D);   %cm-1

nt = 500;
nb = 100;

E = zeros(1,nb);
RD = 2.25; % cm
th = linspace(0,(1-1/nt)*2*pi,nt);
b = linspace(0,RD,nb);

for i=1:nt
  d = sqrt( (RD*cos(th(i))-b).^2 + (RD*sin(th(i))).^2 );
  phi = interp1(mc321(:,1),mc321(:,3),d,'cubic');
  E = E + phi*2*pi*RD/nt;
end


%keyboard
d = sqrt( X.^2 + Y.^2 )*vcm;
O = interp1(b,E,d,'cubic');
O(d>RD) = 0;
end

%function ends
end
                                            




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (name=='D')
% xs=[];
% ys=[];
% rs=[];
% 
% Rm = 225; % max radius
% c = 30;  % separation
% R30 = [cos(30*pi/180) sin(30*pi/180) ; -sin(30*pi/180) cos(30*pi/180) ];
% N = 30;
% 
% %  half-dimensions of phantom:
% pxlen=250;
% pylen=250;
% pzlen=50;
% 
% %radius
% rr = [0.25 0.5 0.75 1.0 1.5 2]*10;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 1st quadrant
% an = 0;
% anr = an * pi /180;
% r = rr(1);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 2nd quadrant
% an = 60;
% anr = an * pi /180;
% r = rr(2);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 3rd quadrant
% an = 120;
% anr = an * pi /180;
% r = rr(3);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 4th quadrant
% an = 180;
% anr = an * pi /180;
% r = rr(4);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 5th quadrant
% an = 240;
% anr = an * pi /180;
% r = rr(5);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%% 6th quadrant
% an = 300;
% anr = an * pi /180;
% r = rr(6);
% sp = 4*r;
% 
% R = [cos(anr) sin(anr) ; -sin(anr) cos(anr) ];
% 
% for j=0:(N-1)
% for i=0:(N-j)
% 
% x = (i+j*0.5)*sp;
% y = j*0.5*sqrt(3)*sp;
% 
% xy = R30*[x ; y];
% xy2 = R*[xy(1)+c+r ; xy(2)];
% 
% xs = [xs xy2(1)];
% ys = [ys xy2(2)];
% rs = [rs r];
% 
% end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%% Clipping
% 
% r2 = xs.^2+ys.^2;
% idx = find(r2<(Rm-rs).^2);
% xs = xs(idx);
% ys = ys(idx);
% rs = rs(idx)*2;
% L = length(idx);
% 
% 
% 
% DR=zeros(v(1:2));
% 
% for i=1:L
%    DR = DR + ( (X-xs(i)).^2 + (Y-ys(i)).^2 < (rs(i)/2)^2 );
% end
% 
% O = double(repmat(DR,[1 1 v(3)]));

%end

