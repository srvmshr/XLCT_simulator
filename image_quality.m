R = 225;
vcm = 0.01;
v=[500 500 10];

[X,Y]=meshgrid( linspace(-(v(1)-1)/2,(v(1)-1)/2,v(1)), linspace(-(v(2)-1)/2,(v(2)-1)/2,v(2)));
Z = linspace(-(v(3)-1)/2,(v(3)-1)/2,v(3));

VV=(X.^2 + Y.^2 <= R^2);
h = fspecial('gaussian',5,1);
VV = imfilter(VV,h,'same');                 %%smooth out the jagged edges
W = repmat(VV,[1 1 v(3)]);

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

phantom = imcomplement(O(:,:,5));
lasimg = imread('lasimg.png');

figure(1); imshow(phantom);
figure(2); imshow(lasimg);

sum1 = 0;
sum2 = 0;

for i = 1:500
    for j = 1:500
        sum1 = sum1 + phantom(i, j);
        sum2 = sum2 + double(lasimg(i, j) / 256);
       
    end
end

abar = sum1 / (250000);
bbar = sum2 / (250000);

var1 = 0;
var2 = 0;
var12 = 0;

for i = 1:500
    for j = 1:500
        var1 =  (var1 + (phantom(i, j) - abar) ^ 2);
        var2 = (var2 + (double(lasimg(i, j) / 256) - bbar) ^ 2);
        
        var12 = (var12 + (phantom(i, j) - abar) * ...
        (double(lasimg(i, j) / 256) - bbar));
    end
end

sig_a = sqrt(var1 / 249999);
sig_b = sqrt(var2 / 249999);
sig_ab = var12 / 249999;

Q0 = (sig_ab / (sig_a * sig_b)) * ...
    (2 * abar * bbar / (abar ^ 2 + bbar ^ 2)) * ...
    (2 * sig_a * sig_b / (sig_a ^ 2 + sig_b ^ 2));

Q0