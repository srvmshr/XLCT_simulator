function [D] = raytrace(b,v,w,A,E,vcm);
% raytrace the beam
% b: beam origin position (3D) in cm
% v: horizontal beam direction (x-y). Unit vector, unitless
% w: Beam dimension (w-h) in cm
% A: Attenuation map, binary
% E: Energy (spectrum) of the beam, in keV
% vcm: size of the voxel in cm
% D: Dose deposited by the beam, in keV

load water;

is = size(A);
D2 = zeros(is(1:2));

% half fov size
%s = vcm*(is(1:2)-1)/2;
s = vcm*(is(1:2))/2;

I = compinter(b(1:2),v,0.999*s);

sl = floor(b(3)+(is(3)-1)/2)+1;

if (length(I)<1)
    disp('Beam does not intersect');
    return;
end

v = v / norm(v,2);

% find main direction
n1 = [abs(v(1))>abs(v(2)),abs(v(1))<=abs(v(2))];

% complementary direction
n2 = 1 - n1;

l = floor(w(1)/(vcm*abs(v(n1==1))))+1;

% calculate attenuation profile
J = floor((I+s)/vcm)+1;
step = 0.001;   %in cm
mu = 0;m=-1000;
k=0;
while (all(J>0) && all(J<=is(1:2)) )
    mu = [mu A(J(1),J(2),sl)];
    m = [m k*step];
    k = k+1;
    J = floor( (I + s + k*step*v)/vcm)+1;
end
mul = interp1(water(:,1),water(:,2),E/1000,'cubic');
mue = interp1(water(:,1),water(:,3),E/1000,'cubic');
cmu = cumsum(mu)*step*mul; % cumulative

% back a little
P = I - v*l*vcm/abs(v(n1==1));
V = floor((P+s)/vcm)+1;
P = (V-0.5)*vcm-s;

% beam profile, in cm
r = linspace(-4*w(1),4*w(1),200);
B = double(abs(r)<=w(1));
B = conv(B,gausswin(12),'same');


f1 = 0;
f2 = 0;
while (f1>0 || f2==0)
    f1 = 0;
    for i=-l:l
        Q = P+i*n2*vcm;               % position in cm
        V = floor((Q+s)/vcm)+1;       % voxel position
        R = (V-0.5)*vcm-s;            % bin to the nearest voxel in cm
        IR = R-I;
        ds = dot(IR,v);         % distance in cm
        dr = sqrt(abs(norm(IR,2)^2-ds^2));
        if (all(V>0) && all(V<=is(1:2)))
            id1 = find(r<=dr,1,'last');
            if (id1<length(r))
                BB = ( B(id1)*(r(id1+1)-dr) + B(id1+1)*(dr-r(id1)) ) / (r(id1+1)-r(id1));
            else
                BB = B(id1);
            end
                
            id2 = find(m<=ds,1,'last');
            if (id2<length(m))
                CC = ( exp(-cmu(id2))*(m(id2+1)-ds) + exp(-cmu(id2+1))*(ds-m(id2)) ) / (m(id2+1)-m(id2));
            else
                CC = exp(-cmu(id2));
            end
            D2(V(1),V(2)) = BB*CC;

            %interp1(r,B,dr,'linear')*exp(-interp1(m,cmu,ds,'nearest','extrap'))*interp1(m,mu,ds,'nearest','extrap')*vcm;
            f1=1;
            f2=1;
        end
                 
    end
    P = P + vcm*v / abs(v(n1==1));
end

D = mue*vcm*repmat(D2,[1 1 is(3)]).*A;
