function [P,v] = genlines(R,z,nr,na);

n = nr*na;

Q = zeros(2,nr);
Q(2,:) = -1.5*R;
Q(1,:) = linspace(-R,R,nr);
V=[0; 1];

A = linspace(0, 2*pi, na);

P = zeros(nr*na,3);
P(:,3) = z;

v = zeros(nr*na,2);

for a=1:na
    for r=1:nr
        R = [cos(A(a)) sin(A(a)); -sin(A(a)) cos(A(a))];
        P(r+(a-1)*nr,1:2) = R * Q(:,r);
        v(r+(a-1)*nr,:) = R * V;
    end
end
