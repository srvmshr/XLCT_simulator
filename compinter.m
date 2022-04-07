function [I] = compinter( P, v, fov )
% calculates the intersection between the line defined by point P and
% vector v, and the rectangle defined by fov (2D)
% output I, intersection point. If I=[], no intersection exists

	l = ([fov -fov] - [P P]) ./ [v v];
	
    % Check which intersections are on the fov boundaries
    hit = zeros(1,4);
    for i=1:4
        hit(i) = all(abs(P + l(i)*v)<=1.001*fov);            
    end
    
    g = l(l>0 & hit);
    
    % find closest intersection
    [junk,idx]=min(g);

    if (length(idx) == 1)
        I =  P + g(idx)*v;
    else
        I = [];
    end
    
    
