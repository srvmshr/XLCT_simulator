function I = recon(P,v,S,it,display,W,E,vcm,A,bw)
% Reconstruct the sinogram S, defined according to P and v
% Based on ML-EM

[n,junk] = size(P);
s = size(W);
sens = zeros(s);
I = zeros(s(1),s(2),it);
J = ones(s);

[X,Y]=meshgrid( linspace(-(s(1)-1)/2,(s(1)-1)/2,s(1)), linspace(-(s(2)-1)/2,(s(2)-1)/2,s(2)));
R = 225;
VV=(X.^2 + Y.^2 <= R^2);

%sensitivity map
 parfor j=1:n
     sens = sens + A.*raytrace_recon(P(j,:),v(j,:),bw,W,E,vcm);
 end
 imwrite(1-sens/max(sens(:)),'sens.png');
 
%sens = ones(s);

for i=1:it
   if (display)
       imagesc(J);axis off;axis image;colormap(1-gray);
       drawnow;
   end
   U = zeros(s);
   parfor j=1:n
       M = A.*raytrace_recon(P(j,:),v(j,:),bw,VV,E,vcm);
       f = sum(J(:).*M(:));
       if (f>0)
           U = U + S(j)/f*M;
       end
   end
   
   idx = find(sens>eps);
   J(idx) = J(idx).*U(idx)./sens(idx);
   J(sens<=eps) = 0;
   J(VV<0.1) = 0;
   
   I(:,:,i) = J;
   i                                            %%display the MLEM iteration count
   imwrite(1-J/max(J(:)),'lasimg.png');
end

