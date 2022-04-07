function [sino] = runsim(P,v,vcm,W,C,E,bw,OA)
% Runs a simulation for the beams described by P,V
% uses the water phantom W and nanophosphor distribution C
% the voxel size in cm is vcm
% the beam energy is E
% the beam width is bw
% The optical attenuation map is OA
% output is sino, the amount of light for each beam
% photon created assuming 60 Photon per keV

conc = 1e-6;            % concentration in g.cm-3
rho_solid_GOS = 7.44;   % density in g.cm-3
phperkv = 60 * conc / rho_solid_GOS;

n = size(P,1);

sino = zeros(n,1);

total_dose = 0;

% Main loop
parfor k=1:n

    D = raytrace(P(k,:),v(k,:),bw,W,E,vcm);

    % uncomment for live display
    imagesc(log(D(:,:,5)));axis off;axis image;drawnow;

    sD = sum(D(:));
    total_dose = total_dose+sD;

    vv = size(W);
    
    % light distribution
    L = phperkv*D.*C.*repmat(OA,[1 1 vv(3)]);

    sino(k)=sum(L(:));

end

%%Add poisson noise 5%

sino = sino + (0.05*sino-0.1*poissrnd(sino));   % 5p.c poisson noise added to simulate
                                     % Poisson fluctuation at detector.

% convert dose from keV to cGy
% Gy = J / Kg 
% keV * 1000 * 1.6e-19 -> J
tdg = (total_dose * 1000) * 1.6e-19 / (pi*2.25^2*0.1/1000) * 100; % cGy

% normalize to get 1 cGy total dose to water
sino = sino / tdg;