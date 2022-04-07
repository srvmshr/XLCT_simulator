clear all
close all
clc
matlabpool open

W=gen_phan('W');    %Water Phantom
OA=gen_phan('A');   %Optical attenuation map
D=gen_phan('S');    %Sphere phantom


na = 30;            % angular bins
vcm = 0.01;         % voxel size
nr = 500;            % radial bins
E = 100;            % energy

bw = [5.0 / nr, 5.0 / nr];

%recon
it = 30;                                                 %%iteration
vcmr = vcm;                                              %%vcm set to vcmr
%vcmr = 0.01;                                            %%original 0.05

% resample water phantom
W2 = imresize(W(:,:,5),1,'bilinear');                    %%originally 0.2
OA2 =imresize(OA,1,'bilinear');

[P,v] = genlines(2.4,0,nr,na);

% simulation
sino = runsim(P,v,vcm,W,D,E,bw,OA);            %%originally 0.05

% reconstruction
I = recon(P,v,sino,it,0,W2,E,vcmr,OA2,bw);
matlabpool close

%%Lines to change
% testrun.m line 13,16
