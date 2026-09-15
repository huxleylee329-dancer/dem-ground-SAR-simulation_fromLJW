function out = SAR_simulation_temporal_review(projectDir,cfg)
% Review-only example. Copy this file into your project when ready.
% Uses existing generateSLC_ground, registration_pixel, regis_subpixel.
% Four acquisitions: A1, B1, A2, B2. No raw echoes are generated.
% All noise levels are scenario assumptions, NOT calibrated device specs.
if nargin<1 || isempty(projectDir), projectDir=pwd; end
if nargin<2, cfg=struct; end
cfg=defaults(cfg);
addpath(projectDir);
assert(nargout('regis_subpixel')==4,'Use the four-output block registration version.');
lambda=299792458/16.7e9;
rho_r=299792458/(2*1e9); rho_u=lambda/(2*1.4);
dr=rho_r/2; du=rho_u/8; K=4*pi/lambda;
positions=[0 -3000 20;0 -3000 21.4];
x=(-160+cfg.dx/2):cfg.dx:(160-cfg.dx/2);
y=(-cfg.halfY+cfg.dy/2):cfg.dy:(cfg.halfY-cfg.dy/2);
[X,Y]=meshgrid(x,y);
h1=100+0.3*(Y+5)+4*exp(-(X/18).^2); % Preserve current project terrain.
q=sqrt(((X-cfg.patchX)/cfg.radiusX).^2+((Y-cfg.patchY)/cfg.radiusY).^2);
patch=zeros(size(q)); inside=q<1;
patch(inside)=0.5*(1+cos(pi*q(inside))); % Smooth, exactly zero outside ellipse.
dh=cfg.dhMax*patch; h2=h1+dh;
[gx,gy]=gradient(h1,cfg.dx,cfg.dy);
area=cfg.dx*cfg.dy*sqrt(1+gx.^2+gy.^2);
clear gx gy q inside
rng(cfg.seed);
alpha1=sqrt(area).*exp(1j*2*pi*rand(size(h1))); % Same t1 model as original.
rho=cfg.rhoStable+(cfg.rhoChanged-cfg.rhoStable)*patch;
innovation=sqrt(area).*(randn(size(h1))+1j*randn(size(h1)))/sqrt(2);
alpha2=rho.*alpha1+sqrt(max(0,1-rho.^2)).*innovation;
% Preserve facet scattering strength for small deformation.
% One common alpha per epoch, shared by both spatial observations.
clear innovation rho area
[R1,~]=geometry(X,Y,h1,positions(1,:));
[R2,~]=geometry(X,Y,h2,positions(1,:));
% Difference of squares avoids subtraction of two almost equal long ranges.
deltaR=dh.*(2*(h1-positions(1,3))+dh)./(R1+R2);
Rref=mean(R1(:));
truth.dh=dh; truth.deltaR_A=deltaR; truth.phi_A=K*deltaR;
truth.patch=patch; truth.x=x; truth.y=y;
fprintf('Scatterers per image: %d; time baseline: %.2f hours\n',numel(h1),cfg.dtSeconds/3600);
fprintf('Height change [mm]: %.6f .. %.6f\n',1e3*min(dh(:)),1e3*max(dh(:)));
fprintf('True temporal LOS [mm]: %.6f .. %.6f; phase [rad]: %.6f .. %.6f\n', ...
    1e3*min(deltaR(:)),1e3*max(deltaR(:)),min(truth.phi_A(:)),max(truth.phi_A(:)));
clear R1 R2 deltaR dh patch
% Differential PATH-AVERAGED refractivity, N=(n-1)*1e6.
% Low-order spatial screen, not an independent Gaussian phase per pixel.
% t1 is the differential atmospheric reference (screen zero).
dN=cfg.dNMean+cfg.dNX*(X/160)+cfg.dNY*(Y/cfg.halfY);
% Determine ONE common image grid from all four geometries.
rlo=Inf;rhi=-Inf;ulo=Inf;uhi=-Inf;
for epoch=1:2
    if epoch==1,h=h1;else,h=h2;end
    for p=1:2
        [R,U]=geometry(X,Y,h,positions(p,:));
        rlo=min(rlo,min(R(:)));rhi=max(rhi,max(R(:)));
        ulo=min(ulo,min(U(:)));uhi=max(uhi,max(U(:)));
    end
end
rAxis=(Rref+floor((rlo-Rref)/dr)*dr-4*rho_r):dr: ...
      (Rref+ceil((rhi-Rref)/dr)*dr+4*rho_r);
uAxis=(floor(ulo/du)*du-4*rho_u):du:(ceil(uhi/du)*du+4*rho_u);
images=cell(1,4); clean=cell(1,4); noiseSNR=nan(1,4);
t=tic;
for epoch=1:2
    if epoch==1,h=h1;alpha=alpha1;else,h=h2;alpha=alpha2;end
    for p=1:2
        k=2*(epoch-1)+p;
        [R,U]=geometry(X,Y,h,positions(p,:));
        if epoch==2
            % e_atm = R * dN * 1e-6 is a one-way equivalent path change [m].
            a=alpha.*exp(-1j*K*(R.*dN*1e-6));
        else
            a=alpha;
        end
        clean{k}=generateSLC_ground(R,U,a,lambda,rAxis,uAxis,rho_r,rho_u,Rref);
        fprintf('Generated acquisition %d/4; cumulative %.2f s\n',k,toc(t));
    end
end
% SNR is defined over the same reference echo ROI, excluding empty padding.
% This threshold-based ROI is an explicit convention, not calibrated sigma0.
roi=abs(clean{1})>0.05*max(abs(clean{1}(:)));
Psignal=mean(abs(clean{1}(roi)).^2);
Pnoise=Psignal*10^(-cfg.snrDb/10);
for k=1:4
    n=focused_noise(size(clean{k}),rho_r/dr,rho_u/du,Pnoise);
    images{k}=clean{k}+n;
    if Pnoise>0
        noiseSNR(k)=10*log10(Psignal/mean(abs(n(roi)).^2));
    else
        noiseSNR(k)=Inf;
    end
end
out.cfg=cfg;out.lambda=lambda;out.rAxis=rAxis;out.uAxis=uAxis;
out.Rref=Rref;out.positions=positions;out.truth=truth;
out.images=images;out.noiseSNR=noiseSNR;
out.cleanTemporalA=angle(clean{1}.*conj(clean{3}));
out.cleanTemporalB=angle(clean{2}.*conj(clean{4}));
if cfg.keepClean,out.clean=clean;end
if ~isempty(cfg.checkpointDir)
 save(fullfile(cfg.checkpointDir,'acquisitions_checkpoint.mat'),'out','-v7.3');
 fprintf('Four-image checkpoint saved.\n');
end
clear clean alpha1 alpha2 alpha a R U h h1 h2 X Y dN
fprintf('Measured reference-ROI noise SNR [dB]: %s\n',mat2str(noiseSNR,4));
% Registration always uses the images, not a terrain-derived warp.
if cfg.runPairs
    fprintf('Registering temporal pair A1/A2...\n');
    out.temporalA=register_pair(images{1},images{3},rAxis,uAxis,lambda,true);
    fprintf('Registering temporal pair B1/B2...\n');
    out.temporalB=register_pair(images{2},images{4},rAxis,uAxis,lambda,true);
    if cfg.runSpatial
        fprintf('Registering spatial pair A1/B1...\n');
        out.spatial1=register_pair(images{1},images{2},rAxis,uAxis,lambda,false);
        fprintf('Registering spatial pair A2/B2...\n');
        out.spatial2=register_pair(images{3},images{4},rAxis,uAxis,lambda,false);
    end
end
if cfg.plot
    figure;
    names={'A1: main t1','B1: slave t1','A2: main t2','B2: slave t2'};
    for k=1:4
        ax=subplot(2,2,k);imagesc(ax,abs(images{k}));colormap(ax,gray(256));
        xlabel('Range pixel');ylabel('Direction pixel');title(names{k});colorbar;
    end
    figure;
    subplot(1,2,1);imagesc(truth.x,truth.y,1e3*truth.dh);axis xy;colorbar;
    xlabel('Ground X [m]');ylabel('Ground Y [m]');title('Imposed vertical change [mm]');
    subplot(1,2,2);imagesc(truth.x,truth.y,1e3*truth.deltaR_A);axis xy;colorbar;
    xlabel('Ground X [m]');ylabel('Ground Y [m]');title('True temporal LOS change [mm]');
    if cfg.runPairs
        figure;
        subplot(1,3,1);show_masked(out.temporalA.phase,out.temporalA.valid);
        caxis([-pi pi]);title('A1 * conj(A2): raw temporal phase [rad]');
        subplot(1,3,2);show_masked(out.temporalA.losRelativeMM,out.temporalA.valid);
        title('LOS change, stable-area referenced [mm]');
        subplot(1,3,3);show_masked(out.temporalA.coherence,out.temporalA.valid);
        caxis([0 1]);title('Local coherence (correlated, oversampled pixels)');
    end
end
end

function cfg=defaults(cfg)
d=struct('dx',0.5,'dy',0.01,'halfY',75,'patchX',50,'patchY',15, ...
    'radiusX',70,'radiusY',35,'dhMax',-0.020,'dtSeconds',86400, ...
    'rhoStable',0.995,'rhoChanged',0.98,'snrDb',25, ...
    'dNMean',0.2,'dNX',0.05,'dNY',0.05,'seed',20260915, ...
    'runPairs',true,'runSpatial',true,'plot',true,'keepClean',false,'checkpointDir','');
names=fieldnames(d);
for k=1:numel(names)
    if ~isfield(cfg,names{k}),cfg.(names{k})=d.(names{k});end
end
assert(cfg.dx>0&&cfg.dy>0&&cfg.halfY>0);
assert(cfg.rhoStable>=0&&cfg.rhoStable<=1&&cfg.rhoChanged>=0&&cfg.rhoChanged<=1);
assert(cfg.radiusX>0&&cfg.radiusY>0);
end

function [R,U]=geometry(X,Y,H,pos)
R=sqrt((X-pos(1)).^2+(Y-pos(2)).^2+(H-pos(3)).^2);
U=(X-pos(1))./R;
end

function n=focused_noise(sz,osR,osU,Pnoise)
% Circular complex Gaussian noise filtered by the current ideal image PSF.
% Oversampled image noise is spatially correlated, not pixelwise white.
if Pnoise==0,n=complex(zeros(sz));return;end
hr=ceil(4*osR);hu=ceil(4*osU);
kr=nsinc((-hr:hr)/osR);ku=nsinc((-hu:hu)/osU).';
kr=kr/norm(kr);ku=ku/norm(ku);
z=(randn(sz(1)+2*hu,sz(2)+2*hr)+1j*randn(sz(1)+2*hu,sz(2)+2*hr))/sqrt(2);
n=sqrt(Pnoise)*conv2(ku,kr,z,'valid');
end

function y=nsinc(x)
y=ones(size(x));k=x~=0;y(k)=sin(pi*x(k))./(pi*x(k));
end

function result=register_pair(A,B,rAxis,uAxis,lambda,isTemporal)
[nr,nc]=size(A);
[Ac,Bc,mr,mc]=registration_pixel(A,B,nr,nc);
[nr,nc]=size(Ac);
[~,Br,registrationValid,info]=regis_subpixel(Ac,Bc,nr,nc);
r=rAxis((1:nc)+max(0,-mc));u=uAxis((1:nr)+max(0,-mr));
ifg=Ac.*conj(Br);
ampThreshold=0.05*max(abs(Ac(:)));
valid=registrationValid & isfinite(ifg) & ...
    abs(Ac)>ampThreshold & abs(Br)>ampThreshold;
phase=angle(ifg);phase(~valid)=NaN;
% Local coherence using a common finite-data support; zero-fill only for sums.
finite=registrationValid & isfinite(ifg);
At=Ac;Bt=Br;At(~finite)=0;Bt(~finite)=0;
win=ones(17,5);
num=conv2(At.*conj(Bt),win,'same');
den=sqrt(conv2(abs(At).^2,win,'same').*conv2(abs(Bt).^2,win,'same'));
coh=min(1,abs(num)./max(den,eps));
supported=conv2(double(finite),win,'same')>=0.8*numel(win);
coh(~supported|~finite)=NaN;
result.phase=phase;result.valid=valid;result.coherence=coh;
result.rAxis=r;result.uAxis=u;result.shift=[mr mc];result.info=info;
if isTemporal
    % This stable sector is deliberately outside the imposed moving patch.
    % It is a scenario-selected reference, not a geometry registration.
    stable=valid & repmat(u(:)<-0.035,1,nc) & coh>0.8;
    assert(nnz(stable)>=20,'Too few stable reference pixels; inspect coherence and scene coverage.');
    phiRef=angle(sum(ifg(stable)));
    referenced=angle(ifg*exp(-1j*phiRef));
    referenced(~valid)=NaN;
    result.referencePhase=phiRef;result.stableMask=stable;
    result.losRelativeMM=1e3*lambda/(4*pi)*referenced;
    % Principal-branch LOS only: do not use for changes beyond wrapping limits.
    % A constant reference removes common offset, NOT a full spatial APS.
    fprintf('Temporal pair shift=[%g,%g], accepted=%d, reference=%.6f rad\n', ...
        mr,mc,info.accepted,phiRef);
end
end

function show_masked(A,mask)
A(~mask)=0;
h=imagesc(A);set(h,'AlphaData',double(mask));set(gca,'Color','w');
xlabel('Range pixel');ylabel('Direction pixel');colorbar;
end

