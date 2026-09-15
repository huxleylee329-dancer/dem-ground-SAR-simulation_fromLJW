function export_temporal_results(out,resultDir,diagnosticFigures)
% Save reproducible MATLAB figures and PNG images; numerical arrays stay in MAT.
imageDir=fullfile(resultDir,'images');
if ~exist(imageDir,'dir'),mkdir(imageDir);end
set(groot,'defaultFigureVisible','off','defaultAxesFontName','Microsoft YaHei', ...
    'defaultTextFontName','Microsoft YaHei');
% These are the diagnostic figures created by registration_pixel/real_coherent.
if nargin<3
    figs=findall(groot,'Type','figure'); % Standalone batch run.
else
    figs=diagnosticFigures;
end
if ~isempty(figs)
    [~,order]=sort([figs.Number]);figs=figs(order);
    labels={'时间配准 A1-A2','时间配准 B1-B2','空间配准 A1-B1','空间配准 A2-B2'};
    for k=1:numel(figs)
        if k<=numel(labels),name=labels{k};else,name=sprintf('诊断图 %d',k);end
        ax=findall(figs(k),'Type','axes');
        if ~isempty(ax)
            title(ax(1),[name '：互相关峰值图']);
            xlabel(ax(1),'相关矩阵列号');ylabel(ax(1),'相关矩阵行号');
        end
        write_figure(figs(k),imageDir,sprintf('00_correlation_%d',k));
    end
end
% Terrain plots are decimated FOR DISPLAY ONLY, original data are in out.truth.
ix=1:max(1,round(numel(out.truth.x)/320)):numel(out.truth.x);
iy=1:max(1,round(numel(out.truth.y)/300)):numel(out.truth.y);
x=out.truth.x(ix);y=out.truth.y(iy);[X,Y]=meshgrid(x,y);
H=100+.3*(Y+5)+4*exp(-(X/18).^2);
dh=out.truth.dh(iy,ix);
f=newfig([1300 650]);tiledlayout(1,2,'TileSpacing','compact');
nexttile;surf(X,Y,H,'EdgeColor','none');view(35,30);colorbar;
xlabel('地面 X / m');ylabel('地面 Y / m');zlabel('高度 / m');
title('t1 地形：原坡面与山脊');caxis([79 128]);
nexttile;surf(X,Y,H+dh,'EdgeColor','none');view(35,30);colorbar;
xlabel('地面 X / m');ylabel('地面 Y / m');zlabel('高度 / m');
title('t2 地形：局部最多下沉 20 mm');caxis([79 128]);
sgtitle('两时刻地形：厘米变化在米级高度图中不明显，需查看差值图');
write_figure(f,imageDir,'01_terrain_two_epochs');
f=newfig([1300 480]);tiledlayout(1,3,'TileSpacing','compact');
nexttile;imagesc(x,y,1e3*dh);axis xy;colorbar;caxis([min(out.truth.dh(:))*1e3 0]);
xlabel('地面 X / m');ylabel('地面 Y / m');title('竖直形变真值 / mm');
nexttile;imagesc(x,y,1e3*out.truth.deltaR_A(iy,ix));axis xy;colorbar;
xlabel('地面 X / m');ylabel('地面 Y / m');title('主位置视线位移真值 / mm');
nexttile;imagesc(x,y,out.truth.phi_A(iy,ix));axis xy;colorbar;
xlabel('地面 X / m');ylabel('地面 Y / m');title('主位置时间相位真值 / rad');
sgtitle('形变真值：地面坐标；负值表示下沉或距离缩短');
write_figure(f,imageDir,'02_deformation_truth');
names={'A1：t1 主图像','B1：t1 辅图像','A2：t2 主图像','B2：t2 辅图像'};
stems={'03_A1_amplitude','04_B1_amplitude','05_A2_amplitude','06_B2_amplitude'};
ampMax=0;
for k=1:4,ampMax=max(ampMax,max(abs(out.images{k}(:))));end
for k=1:4
    f=newfig([1100 550]);ax=axes(f);imagesc(ax,abs(out.images{k}));
    colormap(ax,gray(256));caxis(ax,[0 ampMax]);colorbar;
    xlabel('距离向像素编号（列）');ylabel('方向向像素编号（行）');
    title([names{k} '，线性幅度（含噪声，统一色标）']);
    write_figure(f,imageDir,stems{k});
end
f=newfig([1300 750]);tiledlayout(2,2,'TileSpacing','compact');
for k=1:4
    ax=nexttile;imagesc(ax,abs(out.images{k}));colormap(ax,gray(256));
    caxis(ax,[0 ampMax]);colorbar;title(names{k});
    xlabel('距离向像素编号');ylabel('方向向像素编号');
end
sgtitle('四幅复数图像的线性幅度：统一色标，无 dB 转换');
write_figure(f,imageDir,'07_four_amplitudes');
if isfield(out,'spatial1')
    save_map(out.spatial1.phase,out.spatial1.valid,[-pi pi], ...
        't1 空间干涉相位 A1 × conj(B1) / rad','08_spatial_t1',imageDir);
    save_map(out.spatial2.phase,out.spatial2.valid,[-pi pi], ...
        't2 空间干涉相位 A2 × conj(B2) / rad','09_spatial_t2',imageDir);
end
if isfield(out,'temporalA')
    a=out.temporalA;b=out.temporalB;
    save_map(a.phase,a.valid,[-pi pi], ...
        '主位置时间干涉 A1 × conj(A2) / rad（含大气）','10_temporal_A_raw',imageDir);
    save_map(b.phase,b.valid,[-pi pi], ...
        '辅位置时间干涉 B1 × conj(B2) / rad（含大气）','11_temporal_B_raw',imageDir);
    save_map(a.losRelativeMM,a.valid,[-1 1], ...
        '主位置相对视线位移 / mm（稳定区参考，显示 -1～1 mm）','12_temporal_A_LOS',imageDir);
    save_map(b.losRelativeMM,b.valid,[-1 1], ...
        '辅位置相对视线位移 / mm（稳定区参考，显示 -1～1 mm）','13_temporal_B_LOS',imageDir);
    save_map(a.coherence,a.valid & isfinite(a.coherence),[0 1], ...
        '主位置时间干涉局部相干性（17×5 像素窗口）','14_temporal_A_coherence',imageDir);
    save_map(b.coherence,b.valid & isfinite(b.coherence),[0 1], ...
        '辅位置时间干涉局部相干性（17×5 像素窗口）','15_temporal_B_coherence',imageDir);
    save_map(a.info.dCol,a.valid,[], ...
        '主位置时间配准：残余列偏移 / 像素','16_temporal_A_column_shift',imageDir);
    save_map(b.info.dCol,b.valid,[], ...
        '辅位置时间配准：残余列偏移 / 像素','17_temporal_B_column_shift',imageDir);
    f=newfig([1400 820]);tiledlayout(2,2,'TileSpacing','compact');
    ax=nexttile;plot_map(ax,a.phase,a.valid,[-pi pi]);title('A1/A2 原始时间相位 / rad（含大气）');
    ax=nexttile;plot_map(ax,a.losRelativeMM,a.valid,[-1 1]);title('稳定区参考后的视线位移 / mm');
    ax=nexttile;plot_map(ax,a.coherence,a.valid & isfinite(a.coherence),[0 1]);title('局部相干性');
    ax=nexttile;imagesc(ax,x,y,1e3*out.truth.deltaR_A(iy,ix));axis(ax,'xy');colorbar(ax);
    xlabel('地面 X / m');ylabel('地面 Y / m');title('地面网格上的视线位移真值 / mm');
    sgtitle('时间干涉结果与真值（坐标网格不同，不能直接逐像素相减）');
    write_figure(f,imageDir,'18_temporal_summary');
end
metrics=struct;
metrics.scatterers=numel(out.truth.dh);
metrics.gridSize=size(out.images{1});
metrics.noiseSNR=out.noiseSNR;
metrics.heightRangeMM=1e3*[min(out.truth.dh(:)),max(out.truth.dh(:))];
metrics.losTruthRangeMM=1e3*[min(out.truth.deltaR_A(:)),max(out.truth.deltaR_A(:))];
metrics.phaseTruthRange=[min(out.truth.phi_A(:)),max(out.truth.phi_A(:))];
if isfield(out,'elapsedSeconds'),metrics.elapsedSeconds=out.elapsedSeconds;end
if isfield(out,'temporalA')
    metrics.temporalA=pair_metrics(out.temporalA);
    metrics.temporalB=pair_metrics(out.temporalB);
end
if isfield(out,'spatial1')
    metrics.spatial1=pair_metrics(out.spatial1);
    metrics.spatial2=pair_metrics(out.spatial2);
end
fid=fopen(fullfile(resultDir,'metrics.json'),'w','n','UTF-8');
assert(fid>0);fprintf(fid,'%s',jsonencode(metrics));fclose(fid);
fprintf('Saved %d PNG files to %s\n',numel(dir(fullfile(imageDir,'*.png'))),imageDir);
end

function f=newfig(sz)
f=figure('Visible','off','Color','w','Position',[50 50 sz]);
end
function save_map(A,mask,limits,label,stem,folder)
f=newfig([1100 560]);ax=axes(f);plot_map(ax,A,mask,limits);title(ax,label);
write_figure(f,folder,stem);
end
function plot_map(ax,A,mask,limits)
mask=mask & isfinite(A);A(~mask)=0;
h=imagesc(ax,A);set(h,'AlphaData',double(mask));set(ax,'Color','w');colorbar(ax);
if ~isempty(limits),caxis(ax,limits);end
xlabel(ax,'距离向像素编号（列）');ylabel(ax,'方向向像素编号（行）');
end
function write_figure(f,folder,stem)
exportgraphics(f,fullfile(folder,[stem '.png']),'Resolution',160);
savefig(f,fullfile(folder,[stem '.fig']));
close(f);
end
function m=pair_metrics(p)
m.coarseShift=p.shift;m.acceptedBlocks=p.info.accepted;m.candidateBlocks=p.info.candidates;
m.validPixels=nnz(p.valid);
q=p.coherence(p.valid & isfinite(p.coherence));m.medianCoherence=median(q);
if isfield(p,'referencePhase')
 m.referencePhase=p.referencePhase;
 v=p.losRelativeMM(p.stableMask);
 m.stableLOSMeanMM=mean(v);m.stableLOSStdMM=std(v);
end
end

