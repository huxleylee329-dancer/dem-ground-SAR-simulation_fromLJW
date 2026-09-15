function resultDir=run_temporal_case()
% Copy this file, SAR_simulation_temporal_review.m and export_temporal_results.m
% into the existing project root. Run: resultDir=run_temporal_case;
projectDir=fileparts(mfilename('fullpath'));
resultDir=fullfile(projectDir,'output',['temporal_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(resultDir,'dir'),'Output directory already exists.');
mkdir(resultDir);mkdir(fullfile(resultDir,'images'));
oldVisible=get(groot,'defaultFigureVisible');
existingFigures=findall(groot,'Type','figure');
restore=onCleanup(@() set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
diary(fullfile(resultDir,'run_log.txt'));
diaryCleanup=onCleanup(@() diary('off')); %#ok<NASGU>
cfg=struct;
cfg.dx=0.5; cfg.dy=0.01; cfg.halfY=75;
cfg.patchX=50; cfg.patchY=15; cfg.radiusX=70; cfg.radiusY=35;
cfg.dhMax=-0.020; cfg.dtSeconds=24*3600;
cfg.snrDb=25; cfg.rhoStable=0.995; cfg.rhoChanged=0.98;
cfg.dNMean=0.2; cfg.dNX=0.05; cfg.dNY=0.05;
cfg.seed=20260915;
cfg.plot=false;cfg.keepClean=true;cfg.runPairs=true;cfg.runSpatial=true;
cfg.checkpointDir=resultDir;
t=tic;
out=SAR_simulation_temporal_review(projectDir,cfg);
out.elapsedSeconds=toc(t);
save(fullfile(resultDir,'simulation_results.mat'),'out','-v7.3');
diagnosticFigures=setdiff(findall(groot,'Type','figure'),existingFigures);
export_temporal_results(out,resultDir,diagnosticFigures);
fprintf('Saved all results to: %s\n',resultDir);
end

