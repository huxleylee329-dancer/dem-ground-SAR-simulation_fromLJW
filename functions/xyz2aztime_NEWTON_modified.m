function azitime = xyz2aztime_NEWTON_modified(satellite_time,coef,xyz)
%------brustindex: xyz belongs to  which brust --------
% init_time =(datenum(strrep(xml.imageInfo.imageTime.start,'T',' '))...
%     - datenum(strrep(xml.imageInfo.imageTime.end,'T',' ')))+...
%     (xml.swathTiming.linesPerBurst/2)*xml.imageAnnotation.imageInformation.azimuthTimeInterval ;
init_time = (satellite_time(floor(end/2)) - satellite_time(1)) * 10000;
sol = 0;
xyz = xyz(:)';
azitime = init_time;
for iter = 0:30;
    [S_xyz] = getXyz(azitime,coef);
    [V_xyz] = getVel(azitime,coef);
    [A_xyz] = getAcc_diffVel(azitime,coef);%jiasudu huoqu
    D_xyz = xyz - S_xyz;
    sol = -(V_xyz*D_xyz')/(A_xyz*D_xyz'-V_xyz*V_xyz');
    azitime = azitime+sol;
    if (abs(sol)<1*exp(-10))
        break;
    end
end
end
