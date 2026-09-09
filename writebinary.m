function writebinary(dst, mat)
[nr, nc] = size(mat);
fd = fopen(dst, 'wb');
fwrite(fd, nr, 'int');
fwrite(fd, nc, 'int');
fwrite(fd, mat', 'double');
fclose(fd);
end
%dstÊÇÊä³öÂ·¾¶
% writebinary('E:\wy\23_ronghe\Output_011\wrapped_phase_short_deflat.bin',wrapped_phase_short_deflat);