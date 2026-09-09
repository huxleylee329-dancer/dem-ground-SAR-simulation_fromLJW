function mat = readbinary(src)
fd = fopen(src, 'rb');
nr = fread(fd, 1, 'int');
nc = fread(fd, 1, 'int');
mat = fread(fd, [nc, nr], 'double');
mat = mat';
fclose(fd);
end
%srcÊÇÊäÈëÂ·¾¶
%unwrapped_phase_long = readbinary('E:\wy\23_ronghe\Output_013\unwrapped_phase_long.bin');