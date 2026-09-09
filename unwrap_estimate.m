function [Gradiennt_jump] = unwrap_estimate(image_phase)

phi_threshold = 1*pi;
%跳变点是指那些相位梯度在【-π，π】之间的不连续点，因此在判断时如果绝对值大于π即为跳变点，而不是2π

[pr,pc] = size(image_phase);
image_left = image_phase(:,1:pc-1);
image_right = image_phase(:,2:pc);
dx = floor(abs(image_left-image_right)/phi_threshold);
image_up = image_phase(1:pr-1,:);
image_down = image_phase(2:pr,:);
dy = floor(abs(image_up-image_down)/phi_threshold);
dx(dx>1)=1;
dy(dy>1)=1;
dz1=dx(2:end,:)+dy(:,1:end-1);
dz2=dx(2:end,:)+dy(:,2:end);
dz3=dx(1:end-1,:)+dy(:,1:end-1);
dz4=dx(1:end-1,:)+dy(:,2:end);
dz = dz1 + dz2 + dz3 + dz4;
dz(dz>1)=1;

[nr_z,nr_c] = size(dz);
for i = 2:nr_z-1
    for j = 2:nr_c-1
%         if (dz(i-1,j) + dz(i+1,j) + dz(i,j-1) + dz(i,j-1)) > 2;
%如果用>2这一条件，按照循环顺序，会导致左侧存在三个跳变点的元素都被标记为跳变点，
%最终导致跳变点不断往右侧延伸，因此这一判断条件不太合理
        if (dz(i-1,j) + dz(i+1,j) + dz(i,j-1) + dz(i,j-1)) == 4;         
            dz(i,j) = 1;
        end
    end
end

[pr,pc,~] = find(dz == 1); 
Gradiennt_jump = nansum(sum(dz));%跳变点个数

figure;
imagesc(image_phase);
hold on
scatter(pc,pr,'r.');%跳变点

end