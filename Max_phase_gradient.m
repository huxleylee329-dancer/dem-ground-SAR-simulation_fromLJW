function MPG = Max_phase_gradient(wrapped_phase)

[nr,nc] = size(wrapped_phase);
MPG = zeros(nr,nc);
delta = zeros(8,1);
for i = 2:nr-1
    for j = 2:nc-1
        delta(1) = abs(angle(exp(1j*(wrapped_phase(i-1,j-1) - wrapped_phase(i,j)))))/sqrt(2);
        delta(2) = abs(angle(exp(1j*(wrapped_phase(i,j-1) - wrapped_phase(i,j)))))/sqrt(1);
        delta(3) = abs(angle(exp(1j*(wrapped_phase(i+1,j-1) - wrapped_phase(i,j)))))/sqrt(2);
        delta(4) = abs(angle(exp(1j*(wrapped_phase(i+1,j) - wrapped_phase(i,j)))))/sqrt(1);
        delta(5) = abs(angle(exp(1j*(wrapped_phase(i+1,j+1) - wrapped_phase(i,j)))))/sqrt(2);
        delta(6) = abs(angle(exp(1j*(wrapped_phase(i,j+1) - wrapped_phase(i,j)))))/sqrt(1);
        delta(7) = abs(angle(exp(1j*(wrapped_phase(i-1,j+1) - wrapped_phase(i,j)))))/sqrt(2);
        delta(8) = abs(angle(exp(1j*(wrapped_phase(i-1,j) - wrapped_phase(i,j)))))/sqrt(1);
        MPG(i,j) = mean(delta);
    end
end
MPG(1,:) = MPG(2,:);
MPG(end,:) = MPG(end-1,:);
MPG(:,1) = MPG(:,2);
MPG(:,end) = MPG(:,end-1);

MPG(1,1) = MPG(2,2);
MPG(1,end) = MPG(2,end-1);
MPG(end,1) = MPG(end-1,2);
MPG(end,end) = MPG(end-1,end-1);

end


