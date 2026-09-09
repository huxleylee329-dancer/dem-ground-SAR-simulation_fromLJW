function SPD = Sum_of_phase_difference(wrapped_phase)

[nr,nc] = size(wrapped_phase);
APD = zeros(nr,nc);
delta = zeros(8,1);
for i = 2:nr-1
    for j = 2:nc-1
%         delta(1) = abs(angle(exp(1j*(wrapped_phase(i-1,j-1) - wrapped_phase(i,j)))));
%         delta(2) = abs(angle(exp(1j*(wrapped_phase(i,j-1) - wrapped_phase(i,j)))));
%         delta(3) = abs(angle(exp(1j*(wrapped_phase(i+1,j-1) - wrapped_phase(i,j)))));
%         delta(5) = abs(angle(exp(1j*(wrapped_phase(i+1,j+1) - wrapped_phase(i,j)))));
%         delta(6) = abs(angle(exp(1j*(wrapped_phase(i,j+1) - wrapped_phase(i,j)))));
%         delta(7) = abs(angle(exp(1j*(wrapped_phase(i-1,j+1) - wrapped_phase(i,j)))));
%         delta(8) = abs(angle(exp(1j*(wrapped_phase(i-1,j) - wrapped_phase(i,j)))));
        delta(1) = abs(wrapped_phase(i-1,j-1) - wrapped_phase(i,j));
        delta(2) = abs(wrapped_phase(i,j-1) - wrapped_phase(i,j));
        delta(3) = abs(wrapped_phase(i+1,j-1) - wrapped_phase(i,j));
        delta(5) = abs(wrapped_phase(i+1,j+1) - wrapped_phase(i,j));
        delta(6) = abs(wrapped_phase(i,j+1) - wrapped_phase(i,j));
        delta(7) = abs(wrapped_phase(i-1,j+1) - wrapped_phase(i,j));
        delta(8) = abs(wrapped_phase(i-1,j) - wrapped_phase(i,j));

        APD(i,j) = mean(delta);
    end
end
APD(1,:) = APD(2,:);
APD(end,:) = APD(end-1,:);
APD(:,1) = APD(:,2);
APD(:,end) = APD(:,end-1);

APD(1,1) = APD(2,2);
APD(1,end) = APD(2,end-1);
APD(end,1) = APD(end-1,2);
APD(end,end) = APD(end-1,end-1);

SPD = sum(sum(APD));

end


