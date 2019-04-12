function[value,isterminal,direction] = Tyson_Novak_det_variables_checkpoints(t,x)
keq=1000;
%% variables
m = x(1);
CycBt= x(2);
CKIt = x(7);
Cyc_B = CycBt - 2*CycBt*CKIt/(CycBt + CKIt + keq^-1)*sqrt((CycBt + CKIt + keq^-1)^2 -4*CycBt*CKIt);

%parameter
Cyc_B_TH = 0.1; %Cyc_B threshold for cell division (Tyson and Novak, 2001)

%CELL DIVISION TRIGGER
value = Cyc_B*m - Cyc_B_TH*m;
isterminal = 1;
direction = -1;
end