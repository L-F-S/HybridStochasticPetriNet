% Reproduction of  budding yeast model by Tyson and Novak, 2002.
%
% Author: Lorenzo Federico Signorini
%

function xdt= Tyson_Novak_det (t, x)

% parameters:

mu = 0.01;   
ms = 10.0;
k1 = 0.04;   
k21 = 0.04;
k22 = 1;    
k23 = 1;
k31 = 1;
k32 = 10;
J_three = 0.04;
k41 = 35; % called k4 in original papers (both Tyson and Novak 2001 and Mura 2008)
k42 = 2;  % called k4' in original paper (both Tyson and Novak 2001 and Mura 2008)
J_four = 0.04;
k51 = 0.005;
k52 = 0.2;
n_exp = 4;
J_five = 0.3;
k6 = 0.1;
k7 = 1;
J_seven = 0.001;
k8 = 0.5;
J_eight = 0.001;
k9 = 0.1;
k10 = 0.02;
k11 = 1;  % CKI- mutant: k11 = 0 
k121 = 0.2;
k122 = 50;
k123 = 100;
k131 = 0;
k132 = 1;  %SK- mutant: k131, k132 = 0 (k131 è già = 0) funge
k14 = 1;
k151 = 1.5;
k152 = 0.05;
J_fifteen = 0.01;
k161 = 1;
k162 = 3;
J_sixteen = 0.01;
keq = 1000;

% variables:

m = x(1);
CycBt = x(2);
Cdh1a= x(3);
Cdc20t = x(4);
Cdc20a = x(5);
IEP = x(6);
CKIt = x(7);
SK = x(8);
TF = x(9);
Cyc_B = CycBt - 2*CycBt*CKIt/((CycBt + CKIt + keq^-1)+sqrt((CycBt + CKIt + keq^-1)^2 -4*CycBt*CKIt));

%%%%%%%%%%%%%%%%%%%%%%%%%
%  DIFFERENTIAL EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%

% MASS EQUATION:

m_dt = mu*m*(1-m/ms);

% CELL NETWORK EQUATIONS:

CycBt_dt = k1 - (k21 + k22*Cdh1a + k23*Cdc20a)*CycBt;
Cdh1a_dt = (k31 + k32*Cdc20a)*(1-Cdh1a)/(J_three +1 -Cdh1a) - (k41*m*Cyc_B + k42*SK)*Cdh1a/(J_four+ Cdh1a);
Cdc20t_dt = k51 + k52*(m*Cyc_B)^n_exp/(J_five^n_exp + (m*Cyc_B)^n_exp)-k6*Cdc20t;
Cdc20a_dt = (k7*IEP*(Cdc20t-Cdc20a))/(J_seven+Cdc20t-Cdc20a)-k8*Cdc20a/(J_eight+Cdc20a)-k6*Cdc20a;
IEP_dt = k9*m*Cyc_B*(1-IEP)-k10*IEP;
CKIt_dt = k11 - (k121+k122*SK+k123*m*Cyc_B)*CKIt;
SK_dt = k131 + k132*TF -k14*SK;
TF_dt = (k151*m +k152*SK)*(1-TF)/(J_fifteen+1-TF)-(k161+k162*m*Cyc_B)*TF/(J_sixteen+TF); 
%------------------------------------------------------------------------------------------------------------
xdt = [m_dt;CycBt_dt;Cdh1a_dt;Cdc20t_dt;Cdc20a_dt;IEP_dt;CKIt_dt;SK_dt;TF_dt]; 
end
