
% 29/12/2017, 01:40%    
% Reproduction of  budding yeast model by Tyson and Novak, 2002.
%
% Author: Lorenzo Federico Signorini
%

clc
clear all
close all

%%% INITIAL VALUES

% deterministic initial values will be convertedd to SPN initial values of tokens with the reverse of the formula
% they used to get #tokens: [species] = #species/a,
% where a=(AvgCellVolume*Navogadro*10^-6)^1 = 0.00236012 (10^-6 is to
% scale back from microMolar)

% Initial conditions
m_int = 0.704045; 
CycBt_int = 0.228559; 
Cdh1a_int= 0.011343;
Cdc20t_int = 0.056904;
Cdc20a_int = 2.26E-4;
IEP_int = 0.094007;
CKIt_int = 0.059228;
SK_int = 0.093081;
TF_int = 0.034886;


%parameters: specified inside function

%%%%%%%
%Solve ODEs:
%%%%%%%

%initial conitions:
% LITTLE DIFFERENCE : xint is a column vector now. works better for later outputs this
%way. (BEWARE THOUGH THAT the odefunc still wants a row vector output,
%but accepts this input, AND oode15s  outputs xe as column vector (see below)
xint = [m_int,CycBt_int,Cdh1a_int,Cdc20t_int,Cdc20a_int,IEP_int,CKIt_int,SK_int,TF_int];
%time
tinit = 0;
%initialize vectors to store all values of solutions and time
x = xint;
t = [tinit];
for i=1:6 %iterates the ode15s 13 times. Every iteration stops when Cyc_B_TH is reached.
options = odeset('RelTol',1e-10,'AbsTol', 1e-8,'Event',@Tyson_Novak_det_variables_checkpoints); 
%solver
[tempt,tempx,te,xe,ie] = ode15s(@Tyson_Novak_det,[tinit tinit+300],xint,options);
te %te = column vector of the times at which events occurred.
xe % xe = contains the solution value at each of the event times in te. 
ie %ie = contains indices into the vector returned by the event function.

%UPDATE INITIAL CONDITIONS
sz = size(xe); %DEBUGS
if sz(1) > 1
    display('xe contains more than one line, only keeping last one.')
    last_event = xe(sz(1),:) %a volte, quando le initial conditions sono quelle che triggerano l'evento, comunque parte e fa un giro della ode, ma salva comunque l'evento in una riga di xe, quindi prendere sempre l'ultima riga)
elseif sz(1)<1
    display('ERROR: empty array of final values')
    break
else
    last_event = xe;
end
last_event(1) = last_event(1)/2; % halve mass! 
xint = last_event; % updates initial conditions with final conditions of previous iteration
tinit = te; % updates initial time with final time of previous iteration

% Appends solutions and times
x = [x;tempx]; %appends the solutions I WONDER HOW ThIS WORKS WITH MATRICES
t = [t ; tempt]; %appends the times and then removes first timestep  (apart from first iteration).
if i > 1
    x((1),:) = []; %replaces last row with an empty array
    t(1)=[];
end
display('successfully completed one iteration!!')
end


%solutions
m = x(:,1);
CycBt = x(:,2);
Cdh1a= x(:,3);
Cdc20t = x(:,4);
Cdc20a = x(:,5);
IEP = x(:,6);
CKIt = x(:,7);
SK = x(:,8);
TF = x(:,9);

%PLOTS
figure(1)

plot(t,m,'black',t,CycBt,'r',t,Cdh1a,'b')
xlabel('Time')
legend('Mass m','CycBt','Cdh1a')
