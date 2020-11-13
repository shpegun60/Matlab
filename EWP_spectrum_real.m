%%%%% EWP spectrum analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% by shpegun60%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; % clear command line and workspace
syms x n;
% definition symbols
assume(n,'integer');
assumeAlso(n >= 0);
assume(x >= 0);

%define variables for deep understand
notParityConstant = 2*(n+1);                      % not parity constant
mulConstant = sqrt(notParityConstant/(pi*(n+2))); % multiplication constant for parity and non parity parameters
cosNotParityArgument = cos(x*(2*n + 3));          % cos with non parity argument
sinParityArgument = sin(x*(2*n + 2));             % sin with parity argument
sinNotParityArgument = sin(x*(2*n + 3));          % sin with non parity argument
cosParityArgument = cos(x*(2*n + 4));             % cos with parity argument

% non parity EWP parameter definition ****************************************************
U2n1 = mulConstant * ((sinParityArgument/(notParityConstant*sin(x))) - cosNotParityArgument);
% parity EWP parameter definition ********************************************************
U2n2 = mulConstant * ((1/notParityConstant)*(sinNotParityArgument/sin(x) - 1) - cosParityArgument);


%analysis main function ***************************************************
deltaT =  0.000000001; % 1 ns from oscilloscope
V_Trigger = 100;       % trigger for find impulse (mV)
start_index = 0;
stop_index = 0;

%------------ my signal ------------------------------
SIGNAL = dlmread('signal2.txt');%<-----------------------
%-------------------------------------------------------

mySignalDX = diff(SIGNAL(:,2));
mySignalPowABS = abs(SIGNAL(:,2)/1000).^2;
ENERGY = 0;
% find impulse, integral, energy  ----------------
myFunctionIntegral = 0;
flag_start_find = 0;
for i=1:length(SIGNAL(:,2))
    if (flag_start_find == 0)
        if abs(SIGNAL(i,2)) > V_Trigger
            start_index = i-2;
            flag_start_find = 1;
        end
    elseif (flag_start_find == 1)
        if abs(SIGNAL(i,2)) < V_Trigger
            stop_index = i+2;
            flag_start_find = 2;
        end
    end
    myFunctionIntegral = myFunctionIntegral + SIGNAL(i,2) * deltaT;
    ENERGY = ENERGY + (mySignalPowABS(i) * deltaT); % for parseval eq
end



figure('Name','my signal','NumberTitle','off');
subplot(2,1,1);
plot(SIGNAL(:,1), SIGNAL(:,2));
title('my signal');
set(gca, 'XTick',0:1:length(SIGNAL(:,1)));
set(gca, 'YTick',-6000:1000:6000);
grid on;
xlabel('ns')
ylabel('mV')
line([start_index start_index],[-6000 6000],'Color','red','LineStyle','--');
line([stop_index stop_index],[-6000 6000],'Color','red','LineStyle','--');

subplot(2,1,2);
plot(mySignalDX);
title('my signal differential');
set(gca, 'XTick',0:1:length(SIGNAL(:,1)));
set(gca, 'YTick',-6000:1000:6000);
grid on;
xlabel('ns')
ylabel('f ''(x) mV')

% find elementary wawe coeficient A2n1 and B2n2 for some function***************************
integralSpaceN = 0:1:30;
tau = pi - deltaT;
timeScaleTransformationCoef = tau/(stop_index - start_index);
fprintf('integral -> %.10f \t\n time-scale coefficient -> %.4f \t\n',myFunctionIntegral, timeScaleTransformationCoef);

array_spectrum_non_parity = zeros(length(integralSpaceN),2,'double');
array_spectrum_parity = zeros(length(integralSpaceN),2,'double');

U2n1_N = 0;
U2n2_N = 0;
ParsevalCoef = 0;
myFunctionRecovery = 0;

for i=1:length(integralSpaceN)
    array_spectrum_non_parity(i,1) = integralSpaceN(i);
    array_spectrum_parity(i,1) = integralSpaceN(i);
    
    U2n1_N = subs(U2n1,n,integralSpaceN(i));
    U2n2_N = subs(U2n2,n,integralSpaceN(i));
    
    for j = start_index : stop_index
        array_spectrum_non_parity(i,2) =  array_spectrum_non_parity(i,2) + (SIGNAL(j,2) * subs(U2n1_N, x, ((j - start_index + deltaT) * timeScaleTransformationCoef)) * deltaT);
        array_spectrum_parity(i,2) =  array_spectrum_parity(i,2) + (SIGNAL(j,2) * subs(U2n2_N, x, ((j - start_index + deltaT) * timeScaleTransformationCoef)) * deltaT);
    end
   ParsevalCoef = ParsevalCoef + ((array_spectrum_non_parity(i,2)^2) + (array_spectrum_parity(i,2)^2));
   myFunctionRecovery = myFunctionRecovery + (array_spectrum_non_parity(i,2)* U2n1_N + array_spectrum_parity(i)* U2n2_N);
end

% plot coefficients------------------------------------------------------
figure('Name','Non-parity and parity EWP coeficients','NumberTitle','off');
% plot non parity EWP coeficient
subplot(2,1,1);
%plot(integralSpaceN,notParityCoeficients,'LineStyle','none','-s','MarkerIndices');
for i=1:length(integralSpaceN)
    %line([integralSpaceN(i) integralSpaceN(i)], [0 array_spectrum_non_parity(i,2)],'Marker','s','LineWidth', 2);
    line([integralSpaceN(i) integralSpaceN(i)], [0 array_spectrum_non_parity(i,2)],'LineWidth', 2);
end
title('Non-parity EWP coeficient');
set(gca, 'XTick',0:1:integralSpaceN(length(integralSpaceN)));
grid on;
xlabel('n')
ylabel('A2n+1')

% plot parity EXP coeficient
subplot(2,1,2);
%plot(integralSpaceN,parityCoeficients,'-s','MarkerIndices',1:1:length(parityCoeficients));
for i=1:length(integralSpaceN)
    line([integralSpaceN(i) integralSpaceN(i)], [0 array_spectrum_parity(i,2)],'LineWidth', 2);
end
set(gca, 'XTick',0:1:integralSpaceN(length(integralSpaceN)));
title('Parity EWP coeficient');
grid on;
xlabel('n')
ylabel('B2n+2')

%%%%%%% Schedule check Parseval equality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\t\n\t\n---- Parseval equality ------\t\n');
fprintf('integral(abs(function)^2) -> %.10f \t\nPower sum EWP coeficients -> %.10f\t\n',ENERGY, ParsevalCoef);



%%%%%%% Recovery my function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integralSpaceX =deltaT:(pi-deltaT)/(stop_index - start_index):pi-deltaT;
plotRecoverFunction = subs(myFunctionRecovery,x,integralSpaceX);

% rmsValue = 0;
% for i=1:length(integralSpaceX)
%     rmsValue = rmsValue + ((plotRecoverFunction(i) - SIGNAL(integralSpaceX(i),2))^2) / length(integralSpaceX);
% end
% rmsValue = sqrt(rmsValue);

%fprintf('rms different value -> %f \t\n',rmsValue);

figure('Name','Recovered Function','NumberTitle','off');
plot(integralSpaceX,plotRecoverFunction);
title('my function Recovered');
% set(gca, 'XTick',0:1:length(SIGNAL(:,1)));
% set(gca, 'YTick',-6000:1000:6000);
% set(gca, 'XTick',0:0.2:pi);
axis([0 pi -inf inf]);
grid on;
xlabel('x')
ylabel('f(x)')




