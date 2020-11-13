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

%plot parameters **********************************************************
% plotX = linspace(0.001, pi-0.001, 100); % space for x value
% plotN = linspace(0, 3, 4);         % space for n value     
% 
% U2n1X = subs(U2n1,x,plotX);
% U2n2X = subs(U2n2,x,plotX);
% for i=1:length(plotN)
%    U2n1Plot(i,:) = subs(U2n1X,n,plotN(i));
%    U2n2Plot(i,:) = subs(U2n2X,n,plotN(i));
% end
% 
% figure('Name','Non parity and parity EWP parameters','NumberTitle','off');
% % plot non parity EWP parameter
% subplot(2,1,1);
% plot(plotX,U2n1Plot);
% title('Non-parity parameter');
% set(gca, 'XTick',0:0.2:pi);
% axis([0 pi -1 1]);
% grid on;
% xlabel('x')
% ylabel('U2n+1')
% legend({'n = 0','n = 1','n = 2','n = 3'},'Location','southwest');
% legend('boxoff');
% 
% % plot parity EWP parameter
% subplot(2,1,2);
% plot(plotX,U2n2Plot);
% title('Parity parameter');
% set(gca, 'XTick',0:0.2:pi);
% axis([0 pi -1.2 1]);
% grid on;
% xlabel('x');
% ylabel('U2n+2');
% legend({'n = 0','n = 1','n = 2','n = 3'},'Location','southwest');
% legend('boxoff');


%analysis main function ***************************************************
tau = pi;
m = 5;
amplitude = 1;
frequency = 5;
tetaHeaviside = (heaviside(x/tau) - heaviside((x/tau)-1));

%myFunction = sqrt(2/pi)*sin(6*x);           %%
%myFunction =sqrt(2/pi)*((x*(pi-x))^0.1) * sin(6*x);

% finit models
%myFunction = tetaHeaviside*(sin((pi*x)/tau)^2)*sin(2*pi*frequency*x); %1 model kenno
%myFunction = ((-1)^m)*sin((2*pi*m*x)/pi)*tetaHeaviside;
%myFunction = ((-1)^m)*(1-abs(((2*x)/tau)-1))*sin((2*pi*m*x)/tau)*tetaHeaviside;
%myFunction = (1 - abs(((2*x)/tau)-1))*cos((4*pi*m*x)/tau)*tetaHeaviside;


% quazi finit models
%myFunction = exp(-1*x)*sin(2*pi*frequency*x)*heaviside(x);
%myFunction = amplitude*(1 - (x/tau))*exp((-x)/tau)*heaviside(x); %%%% has no effect
%myFunction = amplitude*(exp(-10*x) - exp(-1.5*x))*cos(sqrt(10*1.5)*x)*heaviside(x);
%myFunction = amplitude*(exp(-1.5*x) - exp(-(1.5+3)*x))*cos(2*pi*frequency*x); % has no name
%myFunction = ((-1)^m)*sin((2*pi*m*x)/tau)*exp(-1*((((2*x)/tau) - 1)^2));
%myFunction = amplitude*cos((2*pi*m*x)/tau)*exp(-1*((((2*x)/tau) - 1)^2));


%myFunction = (exp(-1.5*x) - exp(-(1.5+3)*x))*cos(2*pi*sin(sqrt(4*pi)*x));% elite signal 

% fractal
beta = 0.51;
%myFunction = sin(2*pi*x)*tetaHeaviside + 0.25*sin(2*pi*4*x)*tetaHeaviside + (3^(-2))*sin(2*pi*9*x)*tetaHeaviside + (4^(-2))*sin(2*pi*16*x)*tetaHeaviside + (5^(-2))*sin(2*pi*25*x)*tetaHeaviside;
myFunction = (2/(pi^beta))*(   (1^(-2*beta))*sin(2*pi*(1^2)*x)*tetaHeaviside  + (2^(-2*beta))*sin(2*pi*(2^2)*x)*tetaHeaviside + (3^(-2*beta))*sin(2*pi*(3^2)*x)*tetaHeaviside + (4^(-2*beta))*sin(2*pi*(4^2)*x)*tetaHeaviside + (5^(-2*beta))*sin(2*pi*(5^2)*x)*tetaHeaviside);

myFunctionDX = diff(myFunction,x);
myFunctionIntegral = amplitude*int(myFunction ,x,0.000001,pi-0.000001);
disp('--- conditions check---');
fprintf('integral -> %f \t\n f(0) -> %f \t\n f(pi) -> %f \t\n dx(0) -> %f \t\n dx(pi) -> %f',myFunctionIntegral,subs(myFunction,x,0.0000001),subs(myFunction,x,pi-0.0000001),subs(myFunctionDX,x,0.0000001),subs(myFunctionDX,x,pi-0.0000001));

spaceMFunction = 2:2:6;
spaceXFunction = linspace(0.001, pi-0.001, 600);

functionX = subs(myFunction,x,spaceXFunction);
functionDX = subs(myFunctionDX,x,spaceXFunction);

figure('Name','my function','NumberTitle','off');
subplot(2,1,1);
% plot my function
plot(spaceXFunction,functionX);
title('my function');
set(gca, 'XTick',0:0.2:pi);
axis([0 pi -inf inf]);
grid on;
xlabel('x')
ylabel('f(x)')

subplot(2,1,2);
% plot my function
plot(spaceXFunction,functionDX);
title('my function differential');
set(gca, 'XTick',0:0.2:pi);
axis([0 pi -inf inf]);
grid on;
xlabel('x')
ylabel('f ''(x)')

%*********************************************
% find elementary wawe coeficient A2n1 and B2n2 for some function***************************
integralNotParityFunction = myFunction * U2n1;
integralParityFunction = myFunction * U2n2;

deltaT = 0.0104719088;%0.0313959265; %%% 0.062791853 %% 0.0104719088;
integralSpaceN = 0:1:50;
integralSpaceX =0.00001:deltaT:pi-0.00001;

parityCoeficients = zeros(length(integralSpaceN),1,'double');
notParityCoeficients = zeros(length(integralSpaceN),1,'double');

integralNotParityFunctionX = subs(integralNotParityFunction,x,integralSpaceX);
integralParityFunctionX = subs(integralParityFunction,x,integralSpaceX);

ParsevalCoef = zeros(1,1,'double');
myFunctionRecovery = 0;
% find A2n1 and B2n2 in cycle and put to 2-dimensional array
for i=1:length(integralSpaceN)
    for j=1:length(integralSpaceX)
        parityCoeficients(i) = parityCoeficients(i) + (subs(integralParityFunctionX(j),n,integralSpaceN(i))*deltaT);
        notParityCoeficients(i) = notParityCoeficients(i) + (subs(integralNotParityFunctionX(j),n,integralSpaceN(i))*deltaT);
    end
    ParsevalCoef = ParsevalCoef + ((parityCoeficients(i)^2) + (notParityCoeficients(i)^2));
    myFunctionRecovery = myFunctionRecovery + (notParityCoeficients(i)* subs(U2n1,n,i-1) + parityCoeficients(i)* subs(U2n2,n,i-1));
end

figure('Name','Non-parity and parity EWP coeficients','NumberTitle','off');
% plot non parity EWP coeficient
subplot(2,1,1);
plot(integralSpaceN,notParityCoeficients,'-s','MarkerIndices',1:1:length(notParityCoeficients));
title('Non-parity EWP coeficient');
grid on;
xlabel('n')
ylabel('A2n+1')

% plot parity EXP coeficient
subplot(2,1,2);
plot(integralSpaceN,parityCoeficients,'-s','MarkerIndices',1:1:length(parityCoeficients));
title('Parity EWP coeficient');
grid on;
xlabel('n')
ylabel('B2n+2')

%%%%%%% Schedule check Parseval equality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CheckSchedule = zeros(1,1,'double');
myFunctionABS = abs(myFunction)^2;
ParsevalIntegral = subs(myFunctionABS,x,integralSpaceX);
fprintf('\t\n\t\n---- Parseval equality ------\t\n');
for k=1:length(integralSpaceX)
    CheckSchedule = CheckSchedule + (ParsevalIntegral(k)*deltaT);
end
fprintf('integral(abs(function)^2) -> %f \t\nPower sum EWP coeficients -> %f\t\n',CheckSchedule,ParsevalCoef);

%%%%%%% Recovery my function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotRecoverFunction = subs(myFunctionRecovery,x,spaceXFunction);

rmsValue = 0;
for i=1:length(spaceXFunction)
    rmsValue = rmsValue + ((plotRecoverFunction(i) - functionX(i))^2)/length(spaceXFunction);
end
rmsValue = sqrt(rmsValue);

fprintf('rms different value -> %f \t\n',rmsValue);

figure('Name','Recovered Function','NumberTitle','off');
plot(spaceXFunction,plotRecoverFunction);
title('my function Recovered');
set(gca, 'XTick',0:0.2:pi);
axis([0 pi -inf inf]);
grid on;
xlabel('x')
ylabel('f(x)')



