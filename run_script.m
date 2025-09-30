%program which does empirical illustration for chapter 3
%load in the data set. Here use house price data from hprice.txt
load data\house_prices.txt;

n=size(house_prices,1);
y=house_prices(:,1);
x=house_prices(:,2:5);
x=[ones(n,1) x];
k=5;

%% Model run with weakly informative priors

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=2.4*eye(k);
capv0(2,2)=6e-7;
capv0(3,3)=.15;
capv0(4,4)=.6;
capv0(5,5)=.6;
capv0inv=inv(capv0);

%Call script which carries actually does posterior analysis
calc_script;

%save the log of marginal likelihood for later use
lmargun=lmarglik;

%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02
capsd0 = diag(sqrt(capv0* (v0*s02)/(v0-2)));

'Posterior results based on Informative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
lmarglik
ystarm
ystarsd
ystarcapv

%% Model run with uniformative priors

%Hyperparameters for noninformative prior
v0=0;
capv0inv=0*eye(k);


%Call script which carries actually does posterior analysis
calc_script;

%Print out whatever you want
'Posterior results based on Noninformative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
ystarm
ystarsd
ystarcapv

%% Experiment with somewhat informative priors

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=0.864*eye(k);
capv0(2,2)=2.16e-07;
capv0(3,3)=.054;
capv0(4,4)=.216;
capv0(5,5)=.216;
capv0inv=inv(capv0);

% Check if prior variance is as intended
capsd0 = diag(sqrt(capv0* (v0*s02)/(v0-2)));


%Call script which carries actually does posterior analysis
calc_script;

%save the log of marginal likelihood for later use
lmargun=lmarglik;

%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02

'Posterior results based on Informative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
lmarglik
ystarm
ystarsd
ystarcapv

%% Experiment with more informative priors

%Hyperparameters for natural conjugate prior
v0=5;
b0=0*ones(k,1);
b0(2,1)=10;
b0(3,1)=5000;
b0(4,1)=10000;
b0(5,1)=10000;
s02=1/4.0e-8;
capv0=0.096*eye(k);
capv0(2,2)=2.4e-8;
capv0(3,3)=.006;
capv0(4,4)=.024;
capv0(5,5)=.024;
capv0inv=inv(capv0);

% Check if prior variance is as intended
capsd0 = diag(sqrt(capv0* (v0*s02)/(v0-2)));

%Call script which carries actually does posterior analysis
calc_script;

%save the log of marginal likelihood for later use
lmargun=lmarglik;

%Print out whatever you want
'Hyperparameters for informative natural conjugate prior'
b0
capv0
v0
s02
capsd0 = diag(sqrt(capv0* (v0*s02)/(v0-2)));

'Posterior results based on Informative Prior'
b1
bsd
probpos
bhpdi95
bhpdi99
hmean
hsd
lmarglik
ystarm
ystarsd
ystarcapv

%% Do I need analysis for non-conjugate priors and heteroscedasticity??