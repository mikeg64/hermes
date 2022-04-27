


%% Import the data
%%data = xlsread('VALMc_rho_132_test_sac_all.dat','VALMc_rho_2048_test');

%% Initialize variables.
filename = '/Users/mikegriffiths/proj/smaug_jet_master/matlab/createconfig/VALMc_rho_132_test_sac_all.dat';




consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7




%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy
mu=0.6d0;
R=8.31e3;




%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%16f%16f%16f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

nx2=132;
g=287;
rheight=dataArray{1,1};
densg=dataArray{1,3};
presgo=dataArray{1,4};
tempg=dataArray{1,2};


output=zeros(nx2,4);


lambda=tempg.*R./(consts.ggg.*mu);

tempg0=4382.7163;
rho0=2.1754899e-05;
p0=1320.5359;
index0=8;

dx=rheight(1)-rheight(2);
for i=1:nx2-1
    fx=1.0./lambda;
    fxc=fx(nx2:-1:i)
    integ=inte(fxc,dx)
    presg(i)=p0.*exp(-integ);   
end

for i=1:nx2-1
    fx=1.0./lambda;
    fxc=fx(nx2:-1:i)
    integ=inte(fxc,dx)
    ndensg(i)=(tempg0*rho0./tempg(i)).*exp(-integ);   
end

presg1=presg;

 
presg(132)=(presg(131)+presg(130))/2;



ndensg=ndensg';

ndensg(2)=(ndensg(3)+ndensg(4))/2;
ndensg(1)=(ndensg(3)+ndensg(2))/2;
ndensg(131)=(ndensg(130)+ndensg(129))/2;
ndensg(132)=(ndensg(130)+ndensg(131))/2;

output(:,1)=rheight;
output(:,2)=tempg;
output(:,3)=ndensg;
output(:,4)=presg;

writematrix(output,'test.dat');





