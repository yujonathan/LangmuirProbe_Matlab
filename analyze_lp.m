function [pAns,pFit,rFit]=analyze_lp(savedAns,savedFit,rFit)
% analyze_lp Summary     2015_07_07
%       Perform simple linear fits to langmuir probe
%       data. Several prompts need to be filled out
%       to attain the correct plasma parameters.
%
%       Automated fitting assumes an input potential
%       close to -100 to 50 V.
%
%       The last figure contains all three plotted
%       fits as well as the plasma parameters.
%
%       promptAns, promptFit, and repF can be fed back 
%       into savedAns, savedFit, and repF to utilize 
%       values input into analyze_lp from the last run.
%
%       Example for Command Window: 
%       (First fitting with default parameters)
%           >> [pAns,pFit,rFit] = analyze_lp();
%       (Repeat fit and adjust parameters)
%           >> [pAns,pFit,rFit] = analyze_lp(pAns,pFit,rFit);
%       (Repeat until satisfied with fit)
%
%       (Another data set fitting but using default fit parameters)
%           >> [pAns,pFit,rFit] = analyze_lp(pAns);
%       (Only the initial prompt is filled out with the last input)
%
% Commonly used gasses
%       E  = ( Z,   amu)
%       H  = ( 1, 1.008)
%       D  = ( 1, 2.014)
%       He = ( 2, 4.003)
%       Ar = (18,39.948)
%
% Commonly used default
%
% rfpower  =    '1300';       % [W]
% magnet   =      '50';       % [A]
% flow     =      '12';       % [sccm]
% pressure =     '4.3';       % [mT]
% bias     =  '-100.0';       % [V]
% prompt_0 =       '0';       % 0 = false (auto fill), 1 = true (prompt)  
% Res      =    '50.0';       % [Ohms]    Resistor current is measured across
% l_t      = '2.19e-3';       % [m]       length of tip
% d_t      = '6.35e-4';       % [m]       diameter of tip

%% ------------------------------------------------------
Element = 'He';

% default values for prompt
rfpower  =    '1300';       % [W]
magnet   =      '50';       % [A]
flow     =      '10';       % [sccm]
pressure =     '6.3';       % [mT]
bias     =  '-100.0';       % [V]
prompt_0 =       '1';       % 0 = false (auto fill), 1 = true (prompt)  
Res      =    '50.0';       % [Ohms]    Resistor current is measured across
l_t      =  '3.2e-3';       %'2.19e-3';  % [m] length of tip
d_t      =  '4.6e-4';       %'6.35e-4';  % [m] diameter of tip

if     strcmp(Element,'H')
        Z   = '1';
        amu = '1.008';
elseif strcmp(Element,'D')
        Z   = '1';
        amu = '2.014';
elseif strcmp(Element,'He')
        Z   = '2';
        amu = '4.003';
elseif strcmp(Element,'Ar')
        Z   = '18';
        amu = '39.948';
end

defaultAns = {Z,amu,rfpower,magnet,flow,pressure,bias,prompt_0,Res,l_t,d_t};

if exist('savedAns','var')
     pAns = savedAns;
else pAns = defaultAns;
end

%% ------------------------------------------------------
% define fitting defaults as percentages
nCut   = 3;             % Always cut off first and last n points
  
minI_s = 0.05;          % Start 05% between min and floating pot
maxI_s = 0.25;          % End   25%   " "
minN_i = 0.05;          % Start 05%
maxN_i = 0.25;          % End   25%
minT_e = 0.05;          % Start 05%
maxT_e = 0.17;          % End   17%
minE_s = 0.65;          % Start 65% 
maxE_s = 0.80;          % End   80%

%% ------------------------------------------------------
% define physical constants needed
q_e = 1.60217657e-19;                 % electron charge [C]
% m_e = 9.10938291e-31;                 % electron mass   [kg]
m_n = 1.660538921e-27;                % nucleon mass    [kg]

%% ------------------------------------------------------
% plotting variables
pFit.mark1 =  5.0;
pFit.mark2 =  8.5;
pFit.line1 =  1.5;

backColor = [0.95 0.95 0.95];
legColor  = [0.85 0.85 0.85];

%% ------------------------------------------------------
% get the experimental parameters from the user
prompt = {'species [Z]:',...
          'atomic mass [u]:',...
          'rf power [W]:',...
          'main magnet [A]:',...
          'gas flow rate [sccm]:',...
          'chamber pressure [mT]:',...
          'bias voltage [V]:',...
          'prompt = 1, auto = 0:'};
name = 'Experiment Info';
numlines = 1;
pAns(1:8) = inputdlg(prompt,name,numlines,pAns(1:8));

prompt = {'resistor (Ohms):'...
          'probe tip length (m):'...
          'probe tip diameter (m):'};
name = 'Probe Info';
numlines = 1;
pAns(9:11) = inputdlg(prompt,name,numlines,pAns(9:11));

Z        = str2double(pAns(1));
amu      = str2double(pAns(2));
rfpower  = str2double(pAns(3));
magnet   = str2double(pAns(4));
flow     = str2double(pAns(5));
pressure = str2double(pAns(6));
bias     = str2double(pAns(7));
prompt_0 = str2double(pAns(8));
Res      = str2double(pAns(9));
l_t      = str2double(pAns(10));
d_t      = str2double(pAns(11));

%% ------------------------------------------------------
% get the files for analysis and put them in arrays

if exist('rFit','var')
     file_vin     = rFit.file_vin;
     file_vout    = rFit.file_vout;
else [FileNames,PathName] = uigetfile('*.txt',...
                                      'Select 2 files',...
                                      'MultiSelect','on',...
                                      'C:\ProbeData_mjs\');
     file_vin  = strcat(PathName, FileNames{1});
     file_vout = strcat(PathName, FileNames{2});
     rFit.file_vin  = file_vin;
     rFit.file_vout = file_vout;
end

fid = fopen(file_vin);
V_p = fscanf(fid,'%g');
fclose(fid);

fid = fopen(file_vout);
I_p = fscanf(fid,'%g');         % Still in Volts, must convert to current
fclose(fid);

parsed_vin  = parseName(file_vin );
parsed_vout = parseName(file_vout);

%% ------------------------------------------------------
% setup for the figures
figureName  = [  ' V_p = ' parsed_vin ...
               ' | I_p = ' parsed_vout];
hFig = figure(1);
clf
set(hFig,'Name',figureName,...
    'ToolBar',      'none',...
    'units',  'normalized',...
    'outerposition', [0 0.05 0.95 0.95]);

%% ------------------------------------------------------
% Clean up data for plotting 
% cut off ends of the data
V_p =  V_p(nCut:end-nCut);
I_p =  I_p(nCut:end-nCut);

[~,iSort] = sort(V_p);
V_p = V_p(iSort);
I_p = I_p(iSort);
I_p = -I_p/Res;                     % Current [A]

winSize  = round(length(I_p)*0.01);	% determine size of window for running average
winSize2 = round(winSize/2);

I_a = tsmovavg(I_p,'t',winSize,1);
I_a = circshift(I_a,-winSize2);

V_p = V_p(winSize2+1:end-winSize2-1);
I_p = I_p(winSize2+1:end-winSize2-1);
I_a = I_a(winSize2+1:end-winSize2-1);

% I_a = smooth(V_p,I_a,0.01,'rlowess'); % slows too much, skews 3rd fit
%% ---------------------------------------------------------
% plot the full I-V curve
% get min and max Vp to fit to
if exist('savedFit','var')
    pFit = savedFit;       
else
    pFit.v0a = min(V_p);
    pFit.v0b = max(V_p);
end
pFit.i0a = iFind(V_p,pFit.v0a);
pFit.i0b = iFind(V_p,pFit.v0b);

figure(1)
clf
grid on
set(gca,'XMinorTick','on')
plotIV_p(V_p,I_p,pFit);

%min(V_p)
%max(V_p)

choice = prompt_0;
while choice == 1
 vals = getRange(V_p,pFit.i0a,pFit.i0b,...
                 {'Enter the min V_p :',...
                  'Enter the max V_p :'});
 pFit.v0a = vals(1);
 pFit.v0b = vals(2);
 pFit.i0a = iFind(V_p,pFit.v0a);
 pFit.i0b = iFind(V_p,pFit.v0b);
 
 % replot the adjusted I-V curve
 figure(1)
 clf
 grid on
 set(gca,'XMinorTick','on')
 plotIV_p(V_p,I_p,pFit);

 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 1
               'continue');          % results in choice = 2
end

V_p = V_p(pFit.i0a:pFit.i0b);        % save the newly cut vectors
I_p = I_p(pFit.i0a:pFit.i0b);
I_a = I_a(pFit.i0a:pFit.i0b);

%% ------------------------------------------------------
% limits for plotting
Vmin = min(V_p);
Vmax = max(V_p);
Imin = min(I_p);
Imax = max(I_p);
Vmm  = [Vmin,Vmax];
Imm  = [Imin,Imax];
pFit.Vmin = Vmin;
pFit.Vmax = Vmax;
pFit.Imin = Imin;
pFit.Imax = Imax;

%% ------------------------------------------------------
% fit points utilizing floating potential and above percentages
% find V_float and i_float
x0  = intersections(Vmm,[0,0],V_p,I_a);
V_f = mean(x0);
i_f = iFind(V_p,V_f);

pFit.V_f = V_f;
pFit.i_f = i_f;

% Initial guesses for fitting of data are either saved or default
if exist('savedFit','var')
else
    pFit.v1a = Vmin + (V_f-Vmin)*minI_s;
    pFit.v1b = Vmin + (V_f-Vmin)*maxI_s;
    pFit.v2a = Vmin + (V_f-Vmin)*minN_i;
    pFit.v2b = Vmin + (V_f-Vmin)*maxN_i;
    pFit.v3a = V_f  + (Vmax-V_f)*minT_e;
    pFit.v3b = V_f  + (Vmax-V_f)*maxT_e;
    pFit.v4a = V_f  + (Vmax-V_f)*minE_s;
    pFit.v4b = V_f  + (Vmax-V_f)*maxE_s;
end

%% ------------------------------------------------------
% Initial guess for n_ion fitting
% pFit.i2a = iFind(V_p,-90.0);
% pFit.i2b = iFind(V_p,-75.0);

pFit.i2a = iFind(V_p,pFit.v2a);
pFit.i2b = iFind(V_p,pFit.v2b);

figure(1)
clf
grid on
set(gca,'XMinorTick','on')
dI2 = plotFitN_i(V_p,I_p,I_a,pFit);

choice = prompt_0;
while choice == 1 
 % get n_ion fitting parameters from user       
 vals = getRange(V_p,pFit.i2a,pFit.i2b,...
                 {'Enter the min V_p for n_ion fit:',...
                  'Enter the max V_p for n_ion fit:'});
 pFit.v2a = vals(1);
 pFit.v2b = vals(2);
 pFit.i2a = iFind(V_p,pFit.v2a);
 pFit.i2b = iFind(V_p,pFit.v2b);

 figure(1)
 clf
 grid on
 set(gca,'XMinorTick','on')
 dI2 = plotFitN_i(V_p,I_p,I_a,pFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 1
               'continue');          % results in choice = 2
end

%% ------------------------------------------------------
% Initial guess for I_sat fitting
% pFit.i1a = iFind(V_p,-90.0);
% pFit.i1b = iFind(V_p,-75.0);

pFit.i1a = iFind(V_p,pFit.v1a);
pFit.i1b = iFind(V_p,pFit.v1b);

figure(1)
clf
grid on
set(gca,'XMinorTick','on')
I_sat = plotFitI_s(V_p,I_p,I_a,pFit);


choice = prompt_0;
while choice == 1 
 % get I_sat fitting parameters from user       
 vals = getRange(V_p,pFit.i1a,pFit.i1b,...
                 {'Enter the min V_p for I_sat fit:',...
                  'Enter the max V_p for I_sat fit:'});
 pFit.v1a = vals(1);
 pFit.v1b = vals(2);
 pFit.i1a = iFind(V_p,pFit.v1a);
 pFit.i1b = iFind(V_p,pFit.v1b);
 
 figure(1)
 clf
 grid on
 set(gca,'XMinorTick','on')
 I_sat = plotFitI_s(V_p,I_p,I_a,pFit);

 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 1
               'continue');          % results in choice = 2
end

%% ------------------------------------------------------
% begin fit for T_e and E_sat
LI_p = -abs(log(I_p-I_sat));
LI_a = -abs(log(I_a-I_sat));

% initial guesses for Te and Esat fit
% pFit.i3a = iFind(V_p,11.0);
% pFit.i3b = iFind(V_p,14.0);
% pFit.i4a = iFind(V_p,28.0);
% pFit.i4b = iFind(V_p,35.0);

pFit.i3a = iFind(V_p,pFit.v3a);
pFit.i3b = iFind(V_p,pFit.v3b);
pFit.i4a = iFind(V_p,pFit.v4a);
pFit.i4b = iFind(V_p,pFit.v4b);

figure(1)
clf
grid on
set(gca,'XMinorTick','on')
[V_s,T_e] = plotFitTE(V_p,LI_p,LI_a,pFit);

choice = prompt_0;
while choice == 1 
 % get T_e and E_sat fitting parameters from user       
 vals = getRange2(V_p,pFit.i3a,pFit.i3b,...
                      pFit.i4a,pFit.i4b,...
                  {'Enter the min V_p for T_e fit:',...
                   'Enter the max V_p for T_e fit:'...
                   'Enter the min V_p for E_sat fit',...
                   'Enter the max V_p for E_sat fit'});
 pFit.v3a = vals(1);
 pFit.v3b = vals(2);
 pFit.i3a = iFind(V_p,pFit.v3a);
 pFit.i3b = iFind(V_p,pFit.v3b);
 pFit.v4a = vals(3);
 pFit.v4b = vals(4);
 pFit.i4a = iFind(V_p,pFit.v4a);
 pFit.i4b = iFind(V_p,pFit.v4b);
 
 figure(1)
 clf
 set(gca,'XMinorTick','on')
 grid on
 [V_s,T_e] = plotFitTE(V_p,LI_p,LI_a,pFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 1
               'continue');          % results in choice = 2
end

%% ------------------------------------------------------
% calculations and final plotting

% probe area
A_p = pi*(d_t/2)^2 + pi*d_t*l_t;  % [m^2]

% Ion flux and n_e
phi = I_sat/-q_e/A_p;             % [ions/m^2/s]
m_i = amu*m_n;                    % [kg] (Atomic mass)
n_e = phi/sqrt(q_e*T_e/m_i);      % [m^-3] (plasma density n_0) 
% exp(-Z/2)/sqrt(e*Z*T_e/m_i);    % [m^-3]

% Ion energy
E_ion = V_s-bias;                 % [eV]

% Ion density
n_ion = sqrt(pi^2*m_i/2/A_p^2/q_e^3*abs(dI2));  % from OML theory, m-3

% defining info for legends
IsatInfo = ['I_{sat} fit:'...
            '   V_{min} = ' num2str(pFit.v1a,3)...
            '   V_{max} = ' num2str(pFit.v1b,3)];
I2pInfo  = ['n_{ion} fit:'...
            '   V_{min} = ' num2str(pFit.v2a,3)...
            '   V_{max} = ' num2str(pFit.v2b,3)];                                          
TeInfo   = ['T_{e  } fit:'...
            '   V_{min} = ' num2str(pFit.v3a,3)...
            '   V_{max} = ' num2str(pFit.v3b,3)];                 
EsatInfo = ['E_{sat} fit:'...
            '   V_{min} = ' num2str(pFit.v4a,3)...
            '   V_{max} = ' num2str(pFit.v4b,3)];
TeEsatInfo = {TeInfo,EsatInfo};

%% ------------------------------------------------------
%------------------- Plotting -------------------------
%---------------- Final I-V curve ---------------------
%------------------------------------------------------

figure(1)
clf

% Final I_sat curve
subplot(2,2,1)
set(subplot(2,2,1),'Color',backColor)
title(IsatInfo)
plotI_s(V_p,I_p,I_a,pFit);
h=legend('Raw Data',...
         'Run Avg',...
         'I_{sat} fit');
set(h,'Location','northwest','color',legColor)

% Final n_ion curve
subplot(2,2,2)
set(subplot(2,2,2),'Color',backColor)
title(I2pInfo)
plotFitN_i(V_p,I_p,I_a,pFit);
h=legend('Raw Data',...
         'Run Avg',...
         'n_{ion} fit');
set(h,'Location','southwest','color',legColor)

% Final ln(Ip-I_sat) vs Vp curve
subplot(2,2,3)
set(subplot(2,2,3),'Color',backColor)
title(TeEsatInfo)
plotFitTE(V_p,LI_p,LI_a,pFit);
h = legend('Raw Data',...
           'Run Avg',...
           'T_e fit',...
           'E_{sat} fit',...
           'V_{plasma}');
set(h,'Location','southeast','color',legColor)

% ensure the last legend is displayed
h4 = subplot(2,2,4);
set(h4,'visible','off');

% summary of measurements
ah=gca;
axes('position',[0,0,1,1],'visible','off');

phiInfo = '\phi_{ion} [m^{-2} s^{-1}] = ';

hrz1 = 0.68;
hrz2 = 0.86;
hght = 0.47;
delh = 0.045; 

text(hrz1,hght-delh*0,     'species [Z] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*1, 'atomic mass [u] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*2,    'rf power [W] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*3, 'main magnet [A] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*4, 'gas flow [sccm] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*5,   'pressure [mT] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*6,        'bias [V] = ','HorizontalAlignment','right');

text(hrz2,hght-delh*0,   'V_{float} [V] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*1,  'V_{plasma} [V] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*2,        'T_e [eV] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*3,    'n_0 [m^{-3}] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*4,'n_{ion} [m^{-3}] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*5,    'E_{ion} [eV] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*6,             phiInfo ,'HorizontalAlignment','right');

text(hrz1,hght-delh*0, num2str(       Z));
text(hrz1,hght-delh*1, num2str(     amu));
text(hrz1,hght-delh*2, num2str( rfpower));
text(hrz1,hght-delh*3, num2str(  magnet));
text(hrz1,hght-delh*4, num2str(    flow));
text(hrz1,hght-delh*5, num2str(pressure));
text(hrz1,hght-delh*6, num2str(    bias));

text(hrz2,hght-delh*0, num2str(   V_f,4));
text(hrz2,hght-delh*1, num2str(   V_s,4));
text(hrz2,hght-delh*2, num2str(   T_e,4));
text(hrz2,hght-delh*3, num2str(   n_e,4));
text(hrz2,hght-delh*4, num2str( n_ion,4));
text(hrz2,hght-delh*5, num2str( E_ion,4));
text(hrz2,hght-delh*6, num2str(   phi,4));

% location of files used 
text(.60,.120,  'V_p path: ','HorizontalAlignment','right');
text(.60,.075, ' I_p path: ','HorizontalAlignment','right');
text(.60,.125,  parsed_vin ,'interpreter','none');
text(.60,.080,  parsed_vout,'interpreter','none');

axes(ah);

end				% finish function

%% -------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------
%----------- Subroutines for this function -------------
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------

%% get the range for to perform fit
function [Xvalue] = getRange( VP,iMin,iMax,pStr)
 prompt = pStr;
 name = 'SELECT DATA RANGE TO FIT';
 numlines = 1;
 minVi = num2str(VP(iMin),3);
 maxVi = num2str(VP(iMax),3);
 defAnswer = {minVi,maxVi};
 Xvalue = inputdlg(prompt,name,numlines,defAnswer);
 Xvalue = str2double(Xvalue);
end

function [Xvalue] = getRange2(VP,iMin,iMax,jMin,jMax,pStr)
 prompt = pStr;
 name = 'SELECT DATA RANGE TO FIT';
 numlines = 1;
 minVi = num2str(VP(iMin),3);
 maxVi = num2str(VP(iMax),3);
 minVj = num2str(VP(jMin),3);
 maxVj = num2str(VP(jMax),3);
 defAnswer = {minVi,maxVi,minVj,maxVj};
 Xvalue = inputdlg(prompt,name,numlines,defAnswer);
 Xvalue = str2double(Xvalue);
end

%% plot I-V curve
function [] = plotIV_p(VP,IP,fitV)
 iMin = fitV.i0a;
 iMax = fitV.i0b;
 figure(1)
 clf
 title('I-V Curve')
 xlabel('V_p [V]')
 ylabel('I_p [mA]')
 axis([min(VP(iMin:iMax))     max(VP(iMin:iMax))...
       min(IP(iMin:iMax))*1e3 max(IP(iMin:iMax))*1e3])
 grid on
 hold on
 plot( VP(iMin:iMax),IP(iMin:iMax)*1e3,'db','MarkerSize', fitV.mark1)
end

%% plot and fit I_sat curve
function [I_s,I_std] = plotFitI_s(VP,IP,IA,fitV)
iMin = fitV.i1a;
iMax = fitV.i1b;

xMin = fitV.Vmin;
xMax = fitV.V_f;
yMin = fitV.Imin*1e3;
yMax = 0;

p   = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
IF  =  VP*p(1)+p(2);
IPS = -VP*p(1)+IP;
IAS = -VP*p(1)+IA;
I_s =  p(2);
Iss = [I_s,I_s];
VV  = [fitV.Vmin,fitV.V_f];
V1 = [VP(iMin),VP(iMax)];
I1 = [IF(iMin),IF(iMax)];

diffV = IP(iMin:iMax) - IF(iMin:iMax);
I_std = std(diffV);

subplot(2,1,1);
set(gca,'XMinorTick','on')
xlabel('V_p [V]')
ylabel('I_p [mA]')
axis([xMin xMax yMin yMax])
grid on
hold on
plot( VP, IP*1e3,  'xb',...
      VP, IA*1e3,  '.c',...
      VP, IF*1e3,  '-r',...
    'MarkerSize', fitV.mark1,...
     'LineWidth', fitV.line1)     
plot( V1, I1*1e3,  'or',...         % Show fitting limits
     'MarkerFaceColor', 'r',...
     'MarkerSize',fitV.mark2)
 
yMin = (I_s - (I_s-min(IA))*0.20)*1e3;
yMax = (I_s + (I_s-min(IA))*0.20)*1e3;

subplot(2,1,2);
set(gca,'XMinorTick','on')
xlabel('V_p [V]')
ylabel('I_p - I_{slope of sat} [mA]')
axis([xMin xMax yMin yMax])
grid on
hold on
plot( VP, IPS*1e3, 'xb',...
      VP, IAS*1e3, '.c',...
      VV, Iss*1e3,'--r',...
     'MarkerSize', fitV.mark1,...
      'LineWidth', fitV.line1)
end

%% plot only the I_sat curve
function [I_s,I_std] = plotI_s(VP,IP,IA,fitV)
iMin = fitV.i1a;
iMax = fitV.i1b;

xMin = fitV.Vmin;
xMax = fitV.V_f;
yMin = fitV.Imin*1e3;
yMax = 0;

p   = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
IF  =  VP*p(1)+p(2);
I_s =  p(2);
V1 = [VP(iMin),VP(iMax)];
I1 = [IF(iMin),IF(iMax)];

diffV = IP(iMin:iMax) - IF(iMin:iMax);
I_std = std(diffV);

xlabel('V_p [V]')
ylabel('I_p [mA]')
axis([xMin xMax yMin yMax])
hold on
plot( VP, IP*1e3,  'xb',...
      VP, IA*1e3,  '.c',...
      VP, IF*1e3,  '-r',...
    'MarkerSize', fitV.mark1,...
     'LineWidth', fitV.line1)     
plot( V1, I1*1e3,  'or',...         % Show fitting limits
     'MarkerFaceColor', 'r',...
     'MarkerSize',    fitV.mark2)
end

%% plot and fit n_ion curve
function [derI2,n_std] = plotFitN_i(VP,IP,IA,fitV)
iMin = fitV.i2a;
iMax = fitV.i2b;

xMin = fitV.Vmin;
xMax = fitV.V_f;
yMin = 0;
yMax = mean(IP(1:200))^2*1e6;

p   = polyfit(VP(iMin:iMax),IA(iMin:iMax).^2,1);
nF  = VP*p(1)+p(2);
derI2 =   p(1);
diffV = IP(iMin:iMax).^2 - nF(iMin:iMax).^2;
n_std = std(diffV);

V1 = [VP(iMin),VP(iMax)];
I1 = [nF(iMin),nF(iMax)];

xlabel('V_p [V]')
ylabel('I^2_p [mA]^2')
axis([xMin xMax yMin yMax])
hold on
plot( VP, (IP.^2)*1e6,  'xb',...
      VP, (IA.^2)*1e6,  '.c',...
      VP,      nF*1e6,  '-r',...
         'MarkerSize', fitV.mark1,...
          'LineWidth', fitV.line1);
plot( V1,      I1*1e6,  'or',...
    'MarkerFaceColor',   'r',...
         'MarkerSize', fitV.mark2);
end

%% plot and fit T_e and E_sat curve
function [V_s,T_e,T_std] = plotFitTE(VP,IP,IA,fitV)
iMin = fitV.i3a;
iMax = fitV.i3b;
jMin = fitV.i4a;
jMax = fitV.i4b;

xMin = fitV.V_f;
xMax = fitV.Vmax;
yMin = IA(fitV.i_f);
yMax = IA(end);

p1 = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
TF = VP*p1(1)+ p1(2);
T_e   = 1/p1(1);
diffV = IP(iMin:iMax) - TF(iMin:iMax);
T_std = std(diffV);

Vi = [VP(iMin),VP(iMax)];
Ii = [TF(iMin),TF(iMax)];

p2 = polyfit(VP(jMin:jMax),IA(jMin:jMax),1);
EF = VP*p2(1) + p2(2);
diffV = IP(jMin:jMax) - EF(jMin:jMax);
E_std = std(diffV);

Vj = [VP(jMin),VP(jMax)];
Ij = [EF(jMin),EF(jMax)];

% find V_plasma at intersection of T_e and n_ion
V_s = -(p1(2)-p2(2))/(p1(1)-p2(1));
Vss = [V_s,V_s];
Imm = [min(IA),max(IA)];

xlabel('V_p [V]')
ylabel('ln(I_p-I_{sat}) (A)')
axis([xMin xMax yMin yMax])
hold on
plot(     VP,  IP,  'xb',...
          VP,  IA,  '.c',...
          VP,  TF,  '-r',...               % fit with the initial guess
          VP,  EF,  '-g',...
         Vss, Imm, '--k',...
     'MarkerSize', fitV.mark1,...
      'LineWidth', fitV.line1)
plot(     Vi,  Ii,       'or',...
     'MarkerFaceColor',   'r',...
          'MarkerSize', fitV.mark2)
plot(     Vj,  Ij,       'og',...
     'MarkerFaceColor',   'g',...
          'MarkerSize', fitV.mark2)      
end

%% simple loop to find i corresponding to x_i
function [iVal] = iFind(Vi,xVal)
 iVal = 1;
 for i=1:size(Vi,1)
     if Vi(i)>xVal
        iVal = i;
        break
     else 
	iVal = length(Vi)-1;
     end
 end
end

%% shorten the file name for display
function [fileName] = parseName(fileName)
 for i=1:20
    if i==1
        [t,rem] = strtok(fileName,'\');
        fileName = strcat(t,'\...');
    elseif i>2
        [t,rem] = strtok(rem,'\');
        if isempty(t)==1
            break
        else
            fileName = strcat(fileName,'\',t);
        end
    else
        [t,rem] = strtok(rem,'\');
    end
 end

end

%% find where two curves intersect
function [x0,y0,iout,jout] = intersections(x1,y1,x2,y2,robust)
%INTERSECTIONS Intersections of curves.
%   Computes the (x,y) locations where two curves intersect.  The curves
%   can be broken with NaNs or have vertical segments.
%
% Example:
%   [X0,Y0] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% where X1 and Y1 are equal-length vectors of at least two points and
% represent curve 1.  Similarly, X2 and Y2 represent curve 2.
% X0 and Y0 are column vectors containing the points at which the two
% curves intersect.
%
% ROBUST (optional) set to 1 or true means to use a slight variation of the
% algorithm that might return duplicates of some intersection points, and
% then remove those duplicates.  The default is true, but since the
% algorithm is slightly slower you can set it to false if you know that
% your curves don't intersect at any segment boundaries.  Also, the robust
% version properly handles parallel and overlapping segments.
%
% The algorithm can return two additional vectors that indicate which
% segment pairs contain intersections and where they are:
%
%   [X0,Y0,I,J] = intersections(X1,Y1,X2,Y2,ROBUST);
%
% For each element of the vector I, I(k) = (segment number of (X1,Y1)) +
% (how far along this segment the intersection is).  For example, if I(k) =
% 45.25 then the intersection lies a quarter of the way between the line
% segment connecting (X1(45),Y1(45)) and (X1(46),Y1(46)).  Similarly for
% the vector J and the segments in (X2,Y2).
%
% You can also get intersections of a curve with itself.  Simply pass in
% only one curve, i.e.,
%
%   [X0,Y0] = intersections(X1,Y1,ROBUST);
%
% where, as before, ROBUST is optional.

% Version: 1.12, 27 January 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Theory of operation:
%
% Given two line segments, L1 and L2,
%
%   L1 endpoints:  (x1(1),y1(1)) and (x1(2),y1(2))
%   L2 endpoints:  (x2(1),y2(1)) and (x2(2),y2(2))
%
% we can write four equations with four unknowns and then solve them.  The
% four unknowns are t1, t2, x0 and y0, where (x0,y0) is the intersection of
% L1 and L2, t1 is the distance from the starting point of L1 to the
% intersection relative to the length of L1 and t2 is the distance from the
% starting point of L2 to the intersection relative to the length of L2.
%
% So, the four equations are
%
%    (x1(2) - x1(1))*t1 = x0 - x1(1)
%    (x2(2) - x2(1))*t2 = x0 - x2(1)
%    (y1(2) - y1(1))*t1 = y0 - y1(1)
%    (y2(2) - y2(1))*t2 = y0 - y2(1)
%
% Rearranging and writing in matrix form,
%
%  [x1(2)-x1(1)       0       -1   0;      [t1;      [-x1(1);
%        0       x2(2)-x2(1)  -1   0;   *   t2;   =   -x2(1);
%   y1(2)-y1(1)       0        0  -1;       x0;       -y1(1);
%        0       y2(2)-y2(1)   0  -1]       y0]       -y2(1)]
%
% Let's call that A*T = B.  We can solve for T with T = A\B.
%
% Once we have our solution we just have to look at t1 and t2 to determine
% whether L1 and L2 intersect.  If 0 <= t1 < 1 and 0 <= t2 < 1 then the two
% line segments cross and we can include (x0,y0) in the output.
%
% In principle, we have to perform this computation on every pair of line
% segments in the input data.  This can be quite a large number of pairs so
% we will reduce it by doing a simple preliminary check to eliminate line
% segment pairs that could not possibly cross.  The check is to look at the
% smallest enclosing rectangles (with sides parallel to the axes) for each
% line segment pair and see if they overlap.  If they do then we have to
% compute t1 and t2 (via the A\B computation) to see if the line segments
% cross, but if they don't then the line segments cannot cross.  In a
% typical application, this technique will eliminate most of the potential
% line segment pairs.


% Input checks.
error(nargchk(2,5,nargin))

% Adjustments when fewer than five arguments are supplied.
switch nargin
	case 2
		robust = true;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 3
		robust = x2;
		x2 = x1;
		y2 = y1;
		self_intersect = true;
	case 4
		robust = true;
		self_intersect = false;
	case 5
		self_intersect = false;
end

% x1 and y1 must be vectors with same number of points (at least 2).
if sum(size(x1) > 1) ~= 1 || sum(size(y1) > 1) ~= 1 || ...
		length(x1) ~= length(y1)
	error('X1 and Y1 must be equal-length vectors of at least 2 points.')
end
% x2 and y2 must be vectors with same number of points (at least 2).
if sum(size(x2) > 1) ~= 1 || sum(size(y2) > 1) ~= 1 || ...
		length(x2) ~= length(y2)
	error('X2 and Y2 must be equal-length vectors of at least 2 points.')
end


% Force all inputs to be column vectors.
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);

% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
xy1 = [x1 y1];
xy2 = [x2 y2];
dxy1 = diff(xy1);
dxy2 = diff(xy2);

% Determine the combinations of i and j where the rectangle enclosing the
% i'th line segment of curve 1 overlaps with the rectangle enclosing the
% j'th line segment of curve 2.
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
	repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
	repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
	repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
	repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
	repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
	repmat(min(y2(1:end-1),y2(2:end)).',n1,1));

% Force i and j to be column vectors, even when their length is zero, i.e.,
% we want them to be 0-by-1 instead of 0-by-0.
i = reshape(i,[],1);
j = reshape(j,[],1);

% Find segments pairs which have at least one vertex = NaN and remove them.
% This line is a fast way of finding such segment pairs.  We take
% advantage of the fact that NaNs propagate through calculations, in
% particular subtraction (in the calculation of dxy1 and dxy2, which we
% need anyway) and addition.
% At the same time we can remove redundant combinations of i and j in the
% case of finding intersections of a line with itself.
if self_intersect
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2)) | j <= i + 1;
else
	remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
end
i(remove) = [];
j(remove) = [];

% Initialize matrices.  We'll put the T's and B's in matrices and use them
% one column at a time.  AA is a 3-D extension of A where we'll use one
% plane at a time.
n = length(i);
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1(i,:).';
AA([2 4],2,:) = dxy2(j,:).';
B = -[x1(i) x2(j) y1(i) y2(j)].';

% Loop through possibilities.  Trap singularity warning and then use
% lastwarn to see if that plane of AA is near singular.  Process any such
% segment pairs to determine if they are colinear (overlap) or merely
% parallel.  That test consists of checking to see if one of the endpoints
% of the curve 2 segment lies on the curve 1 segment.  This is done by
% checking the cross product
%
%   (x1(2),y1(2)) - (x1(1),y1(1)) x (x2(2),y2(2)) - (x1(1),y1(1)).
%
% If this is close to zero then the segments overlap.

% If the robust option is false then we assume no two segment pairs are
% parallel and just go ahead and do the computation.  If A is ever singular
% a warning will appear.  This is faster and obviously you should use it
% only when you know you will never have overlapping or parallel segment
% pairs.

if robust
	overlap = false(n,1);
	warning_state = warning('off','MATLAB:singularMatrix');
	% Use try-catch to guarantee original warning state is restored.
	try
		lastwarn('')
		for k = 1:n
			T(:,k) = AA(:,:,k)\B(:,k);
			[unused,last_warn] = lastwarn;
			lastwarn('')
			if strcmp(last_warn,'MATLAB:singularMatrix')
				% Force in_range(k) to be false.
				T(1,k) = NaN;
				% Determine if these segments overlap or are just parallel.
				overlap(k) = rcond([dxy1(i(k),:);xy2(j(k),:) - xy1(i(k),:)]) < eps;
			end
		end
		warning(warning_state)
	catch err
		warning(warning_state)
		rethrow(err)
	end
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) <= 1 & T(2,:) <= 1).';
	% For overlapping segment pairs the algorithm will return an
	% intersection point that is at the center of the overlapping region.
	if any(overlap)
		ia = i(overlap);
		ja = j(overlap);
		% set x0 and y0 to middle of overlapping region.
		T(3,overlap) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
			min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1)))).'/2;
		T(4,overlap) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
			min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1)))).'/2;
		selected = in_range | overlap;
	else
		selected = in_range;
	end
	xy0 = T(3:4,selected).';
	
	% Remove duplicate intersection points.
	[xy0,index] = unique(xy0,'rows');
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		sel_index = find(selected);
		sel = sel_index(index);
		iout = i(sel) + T(1,sel).';
		jout = j(sel) + T(2,sel).';
	end
else % non-robust option
	for k = 1:n
		[L,U] = lu(AA(:,:,k));
		T(:,k) = U\(L\B(:,k));
	end
	
	% Find where t1 and t2 are between 0 and 1 and return the corresponding
	% x0 and y0 values.
	in_range = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
	x0 = T(3,in_range).';
	y0 = T(4,in_range).';
	
	% Compute how far along each line segment the intersections are.
	if nargout > 2
		iout = i(in_range) + T(1,in_range).';
		jout = j(in_range) + T(2,in_range).';
	end
end

% Plot the results (useful for debugging).
% plot(x1,y1,x2,y2,x0,y0,'ok');

end

