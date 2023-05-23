clc; 
clear variables; 
close all;
warning off;

% Script to create ASWING file from Prosperity Excel tool

% 'y' used as input for "yes" in all applicable variables

%% Load data from Excel 

filename = '20230502_ASWING_Prosperity_MASTER_Wing.xlsx';
structure_sheet = 'ASWING STRUCT';
aero_sheet = 'ASWING AERO';

%Inputs:
Outputfilename=['Ptest']; % Define .run File Name

opts = detectImportOptions(filename,'Sheet',structure_sheet,'Range','D1:F27');
opts = setvartype(opts,'Inputs','string');
constants = readtable(filename,opts,'Sheet',structure_sheet);
wingdata = readtable(filename,'Sheet',structure_sheet,'Range','AV49:CL100');
wingaero = readtable(filename,'Sheet',aero_sheet,'Range','A23:R50');
canarddata = readtable(filename,'Sheet',structure_sheet,'Range','CM49:EC100');
canardaero = readtable(filename,'Sheet',aero_sheet,'Range','S23:AL50');
vertstabilizer1data = readtable(filename,'Sheet',structure_sheet,'Range','ED49:FT100');
verticalstabilizer1aero = readtable(filename,'Sheet',aero_sheet,'Range','AM23:BC50');
vertstabilizer2data = readtable(filename,'Sheet',structure_sheet,'Range','FU49:HK100');
verticalstabilizer2aero = readtable(filename,'Sheet',aero_sheet,'Range','BD23:BT50');
boomdata{1} = readtable(filename,'Sheet',structure_sheet,'Range','Q49:U100');
boomdata{2} = readtable(filename,'Sheet',structure_sheet,'Range','V49:Z100');
boomdata{3} = readtable(filename,'Sheet',structure_sheet,'Range','AA49:AE100'); 
boomdata{4} = readtable(filename,'Sheet',structure_sheet,'Range','AF49:AJ100');
boomdata{5} = readtable(filename,'Sheet',structure_sheet,'Range','AK49:AO100');
boomdata{6} = readtable(filename,'Sheet',structure_sheet,'Range','AP49:AT100');
jointdata = readtable(filename,'Sheet',structure_sheet,'Range','K49:O59');
massdata = readtable(filename,'Sheet',structure_sheet,'Range','A51:J100');

%Note - Mass-less engines are hard-coded in this MATLAB file for flutter runs.  
%       This may be moved to Excel INPUT in the future.  See rows ~400.


% eliminating NaN values
jointdata=rmmissing(jointdata);
wingdata=rmmissing(wingdata);
wingaero=rmmissing(wingaero);
canarddata=rmmissing(canarddata);
canardaero=rmmissing(canardaero);
vertstabilizer1data=rmmissing(vertstabilizer1data);
verticalstabilizer1aero=rmmissing(verticalstabilizer1aero);
vertstabilizer2data=rmmissing(vertstabilizer2data);
verticalstabilizer2aero=rmmissing(verticalstabilizer2aero);
boomdata{1}=rmmissing(boomdata{1});
boomdata{2}=rmmissing(boomdata{2});
boomdata{3}=rmmissing(boomdata{3});
boomdata{4}=rmmissing(boomdata{4});
boomdata{5}=rmmissing(boomdata{5});
boomdata{6}=rmmissing(boomdata{6});
massdata=rmmissing(massdata);

%%
stiffness_factor = 1; % factor is applied to all values of EI and GJ
% and is used for flutter debugging. Nominally set this to 1.

%------------------------------------------------------------------------
default=zeros(1000,1); %default vector for values

%Name
aircraftname=string(constants.Inputs(find(contains(constants.Labels,"Name")))); %Name of configuration

%Units
lengthU='1.0 ' + string(constants.Inputs(find(contains(constants.Labels,"Length")))); %Length reference unit
timeU='1.0 ' + string(constants.Inputs(find(contains(constants.Labels,"Time")))); %time reference unit
massU='1.0 ' + string(constants.Inputs(find(contains(constants.Labels,"Mass")))); %mass reference unit

%Constants
g= double(constants.Inputs(find(contains(constants.Labels,"gravity")))); %gravity constant
rhoSL=string(constants.Inputs(find(contains(constants.Labels,"air density")))); %air density sea level
a_SL=string(constants.Inputs(find(contains(constants.Labels,"speed of sound")))); % speed of sound sea level

%Reference
S=string(constants.Inputs(find(contains(constants.Labels,"wing area")))); % wing area
c=string(constants.Inputs(find(contains(constants.Labels,"reference chord")))); %reference chord
b=string(constants.Inputs(find(contains(constants.Labels,"wing span")))); %wing span
referencelocation='y';
momentXYZ(1) = double(constants.Inputs(find(contains(constants.Labels,"reference x loc.")))); % reference x location
momentXYZ(2) = double(constants.Inputs(find(contains(constants.Labels,"reference y loc.")))); % reference y location
momentXYZ(3) = double(constants.Inputs(find(contains(constants.Labels,"reference z loc.")))); % reference z location
accXYZ=str2num(char(constants.Inputs(find(contains(constants.Labels,"ref. Accelerations"))))); %reference location for accelerations [x,y,z], default [0,0,0]
vXYZ=str2num(char(constants.Inputs(find(contains(constants.Labels,"ref. Velocities"))))); %reference location for velocities [x,y,z], default [0,0,0]

%Engine
enginevari = 'y'; %include engines? 'y' for yes

%Ground
NbeamG=[1]; %index beams to ground
tG=[0]; %points of beams to ground
KGtype=[0]; %0=rigid (default), 1=pinned, 2=roller

%Strut
strutvari='n'; %include struts? 'y' for yes
Nbeamstr=[1,1]; %index beams for struts
tstr=[-28,28]; %location on indexed beams
Xpstr=[0,0]; %x location for other pylon endpoint
Ypstr=[0,0]; %y location for other pylon endpoint
Zpstr=[0,0]; %z location for other pylon endpoint
Xwstr=[0,0]; %x location for anchor endpoint
Ywstr=[0,0]; %y location for anchor endpoint
Zwstr=[0,0]; %z location for anchor endpoint
dLstr=[0,0]; %strut slack
EAwstr=[0,0]; %extensional stiffness

%Wing
wingvari='y'; %include wing? 'y' for yes
tgeow=100 * (wingaero.Yle - min(wingaero.Yle))/(max(wingaero.Yle)- min(wingaero.Yle)); %geoemtry indexing
xgeow=wingaero.X_S; %x geometry values
ygeow=wingaero.Yle; %y geometry values
zgeow=wingaero.Zle; %z geometry values
chordgeow=wingaero.Chord_asw; %chord values
twistgeow=wingaero.Twist; %twist values
Xaxgeow=wingaero.Xax; %distance/chord of s-axis from leading edge (Default=0.5)

tshellw=wingdata.t; %evaluation indexing
Cshellw=wingdata.Cshell; %c-Cea where strain is evaluated
Nshellw=wingdata.Nshell; %n-Nea where strain is evaluated
Atshellw=wingdata.Atshell; %A*T for shear stress evaluation
Cmw=wingdata.Cm; % section pitching moment coefficient about chord/4
dCLdaw = wingdata.dCLda;

tmatw=wingdata.t; %material indexing
mgw=wingdata.mg; %weight/span
EIccw=wingdata.EIcc*stiffness_factor; %bending stiffness about c-axis (out-of-plane) (Default=Inf)
EInnw=wingdata.EInn*stiffness_factor; %bending stiffness about n-axis (in-plane) (Default=Inf)
GJw=wingdata.GJ*stiffness_factor; %torsional stiffness (Default=Inf)
Ccgw=wingdata.Ccg; %c-location of section masscentroid

taerow=tgeow; %aero indexing
CLmaxw=wingaero.Clmax; %section CLmax
CLminw=wingaero.Clmin; %section CLmin
Cdfw=wingaero.Cdf; %section friction drag
Cdpw=wingaero.Cdp; %section pressure drag
alphaw=wingaero.LocalAlpha; %angle of zero-lift above c-axis

% Canard
canardvari='y'; %include canard? 'y' for yes
tgeoc=100 * (canardaero.Yle - min(canardaero.Yle))/(max(canardaero.Yle)- min(canardaero.Yle)); %geoemtry indexing
xgeoc=canardaero.X_S; %x geometry values
ygeoc=canardaero.Yle; %y geometry values
zgeoc=canardaero.Zle; %z geometry values
chordgeoc=canardaero.Chord_asw; %chord values
twistgeoc=canardaero.Twist; %twist values
Xaxgeoc=canardaero.Xax; %distance/chord of s-axis from leading edge (Default=0.5)

tshellc=canarddata.t; %evaluation indexing
Cshellc=canarddata.Cshell; %c-Cea where strain is evaluated
Nshellc=canarddata.Nshell; %n-Nea where strain is evaluated
Atshellc=canarddata.Atshell; %A*T for shear stress evaluation
Cmc=canarddata.Cm; % section pitching moment coefficient about chord/4
dCLdac = canarddata.dCLda;

tmatc=canarddata.t; %material indexing
mgc=canarddata.mg; %weight/span
EIccc=canarddata.EIcc*stiffness_factor; %bending stiffness about c-axis (out-of-plane) (Default=Inf)
EInnc=canarddata.EInn*stiffness_factor; %bending stiffness about n-axis (in-plane) (Default=Inf)
GJc=canarddata.GJ*stiffness_factor; %torsional stiffness (Default=Inf)
Ccgc=canarddata.Ccg; %c-location of section masscentroid

taeroc=tgeoc; %aero indexing
CLmaxc=canardaero.Clmax; %section CLmax
CLminc=canardaero.Clmin; %section CLmin
Cdfc=canardaero.Cdf; %section friction drag
Cdpc=canardaero.Cdp; %section pressure drag
alphac=canardaero.LocalAlpha; %angle of zero-lift above c-axis

tcontrolc = canardaero.Yle/(max(canardaero.Yle)-min(canardaero.Yle))*100;
CLdF1c =canardaero.dCLdf;
CMdF1c =canardaero.dCMdf;

% Vertical Stabilizer 1
vs1vari='y'; %include vertical stabilizer 1? 'y' for yes
tgeovs1= 1 +99 * (max(verticalstabilizer1aero.Zle) - verticalstabilizer1aero.Zle)/(max(verticalstabilizer1aero.Zle)- min(verticalstabilizer1aero.Zle)); %geoemtry indexing
xgeovs1=verticalstabilizer1aero.X_S; %x geometry values
ygeovs1=verticalstabilizer1aero.Yle; %y geometry values
zgeovs1=verticalstabilizer1aero.Zle; %z geometry values
chordgeovs1=verticalstabilizer1aero.Chord_asw; %chord values
twistgeovs1=verticalstabilizer1aero.Twist; %twist values
Xaxgeovs1=verticalstabilizer1aero.Xax; %distance/chord of s-axis from leading edge (Default=0.5)

tshellvs1=vertstabilizer1data.t; %evaluation indexing
Cshellvs1=vertstabilizer1data.Cshell; %c-Cea where strain is evaluated
Nshellvs1=vertstabilizer1data.Nshell; %n-Nea where strain is evaluated
Atshellvs1=vertstabilizer1data.Atshell; %A*T for shear stress evaluation
Cmvs1=vertstabilizer1data.Cm; % section pitching moment coefficient about chord/4
dCLdavs1 = vertstabilizer1data.dCLda;

tmatvs1=vertstabilizer1data.t; %material indexing
mgvs1=vertstabilizer1data.mg; %weight/span
EIccvs1=vertstabilizer1data.EIcc*stiffness_factor; %bending stiffness about c-axis (out-of-plane) (Default=Inf)
EInnvs1=vertstabilizer1data.EInn*stiffness_factor; %bending stiffness about n-axis (in-plane) (Default=Inf)
GJvs1=vertstabilizer1data.GJ*stiffness_factor; %torsional stiffness (Default=Inf)
Ccgvs1=vertstabilizer1data.Ccg; %c-location of section masscentroid

taerovs1=tgeovs1; %aero indexing
CLmaxvs1=verticalstabilizer1aero.Clmax; %section CLmax
CLminvs1=verticalstabilizer1aero.Clmin; %section CLmin
Cdfvs1=verticalstabilizer1aero.Cdf; %section friction drag
Cdpvs1=verticalstabilizer1aero.Cdp; %section pressure drag
alphavs1=verticalstabilizer1aero.LocalAlpha; %angle of zero-lift above c-axis

% Vertical Stabilizer 2
vs2vari='y'; %include vertical stabilizer 2? 'y' for yes
tgeovs2= 1 +99 * (max(verticalstabilizer2aero.Zle) - verticalstabilizer2aero.Zle)/(max(verticalstabilizer2aero.Zle)- min(verticalstabilizer2aero.Zle)); %geoemtry indexing
xgeovs2=verticalstabilizer2aero.X_S; %x geometry values
ygeovs2=verticalstabilizer2aero.Yle; %y geometry values
zgeovs2=verticalstabilizer2aero.Zle; %z geometry values
chordgeovs2=verticalstabilizer2aero.Chord_asw; %chord values
twistgeovs2=verticalstabilizer2aero.Twist; %twist values
Xaxgeovs2=verticalstabilizer2aero.Xax; %distance/chord of s-axis from leading edge (Default=0.5)

tshellvs2=vertstabilizer2data.t; %evaluation indexing
Cshellvs2=vertstabilizer2data.Cshell; %c-Cea where strain is evaluated
Nshellvs2=vertstabilizer2data.Nshell; %n-Nea where strain is evaluated
Atshellvs2=vertstabilizer2data.Atshell; %A*T for shear stress evaluation
Cmvs2=vertstabilizer2data.Cm; % section pitching moment coefficient about chord/4   
dCLdavs2 = vertstabilizer2data.dCLda;

tmatvs2=vertstabilizer2data.t; %material indexing   
mgvs2=vertstabilizer2data.mg; %weight/span
EIccvs2=vertstabilizer2data.EIcc*stiffness_factor; %bending stiffness about c-axis (out-of-plane) (Default=Inf)
EInnvs2=vertstabilizer2data.EInn*stiffness_factor; %bending stiffness about n-axis (in-plane) (Default=Inf)
GJvs2=vertstabilizer2data.GJ*stiffness_factor; %torsional stiffness (Default=Inf)
Ccgvs2=vertstabilizer2data.Ccg; %c-location of section masscentroid

taerovs2=tgeovs2; %aero indexing
CLmaxvs2=verticalstabilizer2aero.Clmax; %section CLmax
CLminvs2=verticalstabilizer2aero.Clmin; %section CLmin
Cdfvs2=verticalstabilizer2aero.Cdf; %section friction drag
Cdpvs2=verticalstabilizer2aero.Cdp; %section pressure drag
alphavs2=verticalstabilizer2aero.LocalAlpha; %angle of zero-lift above c-axis

%Booms
for i = 1:6
    bvari(i) = 'y';
    xgeob(i,:)=boomdata{i}.X; %x geometry values
    ygeob(i,:)=boomdata{i}.Y; %y geometry values
    zgeob(i,:)=boomdata{i}.Z; %z geometry values
    radiusb(i,:)=boomdata{i}.Radius; %radius geometry values
    twistgeob(i,:)=boomdata{i}.Twist; %twist geometry value
    tb(i,:) = boomdata{i}.X;
    EInnb(i,:) = [1,1]*double(constants.Inputs(find(contains(constants.Labels,"EInn"))))*stiffness_factor;
    EIccb(i,:) = [1,1]*double(constants.Inputs(find(contains(constants.Labels,"EIcc"))))*stiffness_factor;
    GJb(i,:) = [1,1]*double(constants.Inputs(find(contains(constants.Labels,"GJ"))))*stiffness_factor;
    mgb(i,:) = [1,1]*double(constants.Inputs(find(contains(constants.Labels,"Mg"))));
end

%Joint
jointvari='y'; %specify joints? 'y' for yes
Nbeam1=jointdata.Boom1; %index of beam 1
Nbeam2=jointdata.Boom2; %index of beam 2
t1j=jointdata.t1; %joint location on beam 1
t2j=jointdata.t2; %joint location on beam 2
KJtype=jointdata.Type; %joint type 0=rigid (default), 1=pinned, 2=roller, 3=sprung-hinge

%Jangle parameters are required for joint type 3's - otherwise ignored

%Jangle
Njoint=[2,3]; %index joints for joint angle (sequential)
hx=[0,0]; %hinge axis vector (x) (sequential)
hy=[0,0]; %hinge axis vector (y) (sequential)
hz=[0,0]; %hinge axis vector (z) (sequential)
Momh=[0,0]; %hinge moments
Angh=[1,1]; %hinge angles

%Weight
weightvari='y'; %include weights? 'y' for yes
Xow=massdata.X ; %x location of point mass
Yow=massdata.Y ; %y location of point mass
Zow=massdata.Z ; %z location of point mass
Weight=massdata.M ; % weight of point mass
Nm = massdata.Number; % connection beam
tw = massdata.tloc; % connection location on boom
for i = 1:length(Xow)
    CDAw(i)=0; %drag area of point mass
    Volw(i)=0; %volume of point mass
    Hxow(i)=0; %x angular momentum
    Hyow(i)=0; %y angular momentum
    Hzow(i)=0; %z angular momentum
    Ixxw(i)=0; %MoI xx
    Iyyw(i)=0; %MoI yy
    Izzw(i)=0; %MoI zz
    Ixyw(i)=0; %MoI xy
    Ixzw(i)=0; %MoI xz
    Iyzw(i)=0; %MoI yz
end

% %Fuselage
% fusevari='y'; %include H-stab? 'y' for yes
% tgeof=[0,0]; %geoemtry indexing
% xgeof=[0,0]; %x geometry values
% radiusf=[0,0]; %radius of fuselage
% 
% tmatf=[0,0]; %material indexing
% mgf=[0,0]; %weight/span
% EIccf=[0,0]; %bending stiffness about c-axis (out-of-plane) (Default=Inf)
% EInnf=[0,0]; %bending stiffness about n-axis (in-plane) (Default=Inf)
% GJf=[0,0]; %torsional stiffness (Default=Inf)
% 
% taerof=[0,0]; %aero indexing
% Cdff=[0,0]; %section friction drag
% Cdpf=[0,0]; %section pressure drag
%-------------------------------------------------------------------------

% accounting for mass inputs

if massU == "1.0 kg"
    Weight = Weight * g;
    mgb = mgb * g;
    mgw = mgw * g;
    mgc = mgc * g;
    mgvs1 = mgvs1 * g;
    mgvs2 = mgvs2 * g;
end

fid_out = fopen(append(Outputfilename,'.asw'),'w');

fprintf(fid_out, '#===========================\n');
fprintf(fid_out, 'Name\n');
fprintf(fid_out, '%s\n',aircraftname);
fprintf(fid_out, 'End\n');
fprintf(fid_out, '#===========================\n');
fprintf(fid_out, 'Units\n');
fprintf(fid_out, 'L %s\n',lengthU);
fprintf(fid_out, 'T %s\n',timeU);
fprintf(fid_out, 'M %s\n',massU);
fprintf(fid_out, 'End\n');
fprintf(fid_out, '#===========================\n');
fprintf(fid_out, 'Constant\n');
fprintf(fid_out, '#   g    rho_SL     V_sl_sound\n');
fprintf(fid_out, '%f %f %f\n',g,rhoSL,a_SL);
fprintf(fid_out, 'End\n');
fprintf(fid_out, '#===========================\n');
fprintf(fid_out, 'Reference\n');
fprintf(fid_out, '# Sref  Cref  Bref\n');
fprintf(fid_out, '%f %f %f\n',S,c,b);
if referencelocation=='y'
    fprintf(fid_out, '#  Xmom  Ymom  Zmom\n');
    fprintf(fid_out, '%f %f %f\n',momentXYZ(1),momentXYZ(2),momentXYZ(3));
%     fprintf(fid_out, '%f %f %f\n',accXYZ(1),accXYZ(2),accXYZ(3));
%     fprintf(fid_out, '%f %f %f\n',vXYZ(1),vXYZ(2),vXYZ(3));
end
fprintf(fid_out, 'End\n');
fprintf(fid_out, '#===========================\n');
if jointvari=='y'
    fprintf(fid_out, 'Joint\n');
    fprintf(fid_out, '#  Nbeam1  Nbeam2    t1     t2    [ KJtype ]\n');
    for i=1:length(Nbeam1)
        fprintf(fid_out, '%f %f %f %f %f\n',Nbeam1(i),Nbeam2(i),t1j(i),t2j(i),KJtype(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
    if any(KJtype==3)
        for i=1:length(Njoint)
        fprintf(fid_out, 'Jangle\n');
        fprintf(fid_out, '# Njoint  hx   hy   hz\n');
        fprintf(fid_out, '%f %f %f %f\n',Njoint(i),hx(i),hy(i),hz(i));
        fprintf(fid_out, '#  Momh    Angh\n');
            for j=1:length(Momh)
                fprintf(fid_out, '%f   %f\n',Momh(j),Angh(j));
            end
        fprintf(fid_out, 'End\n');
        fprintf(fid_out, '#===========================\n');
        end
    end
end

fprintf(fid_out, 'Ground\n');
fprintf(fid_out, '#  Nbeam    t    [ KGtype ]\n');
for i=1:length(NbeamG)
    fprintf(fid_out, '%f %f %f\n',NbeamG(i),tG(i),KGtype(i));
end
fprintf(fid_out, 'End\n');
fprintf(fid_out, '#===========================\n');
if strutvari=='y'
    fprintf(fid_out, 'Strut\n');
    fprintf(fid_out, '# Nbeam   t     Xp    Yp    Zp     Xw    Yw    Zw     dL     EAw\n');
    for i=1:length(Nbeamstr)
        fprintf(fid_out, '%f %f %f %f %f %f %f %f %f %f\n',...
            Nbeamstr(i),tstr(i),Xpstr(i),Ypstr(i),Zpstr(i),Xwstr(i),Ywstr(i),Zwstr(i),dLstr(i),EAwstr(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end
if weightvari=='y'
    fprintf(fid_out, 'Weight\n');
    fprintf(fid_out, '# Nbeam  t    Xo   Yo   Zo    Weight   CDA    Vol    Hxo    Hyo    Hzo   [ Ixx   Iyy   Izz   Ixy   Ixz   Ixz ]\n');
    for i=1:length(tw)
        fprintf(fid_out, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',...
            Nm(i),tw(i),Xow(i),Yow(i),Zow(i),Weight(i),CDAw(i),Volw(i),Hxow(i),Hyow(i),...
            Hzow(i),Ixxw(i),Iyyw(i),Izzw(i),Ixyw(i),Ixzw(i),Iyzw(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end
if enginevari == 'y'
    fprintf(fid_out, 'Engine\n');
    fprintf(fid_out, '# KPeng IEtyp Nbeam  t     Xp   Yp   Zp   Tx  Ty  Tz  dF/dPe dM/dPe  Rdisk\n');
    fprintf(fid_out, '    1    0     7     '+string(massdata.X(find(contains(massdata.Var1,"Pusher 2"))))+ ...
    '   '+string(massdata.X(find(contains(massdata.Var1,"Pusher 2"))))+'   '+...
    string(massdata.Y(find(contains(massdata.Var1,"Pusher 2"))))+'  '+ ...
    string(massdata.Z(find(contains(massdata.Var1,"Pusher 2"))))+...
    '  -1.  0.  0.   1.0   0.0     0.5\n');
    fprintf(fid_out, '    1    0     8     '+string(massdata.X(find(contains(massdata.Var1,"Pusher 1"))))+ ...
    '   '+string(massdata.X(find(contains(massdata.Var1,"Pusher 1"))))+'   '+...
    string(massdata.Y(find(contains(massdata.Var1,"Pusher 1"))))+'  '+ ...
    string(massdata.Z(find(contains(massdata.Var1,"Pusher 1"))))+...
    '  -1.  0.  0.   1.0   0.0     0.5\n');
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end
if wingvari=='y'
    fprintf(fid_out, 'Beam 1\n');
    fprintf(fid_out, 'Wing\n');
    fprintf(fid_out, 't x y z chord twist  Xax\n');
    for i=1:length(tgeow)
        fprintf(fid_out, '%f %f %f %f %f %f %f\n',tgeow(i),xgeow(i),ygeow(i),zgeow(i),...
            chordgeow(i),twistgeow(i),Xaxgeow(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't Cshell Nshell Atshell Cm dCLda\n');
    for i=1:length(tshellw)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tshellw(i),Cshellw(i),Nshellw(i),Atshellw(i),Cmw(i),dCLdaw(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't mg EIcc EInn GJ Ccg\n');
    for i=1:length(tmatw)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tmatw(i),mgw(i),EIccw(i),EInnw(i),GJw(i),Ccgw(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't CLmax CLmin Cdf Cdp alpha\n');
    for i=1:length(taerow)
        fprintf(fid_out, '%f %f %f %f %f %f\n',taerow(i),CLmaxw(i),CLminw(i),Cdfw(i),Cdpw(i),alphaw(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end


if canardvari=='y'
    fprintf(fid_out, 'Beam 2\n');
    fprintf(fid_out, 'Canard\n');
    fprintf(fid_out, 't x y z chord twist Xax\n');
    for i=1:length(tgeoc)
        fprintf(fid_out, '%f %f %f %f %f %f %f\n',tgeoc(i),xgeoc(i),ygeoc(i),zgeoc(i),...
            chordgeoc(i),twistgeoc(i),Xaxgeoc(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't Cshell Nshell Atshell Cm dCLda\n');
    for i=1:length(tshellc)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tshellc(i),Cshellc(i),Nshellc(i),Atshellc(i),Cmc(i),dCLdac(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't mg EIcc EInn GJ Ccg\n');
    for i=1:length(tmatc)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tmatc(i),mgc(i),EIccc(i),EInnc(i),GJc(i),Ccgc(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't CLmax CLmin Cdf Cdp alpha\n');
    for i=1:length(taeroc)
        fprintf(fid_out, '%f %f %f %f %f %f\n',taeroc(i),CLmaxc(i),CLminc(i),Cdfc(i),Cdpc(i),alphac(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't dCLdF1 dCMdF1\n');
    for i=1:length(tcontrolc)
        fprintf(fid_out, '%f %f %f\n',tcontrolc(i),CLdF1c(i),CMdF1c(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end
if vs1vari=='y'
    fprintf(fid_out, 'Beam 3\n');
    fprintf(fid_out, 'Vertical Stabilizer 1\n');
    fprintf(fid_out, 't x y z chord twist Xax\n');
    for i=1:length(tgeovs1)
        fprintf(fid_out, '%f %f %f %f %f %f %f\n',tgeovs1(i),xgeovs1(i),ygeovs1(i),zgeovs1(i),...
            chordgeovs1(i),twistgeovs1(i),Xaxgeovs1(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't Cshell Nshell Atshell Cm dCLda\n');
    for i=1:length(tshellvs1)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tshellvs1(i),Cshellvs1(i),Nshellvs1(i),Atshellvs1(i),Cmvs1(i),dCLdavs1(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't mg EIcc EInn GJ Ccg\n');
    for i=1:length(tmatvs1)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tmatvs1(i),mgvs1(i),EIccvs1(i),EInnvs1(i),GJvs1(i),Ccgvs1(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't CLmax CLmin Cdf Cdp alpha\n');
    for i=1:length(taerovs1)
        fprintf(fid_out, '%f %f %f %f %f %f\n',taerovs1(i),CLmaxvs1(i),CLminvs1(i),Cdfvs1(i),Cdpvs1(i),alphavs1(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end
if vs2vari=='y'
    fprintf(fid_out, 'Beam 4\n');
    fprintf(fid_out, 'Vertical Stabilizer 2\n');
    fprintf(fid_out, 't x y z chord twist Xax\n');
    for i=1:length(tgeovs2)
        fprintf(fid_out, '%f %f %f %f %f %f %f\n',tgeovs2(i),xgeovs2(i),ygeovs2(i),zgeovs2(i),...
            chordgeovs2(i),twistgeovs2(i),Xaxgeovs2(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't Cshell Nshell Atshell Cm dCLda\n');
    for i=1:length(tshellvs2)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tshellvs2(i),Cshellvs2(i),Nshellvs2(i),Atshellvs2(i),Cmvs2(i),dCLdavs2(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't mg EIcc EInn GJ Ccg\n');
    for i=1:length(tmatvs2)
        fprintf(fid_out, '%f %f %f %f %f %f\n',tmatvs2(i),mgvs2(i),EIccvs2(i),EInnvs2(i),GJvs2(i),Ccgvs2(i));
    end
    fprintf(fid_out, '\n');
    fprintf(fid_out, 't CLmax CLmin Cdf Cdp alpha\n');
    for i=1:length(taerovs2)
        fprintf(fid_out, '%f %f %f %f %f %f\n',taerovs2(i),CLmaxvs2(i),CLminvs2(i),Cdfvs2(i),Cdpvs2(i),alphavs2(i));
    end
    fprintf(fid_out, 'End\n');
    fprintf(fid_out, '#===========================\n');
end

for i = 1:length(tb)
    if bvari(i)=='y'
        fprintf(fid_out, 'Beam ' + string(i+4) + '\n');
        fprintf(fid_out, 'Structural Beam # ' + string(i) + '\n');
        fprintf(fid_out, 't x y z radius twist\n');
        for j=1:length(tb(i,:))
            fprintf(fid_out, '%f %f %f %f %f %f\n',tb(i,j),xgeob(i,j),ygeob(i,j),zgeob(i,j),...
                radiusb(i,j),twistgeob(i,j));
        end
        fprintf(fid_out, 't mg EIcc EInn GJ\n');
        for j=1:length(tb(i,:))
            fprintf(fid_out, '%f %f %f %f %f\n',tb(i,j),mgb(i,j),EIccb(i,j),EInnb(i,j),...
                GJb(i,j));
        end
    %     fprintf(fid_out, '\n');
    %     fprintf(fid_out, 't mg EIcc EInn GJ\n');
    %     for i=1:length(tmatvs2)
    %         fprintf(fid_out, '%f %f %f %f %f\n',tmatvs2(i),mgvs2(i),EIccvs2(i),EInnvs2(i),GJvs2(i));
    %  
    %     end
    %     fprintf(fid_out, '\n');
    %     fprintf(fid_out, 't CLmax CLmin Cdf Cdp alpha\n');
    %     for i=1:length(taerovs2)
    %         fprintf(fid_out, '%f %f %f %f %f %f\n',taerovs2(i),CLmaxvs2(i),CLminvs2(i),Cdfvs2(i),Cdpvs2(i),alphavs2(i));
    %     end
        fprintf(fid_out, 'End\n');
        fprintf(fid_out, '#===========================\n');
    end

end

fclose(fid_out);

disp(append('Configuration File: ',Outputfilename,'.asw Generated'));