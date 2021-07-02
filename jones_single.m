% Basic Jones matrix routine for ULH phase-modulator model
% Steve J Elston, March 2020, updated July 2020
% based on creating jones-matrices for the quarter-wave plates, the LC
% device, and the mirror and then multiplying them together

% version with normalised tilt and retardance behaviour

% version to do +ve and -ve in continuous sweep.

% version modified so that quarter-wave plates can be varied.

clc
clear
close all

Amax=5.8*pi/180; %Angle amplitude (max tilt)
B=2*pi/180; % cubic factor in tilt vs voltage -
Ret_ini=0.856; %Initial retardance in waves
Ret_red=0.25; %Fractional reduction in retardance

orient=-45*pi/180;%+pi/4; % orientation of device 0-180 

% set orientations and retardances of waveplates
A1=pi/4; % +- pi/2
phi1=pi/2;

A2=-pi/4;% +- pi/2
phi2=pi/2;

% set up Jones matrix for first wave plate at angle A1: WP1
WP1=jones_matrix(A1,phi1);
% negative/mirror version of above.....................
WP1neg = jones_matrix(-1*A1,phi1);

% set up Jones matrix for second wave plate at angle A2: WP2
WP2=jones_matrix(A2,phi2);
% negative/mirror version of above.....................
WP2neg = jones_matrix(-1*A2,phi2);

% set up Jones matrix for the mirror, here the mirror is called: WP3
% this is based on treating the mirror like a half-wave plate with no
% reorientaiton of the opitc axis - the jones-matrices are the same
WP3=jones_matrix(0,pi);

% set up Jones matrix for a half wave plate at angle A (if used): WP4
% WP4=jones_matrix(pi/4,pi);

% now go on to original loop to run through voltage/tilt range...
% that section is same as before, but "orient" for the LC optic axis

count_i=0; % a handy counter for indexing vectors as the model runs....
range=-1:0.001:1;
tic
for V=range% normalise volts
    % range=V;% "range" is normalised to go from -1 to 1 as loop runs
    count_i=count_i+1; % increment counter as loop runs
    
    A=V*(Amax+B)-B*V^3+orient; % LC ULH tilt offset by orient angle
    % M_A=-1*A; % mirror angle of optic axis
    % THIS LINE TO BE CHANGED!
    RET=Ret_ini*(1-Ret_red*((A-orient)/Amax)^2); % retardance scale factor with tilt as loop runs 
    
    % calc Jones matrix for LC layer(s) for positive and mirrored angles
    % angles are "A" and "M_A", matrices are called: WP and M_WP
    phi=2*pi*RET; % this is the retardance of the LC layer, 2pi scaled by "RET"
    
    WP = jones_matrix(A,phi);
    M_WP = jones_matrix(-1*A,phi);
    
    % set up input light polarisation state vector (here in the x-direction)
    Ein = [1;0];
    %Ein = [0;1]; % use this to set input light polarisation in y-direction
    
    %%%%%%%%%% this is the main calulation line %%%%%%%%%%%%%%%%%
    Eout = WP1neg*M_WP*WP2neg*WP3*WP2*WP*WP1*Ein; % calc Eout for full system for +ve angle
    %%%%%%%%%% this is the main calulation line %%%%%%%%%%%%%%%%%
    
    % note: in stuff below use Eout(1) for x-polarisation and Eout(2) for y-polarisation
    tilt(count_i)=A; % load optic axis angle "A" into array called "tilt"
    ret(count_i)=RET; % load retardance scale factor into array called "ret"
    ang(count_i)=angle(Eout(1)); %phase angle of Eout_x into "angle" array
    %load intensities and amplitudes of Eout_x into arrays
    intensity(count_i)=abs(Eout(1)).^2;
    amplitude(count_i)=abs(Eout(1));
    % create handy array of results for data-theory comparisons
    result(count_i,:)=[V,A-orient,RET,ang(count_i),intensity(count_i)];
    % columns are: normalise volts, director tilt in radians, retardance in waves, phase mod angle, intensity mod value
    
end % end of voltage loop

% unwrap angle in results array to prevent weird jumps in the angle
result(:,4)=unwrap(result(:,4));
ang=unwrap(ang);
% plot the ploar plot to compare with spreadsheet model
figure
polarplot(ang,amplitude,'LineWidth',3);

figure
plot(range,ang,'LineWidth',3);
title('Phase');

figure
plot(range,amplitude,'LineWidth',3);
title('Amp');

toc

