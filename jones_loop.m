% Basic Jones matrix routine for ULH phase-modulator model
% Steve J Elston, March 2020, updated July 2020
% based on creating jones-matrices for the quarter-wave plates, the LC
% device, and the mirror and then multiplying them together

% version with normalised tilt and retardance behaviour

% version to do +ve and -ve in continuous sweep.

% version modified so that quarter-wave plates can be varied.

clc
clear
% close all

Timestamp=datestr(now,'yyyy-mm-dd_HH-MM-SS');

Amax=11.75*pi/180; %Angle amplitude (max tilt)
B=3*pi/180; % cubic factor in tilt vs voltage -
Ret_ini=0.856; %Initial retardance in waves
Ret_red=0.25; %Fractional reduction in retardance


count_i=0; % a handy counter for indexing vectors as the model runs....
range=-1:0.01:1;
device_orient=(0:0.5:180)*pi/180;%+pi/4; % orientation of device 0-180
phase_diff=pi();


% set orientations and retardances of waveplates
A1_range=-pi()/2:pi()/20:pi()/2; % +- pi/2
A2_range=-pi()/2:pi()/20:pi()/2;% +- pi/2

tic
for A1=A1_range
    for A2=A2_range
        phi1=pi/2;
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
        
        
        for orient=device_orient
            
            isFind=0;
            count_i=count_i+1;
            clc
            fprintf('%d/%d\t%.3f%%\t%.0fs\n',count_i,size(device_orient,2)*size(A1_range,2)*size(A2_range,2),count_i*100/(size(device_orient,2)*size(A1_range,2)*size(A2_range,2)),toc);
            count_j=0;
            for V=range% normalise volts
                % range=V;% "range" is normalised to go from -1 to 1 as loop runs
                count_j=count_j+1; % increment counter as loop runs
                
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
                tilt(count_j)=A; % load optic axis angle "A" into array called "tilt"
                ret(count_j)=RET; % load retardance scale factor into array called "ret"
                ang(count_j)=unwrap(angle(Eout(1))); %phase angle of Eout_x into "angle" array
                %load intensities and amplitudes of Eout_x into arra);ys
                intensity(count_j)=abs(Eout(1)).^2;
                amplitude(count_j)=abs(Eout(1));
                % create handy array of results for data-theory comparisons
                result(count_j,:)=[V,A-orient,RET,ang(count_j),intensity(count_j)];
                % columns are: normalise volts, director tilt in radians, retardance in waves, phase mod angle, intensity mod value
                
            end % end of voltage loop
            % plot the ploar plot to compare with spreadsheet model
            %%%%%%%%%%%%%%%%%% end of each combination%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            isFind=0;
            ang_max=max(ang);
            ang_min=min(ang);
            if (ang_max-ang_min)>phase_diff
                %%% find the section from ang %%%
                section_i=1;
                for i=1:1:size(ang,2)
                    if ang(i)+phase_diff<=ang_max || ang(i)-phase_diff>=ang_min
                        find_pos=find(ang(i:end)>=ang(i)+pi());
                        find_pos=find_pos+i-1;
                        find_neg=find(ang(i:end)<=ang(i)-pi());
                        find_neg=find_neg+i-1;
                        if ~isempty(find_pos)||~isempty(find_neg)
                            if ~isempty(find_pos)
                                for j=1:1:size(find_pos,2)
                                    section(section_i+j-1,:)=[i,find_pos(j)];
                                end
                            end
                            if ~isempty(find_neg)
                                for j=1:1:size(find_neg,2)
                                    section(section_i+j-1,:)=[i,find_neg(j)];
                                end
                            end
                            section_i=section_i+j;
                        end
                    end
                end
                
                %         if std(section(:,1))<15 || std(section(:,2))<15
                section=nearest(mean(section,1));
                %         end
                
                %%% check if the ang sections are linear (using linear fit, may not accurate) %%%
                for k=1:1:size(section,1)
                    fprintf('\tFitting...\r');
                    [~,ang_gof]=fit([range(section(k,1):section(k,2))]',[ang(section(k,1):section(k,2))]','poly1');
                    [amp_cfit,amp_gof]=fit([range(section(k,1):section(k,2))]',[amplitude(section(k,1):section(k,2))]','poly1');
                    section(k,3)=ang_gof.adjrsquare;
                    section(k,4)=amp_cfit.p1;
                    section(k,5)=amp_gof.adjrsquare;
                end
                section=section(section(:,3)>0.9 & section(:,5)>0.9 & abs(section(:,4))<0.2,:);
                if ~isempty(section)
                    section_found=section;             
                    isFind=1;
                    results_find(count_i,:)=[A1,orient,A2,isFind,range(section_found(1)),range(section_found(2)),section_found(3:5)];
                    fprintf('\tWE FOUND IT!!!!!!!!!!!!!!!!!!!!!!\r');
                    fprintf('\tV1=%.2f, V2=%.2f, ang_gof=%.4f, amp_p1=%.4f, amp_gof=%.4f\r',range(section_found(1)),range(section_found(2)),section_found(3),section_found(4),section_found(5));
                else        
                    isFind=0;
                    results_find(count_i,:)=[A1,orient,A2,isFind,[0,0,0,0,0]];
                end
            else
                isFind=0;
                results_find(count_i,:)=[A1,orient,A2,isFind,[0,0,0,0,0]];
            end
%             results_find(count_i,:)=[A1,orient,A2,isFind];
            
            if isFind==1
                figure(1)
                polarplot(ang,amplitude,'LineWidth',3);
                
                figure(2)
                plot(range,ang,'LineWidth',3);
                title('Phase');
                
                figure(3)
                plot(range,amplitude,'LineWidth',3);
                ylim([0 1]);
                title('Amp');
%                 pause
%                 ;
            end
            clear section;
            % about 60s with no plot display
        end
        % unwrap angle in results array to prevent weird jumps in the angle
    end
end
toc
save(['results_' Timestamp '.mat'],'results_find');

