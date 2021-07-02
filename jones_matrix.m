function [WP] = jones_matrix(A,phi)
%JONES_MATRIX Jones Matrix Generator
%	A   ----    This is the orientation angle of the waveplate
%	phi	----    This is the retardance of the waveplate (pi/2 for 1/4-wave)
RM = [cos(A) -sin(A); sin(A) cos(A)]; % rotation matrix
MRM = [cos(-A) -sin(-A); sin(-A) cos(-A)]; % invesrse rotation matrix
WP = [1 0; 0 exp(1i*phi)]; % unrotated waveplate matrix
WP = MRM*WP*RM; % create jones matrix of waveplate at angle "A"
end

