%
% noe test for 1D potential V8 under microstate Markovian lagtime
%using jumping window
option='noe'; %modify, you can either choose 'noe' or 'bootstrap'
resultdir = 'noe_macro_formula_sliding/';
lagtime=1; %modify
if (exist(resultdir) == 0)
	mkdir(resultdir);
end

trajlist = importdata('trajlist_micro.txt');
traj_len=[];
traj_num  = length(trajlist);

mapping=importdata('Spec_Clus_4_state_mapping.txt')+1; %modify
for j = 1:traj_num
    temp = importdata(trajlist{j})+1; %starting from 1 now
    %mapping
    temp1 = temp;
    for k = 1:length(temp)
        temp1(k) = mapping(temp(k));
    end
    trajMacro{j} = temp1;
end

nMacro = 4;%modify

%for microstates
tCount = zeros(nMacro, nMacro);
tProb = zeros(nMacro, nMacro);
[tCount, tProb] = transCount(trajMacro, traj_num, lagtime, nMacro); %using the msm population as a reference

terminate = 500;%modify
timestep = floor(terminate/lagtime);
realtProb = zeros(nMacro, nMacro, timestep+1);
realtCount = zeros(nMacro, nMacro, timestep+1);
tradtProb = zeros(nMacro, nMacro, timestep+1);

realtProb(:, :, 1)=eye(nMacro);
tradtProb(:, :, 1) = eye(nMacro);

[realtCount(:, :, 2) realtProb(:, :, 2)] = transCount_jump(trajMacro, traj_num, lagtime, nMacro, 1); %using sliding window

error = zeros(nMacro, nMacro, timestep+1);

tradtProb(:, :, 2) = realtProb(:, :, 2);
for time = 3:timestep+1
    tradtProb(:, :, time) = tradtProb(:, :, time-1)*realtProb(:, :, 2);
end

if strcmp(option, 'noe') %choose noe's equation to calculate, this equation is no longer used in pyemma
        %real 
        for time = 3:timestep+1
		time
                [realtCount(:, :, time), realtProb(:, :, time)] = transCount_jump(trajMacro, traj_num, lagtime*(time-1), nMacro, 1);%sliding window
        end

        %error bar
         %reference: ideal case if Markovian
         %CK test error bar, reference: An introduction to Markov State Models and Their Application to Long Timescale Molecular Simulation, Chapter 4.8
         %Formula: (4.52) in Chapter 4.8, page 57
         for time = 2:timestep+1
                for j = 1:nMacro
                        colsum = sum(realtCount(:, j, time));
                        for k = 1:nMacro
                                error(k, j, time) = sqrt((time-1)*(realtProb(k, j, time)-realtProb(k, j, time)^2)/colsum);
                        end
                end
         end
elseif strcmp(option, 'bootstrap') %by default use all the points by jumping window
        for time = 2:timestep+1
                AAA = residence_prob_adv(resultdir, trajMacro, traj_num, traj_len, nMacro, (time-1)*lagtime, 10, 'used_up'); %100 experiments
                realtProb(:, :, time) = reshape(mean(AAA), [nMacro, nMacro]);  %by default, sliding window
                error(:, :, time) = reshape(std(AAA), [nMacro nMacro]);
        end
elseif strcmp(option, 'once') %if specify
        for time = 2:timestep+1
                AAA = residence_prob(resultdir, trajMacro, traj_num, traj_len, nMacro, (time-1)*lagtime, 100, 'once'); %100 experiments
                realtProb(:, :, time) = reshape(mean(AAA), [nMacro, nMacro]);
                error(:, :, time) = reshape(std(AAA), [nMacro nMacro]);
        end
end

disp('begin plotting elements of macrotProb with reference');
x = [0:lagtime:lagtime*(timestep)]*timeunit;  %unit: micro-second

for j = 1:4
        yreal = reshape(realtProb(j, j,:),1,timestep+1);
        yref = reshape(tradtProb(j, j, :),1,timestep+1);
        error_bar = reshape(error(j, j, :), 1, timestep+1);
        errorbar(x, yreal, error_bar, 'o-');
        hold on;
        plot(x, yref, 'x-');
        hold off;
        legend('T(n\tau)', 'T(\tau)^n'');
        title(strcat('Chapamn-Kolmogromov test for element ',' T_{', num2str(j),num2str(j),'}'));
        axis([0 max(x) 0 1]);
        print(strcat(resultdir,'/','T',num2str(j),num2str(j)),'-dpng');
end

x = [0:lagtime:lagtime*(timestep)]*timeunit;%unit: micro-second,

%record new one
AA = zeros(length(x), 17);
AA(:, 1) = x';
number = 1;
for j=1:nMacro
	for k = 1:nMacro
        	temp=reshape(realtProb(j,k,:), 1, timestep+1);
        	AA(:, number+1) = temp';
		number = number+1;
	end
end
dlmwrite(strcat(resultdir, 'CK_real.txt'), AA, 'delimiter', '\t');

AA = zeros(length(x), 17);
AA(:, 1) = x';
number = 1;
for j=1:nMacro
	for k = 1:nMacro
        	temp=reshape(tradtProb(j,k,:), 1, timestep+1);
        	AA(:, number+1) = temp';
		number = number+1;
	end
end
dlmwrite(strcat(resultdir, 'CK_trad.txt'), AA, 'delimiter', '\t');

AA = zeros(length(x), 17);
AA(:, 1) = x';
number = 1;
for j=1:nMacro
	for k = 1:nMacro
        	temp=reshape(error(j,k,:), 1, timestep+1);
        	AA(:, number+1) = temp';
		number = number+1;
	end
end
dlmwrite(strcat(resultdir, 'CK_errorbar.txt'), AA, 'delimiter', '\t');
