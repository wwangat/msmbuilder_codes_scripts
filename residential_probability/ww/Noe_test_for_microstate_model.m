%
% noe test for a system under microstate Markovian lagtime
%using jumping window
option='bootstrap'; %we provide two options: bootstrap, noe's formula 
resultdir = 'noe_micro_bootstrap_sliding/';
lagtime=150;
if (exist(resultdir) == 0)
	mkdir(resultdir);
end

%addpath('/home/wang/Labwork/MOE_correlation/correlation_test/subroutine/');
trajlist = importdata('trajlist_micro.txt');
traj_len=[];
traj_num  = length(trajlist);

for j = 1:traj_num
    trajMicro{j} = importdata(trajlist{j})+1; %starting from 1 now
end

nMicro = 600;


%for microstates
tCount = zeros(nMicro, nMicro);
tProb = zeros(nMicro, nMicro);
[tCount, tProb] = transCount(trajMicro, traj_num, lagtime, nMicro); %using the msm population as a reference
%lagtime=1, to be consistent with siqin
colsum=sum(tCount, 1);
[value, index] = sort(colsum, 'descend');

Noeindex = index(1:20);
population = value/sum(value);
disp('the corresponding population of most popular 20 states are:');
[(Noeindex-1)' (population(1:20))']

terminate = 750;
timestep = floor(terminate/lagtime);
realtProb = zeros(nMicro, nMicro, timestep+1);
realtCount = zeros(nMicro, nMicro, timestep+1);
tradtProb = zeros(nMicro, nMicro, timestep+1);
realtProb(:, :, 1)=eye(nMicro);
tradtProb(:, :, 1) = eye(nMicro);

[realtCount(:, :, 2) realtProb(:, :, 2)] = transCount_jump(trajMicro, traj_num, lagtime, nMicro, 1); %using sliding window

error = zeros(nMicro, nMicro, timestep+1);

tradtProb(:, :, 2) = realtProb(:, :, 2);
for time = 3:timestep+1
    tradtProb(:, :, time) = tradtProb(:, :, time-1)*realtProb(:, :, 2);
end

if strcmp(option, 'noe') %choose noe's equation to calculate, this equation is no longer used in pyemma
        %real 
        for time = 3:timestep+1
                [realtCount(:, :, time), realtProb(:, :, time)] = transCount_jump(trajMicro, traj_num, lagtime*(time-1), nMicro, 1);%sliding window
        end

        %error bar
         %reference: ideal case if Markovian
         %Noe test error bar, reference: An introduction to Markov State Models and Their Application to Long Timescale Molecular Simulation, Chapter 4.8
         %Formula: (4.52) in Chapter 4.8, page 57
         for time = 2:timestep+1
                for j = 1:nMicro
                        colsum = sum(realtCount(:, j, time));
                        for k = 1:nMicro
                                error(k, j, time) = sqrt((time-1)*(realtProb(k, j, time)-realtProb(k, j, time)^2)/colsum);
                        end
                end
         end
elseif strcmp(option, 'bootstrap') %by default use all the points by jumping window
        for time = 2:timestep+1
                AAA = residence_prob_adv(resultdir, trajMicro, traj_num, traj_len, nMicro, (time-1)*lagtime, 100, 'used_up'); %100 experiments
                realtProb(:, :, time) = reshape(mean(AAA), [nMicro, nMicro]);  %by default, sliding window
                error(:, :, time) = reshape(std(AAA), [nMicro nMicro]);
        end
elseif strcmp(option, 'once') %if specify
        for time = 2:timestep+1
                AAA = residence_prob(resultdir, trajMicro, traj_num, traj_len, nMicro, (time-1)*lagtime, 100, 'once'); %100 experiments
                realtProb(:, :, time) = reshape(mean(AAA), [nMicro, nMicro]);
                error(:, :, time) = reshape(std(AAA), [nMicro nMicro]);
        end
end

disp('begin plotting elements of macrotProb with reference');
 x = [0:lagtime:lagtime*(timestep)]*0.1;  %unit: ns

for j = 1:20
        yreal = reshape(realtProb(Noeindex(j), Noeindex(j),:),1,timestep+1);
        yref = reshape(tradtProb(Noeindex(j), Noeindex(j), :),1,timestep+1);
        error_bar = reshape(error(Noeindex(j), Noeindex(j), :), 1, timestep+1);
        errorbar(x, yreal, error_bar, 'o-');
        hold on;
        plot(x, yref, 'x-');
        hold off;
        legend('T(n\tau)', 'T(\tau)^n');
        title(strcat('Chapamn-Kolmogromov test for element ',' T_{', num2str(Noeindex(j)),num2str(Noeindex(j)),'}'));
        axis([0 max(x) 0 1]);
        print(strcat(resultdir,'/','T',num2str(j),num2str(j)),'-dpng');
end

x = [0:lagtime:lagtime*(timestep)]*0.1;
%record new one
AA = zeros(length(x), 21);
AA(:, 1) = x';
for j=1:20
        temp=reshape(realtProb(Noeindex(j),Noeindex(j),:), 1, timestep+1);
        AA(:, j+1) = temp';
end
dlmwrite(strcat(resultdir, 'Noe_150_real.txt'), AA, 'delimiter', '\t');

AA = zeros(length(x), 21);
AA(:, 1) = x';
for j=1:20
        temp=reshape(tradtProb(Noeindex(j),Noeindex(j),:), 1, timestep+1);
        AA(:, j+1) = temp';
end
dlmwrite(strcat(resultdir, 'Noe_150_trad.txt'), AA, 'delimiter', '\t');
AA = zeros(length(x), 21);
AA(:, 1) = x'; 
for j = 1:20
    temp = reshape(error(Noeindex(j), Noeindex(j), :), 1, timestep+1);
    AA(:, j+1) = temp';
end
dlmwrite(strcat(resultdir, 'Noe_150_errorbar.txt'), AA, 'delimiter', '\t');
