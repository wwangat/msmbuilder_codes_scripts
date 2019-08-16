%transition count matrix
%
% ============================================================================
%       Filename:  transCount_jump.m
%    Description:  input a trajectory containing many segment, output transition probability matrix(column normalized)
%                   calculate by jumping window
%       Modified:  2015-12-17 01:04
%          Usage:  [tCount, tProb] = transCount_jump(traj, 100, 10000, 2, 4, jump_step), traj_len is a array, jump_step is to guarantee that counts are independent
%         Author:  WANG Wei        (wwangat@gmail.com)
% ============================================================================
%

function  [tCount, tProb] = transCount_jump(traj, traj_num, lagtime, nStates, jump_step)
  tCount = zeros(nStates, nStates);
  for j = 1:traj_num
    trajectory = traj{j};
    traj_len = length(trajectory);
    for k = 1:jump_step:traj_len-lagtime
      old_index = trajectory(k); %since this time start from 1
      new_index = trajectory(k+lagtime);
	if (old_index>=1 & new_index>=1)
      	    tCount(new_index, old_index) = tCount(new_index, old_index)+1;
	end
    end
  end
  %%definition of Tij: transition from state j to state i
  tCount = (tCount+tCount')/2.0;
  tProb = zeros(nStates, nStates);
  for j = 1:nStates
    colsum = sum(tCount(:, j));
    if colsum~=0
        tProb(:, j) = tCount(:, j)/colsum;
    end
  end
end
