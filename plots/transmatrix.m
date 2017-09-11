%
% ============================================================================
%       Filename:  transmatrix.m
%    Description:  input a trajectory containing many segment, output transition probability matrix(column normalized)
%       Modified:  2015-12-17 01:04
%          Usage:  tProb = transmatrix(traj, 100, 10000, 2, 4), traj_len is a array
%         Author:  WANG Wei        (wwangat@gmail.com)
% ============================================================================
%

function  [tProb] = transmatrix(traj, traj_num, lagtime, nStates)
  tCount = zeros(nStates, nStates);
  for j = 1:traj_num
    trajectory = traj{j};
    traj_len = length(trajectory);
    for k = 1:traj_len-lagtime
      old_index = trajectory(k); %since this time start from 1
      new_index = trajectory(k+lagtime);
      tCount(new_index, old_index) = tCount(new_index, old_index)+1;
    end
  end
  %%definition of Tij: transition from state j to state i
  tCount = (tCount+tCount')/2.0;
  if (length(find(sum(tCount, 1)==0)))~=0
     disp('need to remove disconnected states first, calculating by ignoring disconnected states');
  end
  tProb = zeros(nStates, nStates);
  for j = 1:nStates
    colsum = sum(tCount(:, j));
    if (colsum ~= 0)
        tProb(:, j) = tCount(:, j)/colsum;
    end
  end
end
