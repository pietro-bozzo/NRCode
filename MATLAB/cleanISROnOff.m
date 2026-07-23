function [on,off,intervals,on_off,off_on] = cleanISROnOff(aval,intervals)

if isempty(aval) || isempty(intervals)
  [on,off,intervals,on_off,off_on] = deal([]);
  return
end

% find avalanches which start and stop ISR intervals
[~,start_ind] = Restrict(aval(:,1)-1e-5,intervals);
start_ind = setdiff((1:size(aval,1)).',start_ind);
stop_ind = [start_ind(2:end)-1;size(aval,1)];

indeces = [start_ind,stop_ind];
is_ok = true(size(aval,1),1); % is_ok(i) is true if i-th avalanche must be kept in on
indeces(indeces(:,1)==indeces(:,2),:) = []; % keep intervals containing only one avalanche THIS SHOULD NOT BE POSSIBLE AND IS NOT GUARANTEED TO HAVE AN OFF

% set to 0 one element of is_ok for every row of indeces, randomly excluding either the first or last avalanche of every ISR interval
random_ind = randi(2,size(indeces,1),1);
is_ok(indeces(sub2ind(size(indeces),(1:size(indeces,1)).',random_ind))) = false;

% build ON and OFF intervals
on = aval(is_ok,:);
off = SubtractIntervals(intervals,aval);

% rebuild intervals without excluded avalanches
intervals(random_ind==1,1) = aval(start_ind(random_ind==1),2);
intervals(random_ind==2,2) = aval(stop_ind(random_ind==2),1);

% valid transitions
on_off = Restrict(on(:,2)+1e-7,off);
off_on = Restrict(on(:,1)-1e-7,off);