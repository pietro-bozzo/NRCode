function ISRPopulGLM_(session,regs,events) % WAS  [r_squared,max_delay] = 

arguments
  session (1,:) char
  regs (:,1) string
  events (:,1) string = []
end

R = regions(session,'states','sws','phases',events,'events',["InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"],'verbose',false);

% declare variables
% n_delays = 10;
% for r = regs'
%   r_squared.(r) = nan(1,2*n_delays+1);
%   max_delay.(r) = NaN;
% end

% keep only regions found in data
[~,ind] = ismember(R.ids,regs);
ind = ind(ind~=0);
regs = regs(ind);

% y: NR population firing
window = 0.05; % s
[y,time.nr] = R.firingRate('sws','nr','window',window,'step',4);

% save data for Python analysis
for j = 1 : numel(regs)
  % X: units firing rates
  [data.(regs(j)),time.(regs(j))] = R.unitFiringRates('sws',regs(j),'window',window,'step',4);
end
data.y = y;
% make length uniform
min_len = min(structfun(@(x) size(x,1),data));
data = structfun(@(x) x(1:min_len,:),data,'UniformOutput',false);
save4Python(session,'populGLM',data,'PopulGLM')

% --- tryed to implement in Python ---
% % fit GLM
% window = 0.2; % s
% r_squared.delays = (-n_delays : n_delays) * (time.nr(2)-time.nr(1));
% downs_factor = 10;
% r = zeros(5,2*n_delays+1);
% for j = 1 : numel(regs)
% 
%   % X: units firing rates
%   [X,time.(regs(j))] = R.unitFiringRates('all',regs(j),'window',window,'step',4,'smooth',5);
%   % shorten and downsample response
%   y_downs = decimate(y(n_delays+1:size(X,1)-n_delays),downs_factor);
%   time.nr_downs = time.nr(n_delays+1:size(X,1)-n_delays);
%   time.nr_downs = time.nr_downs(1:downs_factor:end);
%   % divide data in 5 folds
%   ids = repelem((1:5).',ceil(numel(y_downs)/5),1);
%   ids = ids(1:numel(y_downs));
% 
%   for delay = -n_delays : n_delays
% 
%     % apply time delay
%     X_delayed = X(n_delays+1+delay:end+delay-n_delays,:);
%     % downsample predictors
%     X_downs = [];
%     for k = 1 : size(X_delayed,2)
%       X_downs(:,k) = decimate(X_delayed(:,k),downs_factor);
%     end
%     time.(regs(j)+"_downs") = time.(regs(j))(n_delays+1+delay:end+delay-n_delays);
%     time.(regs(j)+"_downs") = time.(regs(j)+"_downs")(1:downs_factor:end);
%     % standardize
%     X_downs = zscore(X_downs,0,1);
% 
%     for i = 1 : 5
% 
%       % fit on this fold
%       is_train = ids ~= i;
%       mdl = fitglm(X_downs(is_train,:),y_downs(is_train),'interactions','Distribution','gamma','Link','log'); % can also try including interactions
%       prediction = predict(mdl,X_downs(~is_train,:));
%       % evaluate error
%       r(i,n_delays+1+delay) = 1 - sum((y_downs(~is_train)-prediction).^2) / sum((y_downs(~is_train)-mean(y_downs(is_train))).^2);
% 
%     end
%   end
% 
%   r_squared.(regs(j)) = mean(r,1);
%   [~,max_delay.(regs(j))] = max(r_squared.(regs(j)));
%   max_delay.(regs(j)) = r_squared.delays(max_delay.(regs(j)));
% 
% end

end

%% ---

function plotInDebug

  figure; hold on
  plot(time.nr_downs,y_downs,'DisplayName','data')
  plot(time.nr_downs(~is_train),prediction,'DisplayName','prediction')
  legend

  %%

  figure; hold on
  fields = fieldnames(r_squared);
  for f = 1 : numel(fields)-1
    hand(f) = plot(r_squared.delays,r_squared.(fields{f}),'DisplayName',fields{f},'Color',myColors(f));
    xline(max_delay.(fields{f}),'Color',myColors(f))
  end
  legend(hand), xlabel('delay (s)'), ylabel('R^2')
  clearvars fields f hand

end