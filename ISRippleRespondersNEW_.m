function [outputArg1,outputArg2] = ISRippleRespondersNEW_(session,regs,events,opt)
% I APPARENTLY LOST RESPONDERS WHEN LIMITING TO SLEEPM, MAYBE DIFFERENCE WRT SLEEPN?

arguments
  session (1,:) char
  regs (:,1) string
  events (:,1) string = []
  opt.shuffle (1,1) {mustBeNumeric,mustBeInteger,mustBeNonnegative} = 500
  opt.save (1,1) {mustBeLogical} = false
  opt.load (1,1) {mustBeLogical} = true
  opt.verbose (1,1) {mustBeLogical} = false
end

R = regions(session,'states','sws','events',["ripples","InfraSlowRhythm/slownr","InfraSlowRhythm/slowavalnr"],'verbose',false);
sws_intervals = R.eventIntervals(events,'sws');
us_intervals = R.eventIntervals(events,'sws','slownr');
us_avals = R.eventIntervals(events,'sws','slowavalnr');

% define PETH function
limits = [-0.5,0.5];
smooth = 4;
n_bins = 251;
pethf = @(x,y,z) PETH(x,y,'group',z,'durations',limits,'smooth',smooth,'nBins',n_bins,'fast','on');

% keep only regions found in data
[~,ind] = ismember(R.ids,regs);
ind = ind(ind~=0);
regs = regs(ind);

% ripples
type = "ripples";
rip.all = R.eventInfo(type);
rip.all = rip.all(:,1);

% restrict in ISR, ON and OFF
[rip.us,ind] = Restrict(rip.all,us_intervals);
is_nus = true(size(rip.all));
is_nus(ind) = false;
rip.nus = rip.all(is_nus);
%[rip.on,ind] = Restrict(rip.us,us_avals);
%is_off = true(size(rip.us));
%is_off(ind) = false;
%rip.off = rip.us(is_off);
names = ["all","us","nus"]; % ,"on","off"];

% shift wrt sws, to allow shuffling baseline
rip = structfun(@(x) Restrict(x,sws_intervals,'shift','on'),rip,'UniformOutput',false);

% shuffle ripples
shuffle.rip = zeros(numel(rip.all),opt.shuffle);
for k = 1 : opt.shuffle
  shuffle.rip(:,k) = shuffleSpikes(rip.all);
end

% subsample ripples
n_rip = floor(min(structfun(@numel,rip))/100) * 100;
rip_subs = structfun(@(x) sort(x(randperm(numel(x),n_rip))),rip,'UniformOutput',false);
shuffle.rip_subs = zeros(n_rip,opt.shuffle);
for k = 1 : opt.shuffle
  shuffle.rip_subs(:,k) = sort(shuffle.rip(randperm(size(shuffle.rip,1),n_rip),k));
end

% region data
spikes = R.spikes('all',regs);
spikes = Restrict(spikes,sws_intervals,'shift','on');
reg_units = R.units(regs);
%t = linspace(limits(1),limits(2),n_bins+1);
%phases = linspace(0,2*pi,n_bins+1);
%phases = (phases(1:end-1) + phases(2:end)) / 2;

% declare variables
for name = names
  peth.(name) = nan(numel(reg_units),n_bins);
  [max_response.(name),min_response.(name),delay.(name)] = deal(nan(numel(reg_units),1));
  p_value.(name) = nan(numel(reg_units),2);
end

% BIG sync MAYBE NO TIME IS GAINED DUE TO SIZE OF ARRAYS, MUST CHECK
% rip_cat = [shuffle.rip_subs(:); structCat(rip_subs,"vert")];
% rip_id = [repelem((1:opt.shuffle).',size(shuffle.rip_subs,1),1); repelem(opt.shuffle+(1:numel(fieldnames(rip_subs))).',structfun(@numel,rip_subs),1)];
% [sync,Ie,Is] = Sync(spikes(:,1),rip_cat,'durations',limits);
% Is = spikes(Is,2);
% Ie = rip_id(Ie);
% d = discretize(sync,t);
% d = phases(d).'; % d(i) is phase of sync(i)
% 
% for n = 1 : numel(names)
% 
%   % assess significance
%   for i = 1 : numel(reg_units)
%     % analyze i-th unit
%     p_shuffled = cellfun(@(x) x(x(:,2)==reg_units(i),1),shuffle.phases,'UniformOutput',false);
%     len = max(cellfun(@numel, p_shuffled));
%     p_sh = nan(len,opt.shuffle);
%     for k = 1 : numel(p_shuffled)
%       p_sh(1:numel(p_shuffled{k}),k) = p_shuffled{k};
%     end
% 
%     p_sh = 
%     stat = circEntrainment(d(Is==reg_units(i) & Ie==500+n),'mode','phase','shuffle',p_sh);
% 
% 
%     p_value.(names(n))(i,1) = stat.p;
%   end
% 
% end

% baseline: response to shuffled ripples (random events) I SHOULD SELECT A HOMOGENEEOUS N OF RIP IN ALL CONDITIONS I WANT TO COMPARE (+ SHUF)
% COULD ALSO SEE IF RIP AT TRANSITIONS RECRUTE MORE RESPONDERS? MAYBE DUH DUE TO COUPLING
for k = 1 : opt.shuffle
  
  [~,~,P] = pethf(spikes(:,1),shuffle.rip_subs(:,k),spikes(:,2));
  shuffle.max_response(:,k) = cellfun(@(x) max(x),P);
  shuffle.min_response(:,k) = cellfun(@(x) min(x),P);

  % WAS USING R, BUT THIS MISSES CASES WHERE FR IS LOWER IN WHOLE WINDOW
  % [sync,~,Is] = Sync(spikes(:,1),shuffle.rip_subs(:,k),'durations',limits,'fast','on'); % I COULD sync TO ALL SORTED RIPS AND REMEMBER THEIR TYPE!! OR AT LEAST TO ALL SHUFFLED RIPS
  % d = discretize(sync,t);
  % p = phases(d).'; % p(i) is phase of sync(i)
  % groups = spikes(Is,2); % groups(i) is unit which produced p(i)
  % shuffle.phases{k,1} = [p,groups];

end

% PETH over all spikes and all conditions REMEMBER TO CHECK N SPIKES
for name = names

  % PETH
  [~,t,P] = pethf(spikes(:,1),rip_subs.(name),spikes(:,2));
  peth.(name) = vertcat(P{:});
  max_response.(name) = cellfun(@(x) max(x),P);
  min_response.(name) = cellfun(@(x) min(x),P);

  p_value.(name)(:,1) = percentRank(shuffle.max_response.',max_response.(name).','center').';
  p_value.(name)(:,2) = percentRank(shuffle.min_response.',min_response.(name).','center').';

  % [sync,~,Is] = Sync(spikes(:,1),rip_subs.(name),'durations',limits,'fast','on');
  % d = discretize(sync,t);
  % p = phases(d).'; % p(i) is phase of sync(i)
  % groups = spikes(Is,2); % groups(i) is unit which produced p(i)

  % assess significance
  % for i = 1 : numel(reg_units)
    % analyze i-th unit
    % p_shuffled = cellfun(@(x) x(x(:,2)==reg_units(i),1),shuffle.phases,'UniformOutput',false);
    % len = max(cellfun(@numel, p_shuffled));
    % p_sh = nan(len,opt.shuffle);
    % for k = 1 : numel(p_shuffled)
    %   p_sh(1:numel(p_shuffled{k}),k) = p_shuffled{k};
    % end
    % stat = circEntrainment(p(groups==reg_units(i)),'mode','phase','shuffle',p_sh);
    % p_value.(name)(i,1) = stat.p;
  % end

end

% find responders
for name = names
  h_max.(name) = holmBonferroni(p_value.(name)(:,1));
  h_min.(name) = holmBonferroni(p_value.(name)(:,2));
  %responders.(name) = reg_units(logical(h.(name)));
end
h_max.cat = structCat(h_max,"horz");
h_min.cat = structCat(h_min,"horz");
h = h_max.cat | h_min.cat;

%    READ LITERATURE BEFORE DOING STUFF
% CHECK IF RESPONDERS CHANGE IN OFF, ON
% RESPONDERS HAVE weaker response in ISR compared with outside! but response in off (and not on) is as high as outside
% what about response delay? for all same groups
% Can check plots of peth or check whether number of responders changes
% can also check if taking only ripples having ON-OFF transition close changes


figure, distPlot(h.cat,[],'withinlines',1);
figure, distPlot([response.all,response.us,response.nus,response.on,response.off],[],'withinlines',1);
dataa = [];
for name = names
  dataa = [dataa;response.(name)];
end
ANOVATests(dataa,'parametric',false);

end

% --- helper functions ---

function [response,delay] = absResponse(peth,t,baseline,window)
  
  % mean over baseline
  ind = find(t<baseline,1,'last');
  m = mean(peth(1:ind)); 
  s = std(peth(1:ind));
  % response in window
  is_window = t >= window(1) & t <= window(2);
  t_window = t(is_window);
  peth_window = peth(is_window);
  [~,delay] = max(abs(peth_window-m));
  response = (peth_window(delay) - m) / s; % standardize w.r.t. baseline
  delay = t_window(delay); % convert to s

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  %% PETH colormap
  i = 1;
  name = 'all';
  makeFigure('colormap',"PETH unit "+i+", "+name);
  spikes_unit = spikes(spikes(:,2)==reg_units(i),1);
  P = pethf(spikes_unit,rip_subs.(name),[]);
  PlotColorMap(P,'x',t)
  clearvars i name spikes_unit P

%% plot PETH
i = 6;
ind_to_plot = [1,2,3];
makeFigure('peth',"PETH unit "+i);
spikes_unit = spikes(spikes(:,2)==reg_units(i),1);
% plot shuffled peth
for k = 1 : 10 : opt.shuffle
  [~,~,P] = pethf(spikes_unit,shuffle.rip_subs(:,k),[]);
  plot(t,P,'Color',0.5*k/(opt.shuffle-1)*[1,1,1]+0.5)
end
% plot
for n = ind_to_plot
  [~,~,P] = pethf(spikes_unit,rip_subs.(names(n)),[]);
  hand(n) = plot(t,P,'Color',myColors(n),'DisplayName',names(n)+", r: "+max_response.(names(n))(i)+", "+h_max.(names(n))(i));
end
%hand(end+1) = plot(nan,'Color',[0.7,0.7,0.7],'DisplayName',"shuffle, \langler\rangle: "+mean(shuffle.response(i,:)));
legend(hand)
xlabel('time from ripple (s)'), ylabel('spike rate (Hz)')

% makeInset(0.1,0.75,0.12,NaN);
% plotDistr(shuffle.response(i,:),'colors',.7*[1,1,1],'nbins',20)
% for n = ind_to_plot
%   xline(response.(names(n))(i),'Color',myColors(n),'LineWidth',1.3)
% end
clear hand P n k ind_to_plot

%% plot bootstrap test
makeFigure('peth',"Bootsrap unit "+i+", responder: "+h(i)+", r: "+response(i));

end