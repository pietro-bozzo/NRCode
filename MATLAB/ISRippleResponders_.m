function [outputArg1,outputArg2] = ISRippleResponders_(session,regs,events,opt)
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
[rip.on,ind] = Restrict(rip.us,us_avals);
is_off = true(size(rip.us));
is_off(ind) = false;
rip.off = rip.us(is_off);
names = ["all","us","nus","on","off"];

% shift wrt sws, to allow shuffling baseline
for name = names
  rip.(name) = Restrict(rip.(name),sws_intervals,'shift','on');
end

% shuffle ripples
for name = "all"
  shuffle.("rip_"+name) = zeros(numel(rip.(name)),opt.shuffle);
  for k = 1 : opt.shuffle
    shuffle.("rip_"+name)(:,k) = shuffleSpikes(rip.(name));
  end
end

% region data
spikes = R.spikes('all',regs);
spikes = Restrict(spikes,sws_intervals,'shift','on');
reg_units = R.units(regs);
window = [-0.2,0.2];
limits = window + [-1.5,0];
smooth = 2;
n_bins = 251;
baseline = -1;

% declare variables
for name = names
  peth.(name) = nan(numel(reg_units),n_bins);
  response.(name) = nan(numel(reg_units),1);
  delay.(name) = nan(numel(reg_units),1);
end

% baseline: response to shuffled ripples (random events)
for k = 1 : opt.shuffle
  [~,t,P] = PETH(spikes(:,1),shuffle.rip_all(:,k),'durations',limits,'smooth',smooth,'nBins',n_bins,'group',spikes(:,2));
  shuffle.response(:,k) = cellfun(@(x) absResponse(x,t,baseline,window),P);
end

% PETH over all spikes and all conditions REMEMBER TO CHECK N SPIKES
for name = names
  [~,t,P] = PETH(spikes(:,1),rip.(name),'durations',limits,'smooth',smooth,'nBins',n_bins,'group',spikes(:,2));
  peth.(name) = vertcat(P{:});
  [response.(name),delay.(name)] = cellfun(@(x) absResponse(x,t,baseline,window),P);
  p.(name) = percentRank(shuffle.response.',response.(name).','center').';
end

% find responders
h.cat = [];
for name = names
  h.(name) = holmBonferroni(p.(name));
  h.cat = [h.cat,h.(name)];
  responders.(name) = reg_units(logical(h.(name)));
end

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

%% plot PETH
i = 1;
ind_to_plot = [1];
makeFigure('peth',"PETH unit "+i);
% compute shuffled peth
spikes_unit = spikes(spikes(:,2)==reg_units(i),1);
l = 1;
for k = 1 : 10 : opt.shuffle
  [~,~,d(:,l)] = PETH(spikes_unit,shuffle.rip_all(:,k),'durations',[limits(1)-diff(limits)*1.5,limits(2)],'smooth',smooth,'nBins',n_bins); l = l + 1;
end
% plot
for k = 1 : l-1
  plot(t,d(:,k),'Color',0.5*k/(l-1)*[1,1,1]+0.5)
end
for n = ind_to_plot
  hand(n) = plot(t,peth.(names(n))(i,:),'Color',myColors(n),'DisplayName',names(n)+", r: "+response.(names(n))(i)+", "+h.(names(n))(i));
end
hand(end+1) = plot(nan,'Color',[0.7,0.7,0.7],'DisplayName',"shuffle, \langler\rangle: "+mean(shuffle.response(i,:)));
legend(hand)
xlabel('time from ripple (s)'), ylabel('spike rate (Hz)')

makeInset(0.1,0.75,0.12,NaN);
plotDistr(shuffle.response(i,:),'colors',.7*[1,1,1],'nbins',20)
for n = ind_to_plot
  xline(response.(names(n))(i),'Color',myColors(n),'LineWidth',1.3)
end
clear hand d n k l ind_to_plot

%% plot bootstrap test
makeFigure('peth',"Bootsrap unit "+i+", responder: "+h(i)+", r: "+response(i));

end