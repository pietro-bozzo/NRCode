%% Analyse entrainment of spikes and ripples by US rhythm

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/IS_states.batch';
saveFlag = false;

%% 1. ripple-spindle coupling

% load data
events = ["ripples","spindles"];
window = 1; % s
task = "sleepm";
duration = 1200;
args = {window,'event',task,'state','sws','duration',duration,'shuffle',1000};
[phase,is_coupl,not_us,us_phase,phase_sh,is_coupl_sh,session] = runBatch(batch_file,@ISRipSpinCoupl_,args);
clearvars args

%%
[~,~,session] = unique(vertcat(session{:}),'stable');
animal = reverseCellStruct(phase,@numel);
animal = structfun(@(x) repelem(session,x,1),animal,'UniformOutput',false);
animal.session = session;
phase = reverseCellStruct(phase);
is_coupl = reverseCellStruct(is_coupl);
not_us = reverseCellStruct(not_us);
phase_sh = reverseCellStruct(phase_sh);
is_coupl_sh = reverseCellStruct(is_coupl_sh);
% reference distribution of ultra-slow phases, smoothing away border effects
[phase_ref,phi] = CircularDistribution(vertcat(us_phase{:}),'nBins',250,'normalize','pdf');
phase_ref(round(numel(phase_ref)*5/12):round(numel(phase_ref)*7/12)) = circularSmooth(phase_ref(round(numel(phase_ref)*5/12):round(numel(phase_ref)*7/12)),'gaussian',15);
phase_ref([round(numel(phase_ref)*11/12):end,1:round(numel(phase_ref)/12)]) = circularSmooth(phase_ref([round(numel(phase_ref)*11/12):end,1:round(numel(phase_ref)/12)]),'gaussian',15);
clearvars session

%% analyze entrainment
for event = events

  [stats.(event),distr.(event),~,shuffle.(event)] = circEntrainment(phase.(event),phase_ref,'mode','phase','shuffle',phase_sh.(event));
  % repeat separating data in coupled / not coupled
  [stats_coupl.(event),distr_coupl.(event),~,shuffle_coupl.(event)] = circEntrainment(phase.(event)(is_coupl.(event)),phase_ref,'mode','phase','shuffle',phase_sh.(event)(is_coupl.(event),:));
  [stats_uncoupl.(event),distr_uncoupl.(event),~,shuffle_uncoupl.(event)] = circEntrainment(phase.(event)(~is_coupl.(event)),phase_ref,'mode','phase','shuffle',phase_sh.(event)(~is_coupl.(event),:));

  % compute separately per animal
  for a = 1 : 2
    a_ind = animal.(event) == a;
    [animal.stats.(event){a},animal.distr.(event){a},~,animal.shuffle.(event){a}] = circEntrainment(phase.(event)(a_ind),phase_ref,'mode','phase','shuffle',phase_sh.(event)(a_ind,:));

    final_ind = a_ind & is_coupl.(event);
    [animal.coupl.stats.(event){a},animal.coupl.distr.(event){a},~,animal.coupl.shuffle.(event){a}] = circEntrainment(phase.(event)(final_ind),phase_ref,'mode','phase', ...
      'shuffle',phase_sh.(event)(final_ind,:));
    final_ind = a_ind & ~is_coupl.(event);
    [animal.uncoupl.stats.(event){a},animal.uncoupl.distr.(event){a},~,animal.uncoupl.shuffle.(event){a}] = circEntrainment(phase.(event)(final_ind),phase_ref,'mode','phase', ...
      'shuffle',phase_sh.(event)(final_ind,:));
  end

end
clearvars shuffled event a_ind final_ind a

%% plot on all data OLD POLAR CODE
% [~,axs] = makeFigure('all','Average entrainment, '+task,[1,2],'polar',true,'size',[1500,700]);
% OnOffAxes(0.4,axs,'ticks',[0.1,0.3])
% for i = 1 : numel(events)
%   polarSemplot(phi,shuffle.(events(i)).distr.',[0.7,0.7,0.7],'ax',axs(i),'legend','shuffled (mean ± sem)');
%   polarplot(axs(i),phi,circularSmooth(distr.(events(i)),[],'gaussian',10),'Color',paperColors(9+i),'LineWidth',1.3,'DisplayName','distribution');
%   polarplot(axs(i),[0,stats.(events(i)).phi],[0,stats.(events(i)).R],'k','DisplayName','center','LineWidth',2);
%   legend(axs(i))
%   center_status = (stats.(events(i)).phi+pi/2) / pi;
%   title(axs(i),events(i)+", R: "+num2str(stats.(events(i)).R,2)+", s: "+num2str(center_status,2)+" (center), p: "+num2str(stats.(events(i)).p),'FontSize',10)
% end
% 
% clear i axs center_status
%saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/USRipplesSpindles/1rip_spin_'+task,["png","svg"],'pause',1);

%% a) entrainment
[~,axs] = makeFigure('all','Pooled distribution, '+task,[1,2],'size',[1000,500]);
OnOffAxes(0.35,axs)
phase_ind = [ceil(250*3/4):250,1:floor(250*3/4)];
for i = 1 : numel(events)

  for a = 1 : 2
    plot(axs(i),[phi;phi+phi(end)],repmat(mean(animal.shuffle.(events(i)){a}.distr(phase_ind,:),2),2,1),'Color',sqrt(0.7)*[1,1,1],'LineWidth',1.3);
    plot(axs(i),[phi;phi+phi(end)],repmat(circularSmooth(animal.distr.(events(i)){a}(phase_ind),[],'gaussian',10),2,1),'Color',sqrt(paperColors(9+i)),'LineWidth',1.3);
  end

  hand(1) = plot(axs(i),[phi;phi+phi(end)],repmat(mean(shuffle.(events(i)).distr(phase_ind,:),2),2,1),'Color',0.7*[1,1,1],'LineWidth',1.3,'DisplayName','shuffled');
  hand(2) = plot(axs(i),[phi;phi+phi(end)],repmat(circularSmooth(distr.(events(i))(phase_ind),[],'gaussian',10),2,1),'Color',paperColors(9+i),'LineWidth',1.3,'DisplayName',"distribution");
  legend(axs(i),hand), title(axs(i),events(i)+", p: "+stats.(events(i)).p)

end
ylabel(axs(1),'\langleevent distribution\rangle'), yticks(axs,[0.1,0.2,0.3])

clear i axs phase_ind hand a
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3a_dist_'+task,["png","svg"]);

%% compute PETHs
[isr_peth,trans,t] = runBatch(batch_file,@RipSpinPETH_,{'limits',[-5,5;-5,5]});
isr_peth = reverseCellStruct(isr_peth);
trans = reverseCellStruct(trans);
t = t{1};

%% b) transition PETHs
[~,axs] = makeFigure('peth','Ripples, spindles transition PETHs',[1,4],'size',[1400,400]);
color_ind = [10,10,11,11];
fields = string(fieldnames(trans))';

for i = 1 : 4

  for a = 1 : 2
    a_ind = animal.session == a;
    plot(axs(i),t.trans,mean(trans.(fields(i))(a_ind,:),1,'omitmissing'),'Color',sqrt(paperColors(color_ind(i))));
  end
  plot(axs(i),t.trans,mean(trans.(fields(i)),1,'omitmissing'),'Color',paperColors(color_ind(i)));
  xline(axs(i),0,'--')
  
end
title(axs([1,3]),'ON ➞ OFF'), title(axs([2,4]),'OFF ➞ ON')
xlabel(axs,'time from transition (s)'), ylabel(axs(1),'\langleevent rate\rangle (Hz)'), set(axs,'XTick',[-5,0,5],'XLim',[-5,5],'YLim',[0.05,0.85])
clearvars fields color_ind axs a_ind a i
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3b_trans',["png","svg"]);

%% c) ripples vs spindles PETH
[~,axs] = makeFigure('rp_peth','Ripples vs spindles PETHs',[1,3],'size',[1200,400]);
fields = string(fieldnames(isr_peth))';
titles = ["all","ISR","nISR"];

j = 1;
for i = [1,3,4]

  for a = 1 : 2
    plot(axs(j),t.ref,mean(isr_peth.(fields(i))(animal.session == a,:),1,'omitmissing'),'Color',sqrt(paperColors(10)));
  end
  plot(axs(j),t.ref,mean(isr_peth.(fields(i)),1,'omitmissing'),'Color',paperColors(10));
  xline(axs(j),[-1,0],'--')
  title(axs(j),titles(j))

  j = j + 1;

end
xlabel(axs,'time from spindle (s)'), ylabel(axs(1),'\langleripple rate\rangle (Hz)'), set(axs,'XTick',[-5,-1,0,5],'YLim',[0.3,0.9])
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3c_rip_spin',["png","svg"]);
clearvars fields titles i j a axs

%% 
[~,axs] = makeFigure('rp_peth','Ripples vs spindles PETHs',[2,4],'size',[1600,800]);
color_ind = [10,11,10*ones(1,6)];
fields = string(fieldnames(isr_peth))';
titles = ["all","all (reverse)","ISR","nISR","ON","OFF","[ON/2,OFF/2]","[OFF/2,ON/2]"];

for i = 1 : 8

  for a = 1 : 2
    plot(axs(i),t.ref,mean(isr_peth.(fields(i))(animal.session == a,:),1,'omitmissing'),'Color',sqrt(paperColors(color_ind(i))));
  end
  plot(axs(i),t.ref,mean(isr_peth.(fields(i)),1,'omitmissing'),'Color',paperColors(color_ind(i)));
  xline(axs(i),0,'--')
  title(axs(i),titles(i))

end
xlabel(axs([1,3:8]),'time from spindle (s)'), xlabel(axs(2),'time from ripple (s)'), ylabel(axs([1,5]),'\langleevent rate\rangle (Hz)'), set(axs([1,3:8]),'YLim',[0.3,0.9]), set(axs(2),'YLim',[0.11,0.25])
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3c_rip_spin_suppl',["png","svg"]);
clearvars axs color_ind fields titles i a

%% d) coupling probability
makeFigure('perc','Coupled events per animal','size',[4,3]);
for a = 1 : 2
  for i = 1 : 2
    n = sum(not_us.(events(i))(animal.session==a,:),1);
    frac_animal((1:2)+2*(i-1),a) = [mean(is_coupl.(events(i))(animal.(events(i))==a)), n(2)/n(1)];
  end
end
hand = distPlot(frac_animal.'*100,'group2',[1,2],'smedian',false,'withinlines',true,'color',[paperColors(10);.7,.7,.7;paperColors(11);.7,.7,.7],'label',["US","nUS","US","nUS"],'ssize',55);
delete(hand.wl(2))
ylabel('percentage (%)'), yticks([15,30,45]), ylim([10,50])
%plotPercent(frac_animal*100,'colors',repelem([paperColors(10);.7,.7,.7;paperColors(11);.7,.7,.7],2,1) .* repmat([.9;1.1]*[1,1,1],4,1),'labels',["US","nUS","US","nUS"])
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3d_perc',["png","svg"]);
clearvars a i n frac_animal hand

%% e) percentage of coupled events as a f of status
bin_width = pi / 15;
phi_x = bin_width/2 : bin_width : 4*pi;
[~,axs] = makeFigure('prob','P of event being coupled, '+task,[1,2],'size',[7,3],'format','poster');
yt = [0,0.55;0.2,0.75];
for i = 1 : numel(events)

  OnOffAxes(yt(i,:),axs(i))

  % shuffle
  prob = zeros(size(phase_sh.(events(i)),2),numel(phi_x)/2);
  for j = 1 : size(phase_sh.(events(i)),2)
    bin_ind = ceil(phase_sh.(events(i))(:,j)/bin_width);
    bin_ind(bin_ind==0) = 1;
    prob(j,:) = accumarray(bin_ind,is_coupl_sh.(events(i))(:,j)) ./ accumarray(bin_ind,1);
  end
  semplot(phi_x,[prob,prob],[.7,.7,.7],'ax',axs(i))

  % animal curves
  for a = unique(animal.session).'
    a_ind = animal.(events(i)) == a;
    prob = accumarray(ceil(phase.(events(i))(a_ind)/bin_width),is_coupl.(events(i))(a_ind)) ./ accumarray(ceil(phase.(events(i))(a_ind)/bin_width),1);
    prob = circularSmooth(prob,'gaussian',3);
    plot(axs(i),phi_x,[prob;prob],'Color',sqrt(paperColors(9+i)))
  end

  % mean
  prob = accumarray(ceil(phase.(events(i))/bin_width),is_coupl.(events(i))) ./ accumarray(ceil(phase.(events(i))/bin_width),1);
  prob = circularSmooth(prob,'gaussian',3);
  plot(axs(i),phi_x,[prob;prob],'Color',paperColors(9+i))

end
ylabel(axs(1),'coupling probability'), yticks(axs,[0.1,0.3,0.5])
clear bin_width prob i axs yt phi_x a a_ind
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3e_prob_phi_'+task,["png","svg"]);
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/PosterFigures/prob_phi_'+task,"svg");

% TO DO: shuffle ripples or spindles preserving their US entrainment, and show coupling disappears: it is not just distributions wrt US that give it

%% plot coupled, uncoupled
[~,axs] = makeFigure('coupl',"Average ripple-spindle couples entrainment, "+task+", w: "+num2str(window,3)+" s",[1,4],'size',[1400,400]);
OnOffAxes(0.4,axs)
phase_ind = [ceil(250*3/4):250,1:floor(250*3/4)];
for i = 1 : numel(events)

  for a = 1 : 2
    plot(axs(2*i-1),[phi;phi+phi(end)],repmat(mean(animal.coupl.shuffle.(events(i)){a}.distr(phase_ind,:),2),2,1),'Color',sqrt(0.7)*[1,1,1],'LineWidth',1.3);
    plot(axs(2*i-1),[phi;phi+phi(end)],repmat(circularSmooth(animal.coupl.distr.(events(i)){a}(phase_ind),[],'gaussian',10),2,1),'Color',sqrt(paperColors(9+i)),'LineWidth',1.3);

    plot(axs(2*i),[phi;phi+phi(end)],repmat(mean(animal.uncoupl.shuffle.(events(i)){a}.distr(phase_ind,:),2),2,1),'Color',sqrt(0.7)*[1,1,1],'LineWidth',1.3);
    plot(axs(2*i),[phi;phi+phi(end)],repmat(circularSmooth(animal.uncoupl.distr.(events(i)){a}(phase_ind),[],'gaussian',10),2,1),'Color',sqrt(paperColors(9+i)),'LineWidth',1.3);
  end

  hand(1) = plot(axs(2*i-1),[phi;phi+phi(end)],repmat(mean(shuffle_coupl.(events(i)).distr(phase_ind,:),2),2,1),'Color',0.7*[1,1,1],'LineWidth',1.3,'DisplayName','shuffled');
  hand(2) = plot(axs(2*i-1),[phi;phi+phi(end)],repmat(circularSmooth(distr_coupl.(events(i))(phase_ind),[],'gaussian',10),2,1),'Color',paperColors(9+i),'LineWidth',1.3,'DisplayName',"distribution");
  title(axs(2*i-1),"coupled, p: "+stats_coupl.(events(i)).p)

  hand(1) = plot(axs(2*i),[phi;phi+phi(end)],repmat(mean(shuffle_uncoupl.(events(i)).distr(phase_ind,:),2),2,1),'Color',0.7*[1,1,1],'LineWidth',1.3,'DisplayName','shuffled');
  hand(2) = plot(axs(2*i),[phi;phi+phi(end)],repmat(circularSmooth(distr_uncoupl.(events(i))(phase_ind),[],'gaussian',10),2,1),'Color',paperColors(9+i),'LineWidth',1.3,'DisplayName',"distribution");
  legend(axs(2*i),hand), title(axs(2*i),"uncoupled, p: "+stats_uncoupl.(events(i)).p)

end
ylabel(axs(1),'\langleevent distribution\rangle'), ylim(axs,[0,0.5]), yticks(axs,[0.1,0.3])
clear i axs phase_ind a hand
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/ISRipplesSpindles/3f_coupl_entr_'+task,["png","svg"]);



%% TRY RIPPLE RESPONDERS
[isr_peth,trans,t] = runBatch(batch_file,@ISRippleRespondersNEW_,{["pfc"],[]},'sessions',15);



%% compare percentages across different sleeps

window = 0.5; % s
tasks = ["sleepm","sleepn","sleeps"];
duration = 0;
for task = tasks
  args = {window,'event',task,'state','sws','duration',duration};
  [phase.(task),is_coupl.(task),not_us.(task),us_phase.(task)] = runBatch(batch_file,@USRipSpin_,args,verbose=true);
end

ind = ~cellfun(@isempty,is_coupl.sleepn);
for task = tasks
  is_coupl_pool.(task) = reverseCellStruct(is_coupl.(task)(ind));
  is_coupl.(task) = reverseCellStruct(is_coupl.(task)(ind),@(x) x,@(x) x);
  not_us_pool.(task) = reverseCellStruct(not_us.(task)(ind));
  not_us.(task) = reverseCellStruct(not_us.(task)(ind),@(x) x,@(x) x);
end
events = string(fieldnames(is_coupl_pool.sleepm))';

for task = tasks
  for i = 1 : numel(events)
    A = cellfun(@(x) numel(x),is_coupl.(task).(events(i)));
    B = cellfun(@(x) x(1),not_us.(task).(events(i)));
    event_f.(task).(events(i)) = A ./ (A+B);
  end
end

clear ind task



%% TRY ANOVA
i=1; % Ripples
data_us = []; data_nus = [];
  for task = tasks
    data_us = [data_us,cellfun(@(x) mean(x),is_coupl.(task).(events(i)))]; % ./ event_f.(task).(events(3-i))];
    data_nus = [data_nus,cellfun(@(x) x(2)/x(1),not_us.(task).(events(i)))]; % ./ (1-event_f.(task).(events(3-i)))];
  end

  [p,tbl,stt] = anova2(data_us,1,'off');
  multcompare(stt,'Display','off')

%% I COULD JUSTIFY DIVIDING BY N OF EVENTS WITH A GILLESPIE SIMULATION WHERE I INCREASE THE N OF RIP AND SEE WHAT P(COUPL) DOES
[~,axs] = makeFigure('coupl',"P of event being coupled, w: "+num2str(window,3)+" s",[1,2]);
p = []; star_coord = [];
for i = 1 : numel(events)
  data_us = []; data_nus = [];
  for task = tasks
    data_us = [data_us,cellfun(@(x) mean(x),is_coupl.(task).(events(i)))]; % ./ event_f.(task).(events(3-i))];
    data_nus = [data_nus,cellfun(@(x) x(2)/x(1),not_us.(task).(events(i)))]; % ./ (1-event_f.(task).(events(3-i)))];
  end

  % generalized linear model
  y_us = data_us.';
  y_nus = data_nus.';
  n = numel(tasks);
  x = [zeros(1,n-1);diag(ones(n-1,1))];
  x = repmat(x,size(y_us,2),1);
  glm_us = fitglm(x,y_us(:),'CategoricalVars',[1,2]);
  glm_nus = fitglm(x,y_nus(:),'CategoricalVars',[1,2]);

  plot(axs(i),data_us.'*100,'Color',[paperColors(9+i),.3],'LineWidth',1)
  plot(axs(i),data_nus.'*100,'Color',[.6,.6,.6,.3],'LineWidth',1)

  coeff = glm_us.Coefficients.Estimate * 100;
  coeff(2:end) = coeff(2:end) + coeff(1);
  hand(1) = errorbar(axs(i),coeff,glm_us.Coefficients.SE*100,'Color',paperColors(9+i),'LineWidth',1.3,'CapSize',0,'DisplayName','in ISR');
  star_coord = [star_coord;[1;2;3],coeff*1.1];
  coeff = glm_nus.Coefficients.Estimate * 100;
  coeff(2:end) = coeff(2:end) + coeff(1);
  hand(2) = errorbar(axs(i),coeff,glm_nus.Coefficients.SE*100,'Color',[.6,.6,.6],'LineWidth',1.3,'CapSize',0,'DisplayName','not ISR');
  p = [p,glm_us.Coefficients.pValue,glm_nus.Coefficients.pValue];
  star_coord = [star_coord;[1;2;3],coeff*1.1];

  xlim(axs(i),[0.5,3.5]), xticks(axs(i),1:numel(tasks)), xticklabels(axs(i),tasks), title(axs(i),events(i)), legend(axs(i),hand)
end
ylim(axs(1),[5,35]),  ylim(axs(2),[25,55]), ylabel(axs(1),'percentage (%)')

h = holmBonferroni(p(:));
hand = scatter(axs(1),star_coord(logical(h(1:6)),1),star_coord(logical(h(1:6)),2),'k*');
hand(2) = scatter(axs(2),star_coord([false(6,1);logical(h(7:end))],1),star_coord([false(6,1);logical(h(7:end))],2),'k*');
RemoveFromLegend(hand)

% for i = 1 : numel(events)
%   data_us = []; data_nus = [];
%   for task = tasks
%     d = is_coupl_pool.(task).(events(i));
%     data_us = [data_us;sum(d),numel(d)];
%     data_nus = [data_nus;fliplr(sum(not_us_pool.(task).(events(i))))];
%   end
%   b = bar(axs(i),[data_us(:,1)./data_us(:,2)*100,data_nus(:,1)./data_nus(:,2)*100],'FaceColor','flat','CData',[0,0,0],'LineWidth',1.1);
%   b(1).FaceColor = paperColors(9+i);
%   b(2).FaceColor = [0.6,0.6,0.6];
%   xticks(axs(i),1:numel(tasks)), xticklabels(axs(i),tasks), ylabel(axs(i),'percentage (%)')
%   %[p,~,z] = twoProportionsZ(data_us(1),data_us(2),data_nus(1),data_nus(2));
% end

clear i d data_us data_nus n p x y_us y_nus h star_coord coeff task hand glm_us glm_nus axs
saveFlag && saveFig(gcf,'ReuSlowRythm/Results/FinalFigures/USRipplesSpindles/4coupling_task_effect',["png","svg"]);