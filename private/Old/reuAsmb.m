% Analyse assemblies activity in Nucleus Reuniens data

% set up batch
batch_file = '/mnt/hubel-data-103/Pietro/Data/BatchFiles/sessions_reu.txt';

% choose parameter values
event = "tachem";
state = "sws";
regs = [12;32;71];
ICA_window = 0.03;
plt = false;
args = {'event',event,'state',state,'regions',regs,'ICA_window',ICA_window,'plot',plt,'verbose',false};

% run batch
[activ_rate,weights,neurons] = runBatch(batch_file,@reuAsmb_,args,verbose=true);

%% pool outputs
contrib = [];
contrib_no_eff = [];
for i = 1 : size(activ_rate,1)
  % activation rates pre- and post-sleep
  session_activ_rate = activ_rate{i};
  % memeber participations per region
  members = abs(weights{i}) > (1 ./ sqrt(sum(weights{i}~=0,1))); % binarize for assembly belonging matrix
  memebers_eff = members(:,session_activ_rate(:,2)>session_activ_rate(:,1)); % keep ICs which activate more after task
  memebers_no_eff = members(:,session_activ_rate(:,2)<=session_activ_rate(:,1));
  reg_ind = zeros(size(neurons{i},1),2);
  session_contrib = zeros(size(neurons{i},1),size(memebers_eff,2));
  session_contrib_no_eff = zeros(size(neurons{i},1),size(memebers_no_eff,2));
  for j = 1 : size(neurons{i},1)
    reg_ind(j,:) = [neurons{i}{j}(1),neurons{i}{j}(end)];
    session_contrib(j,:) = sum(memebers_eff(reg_ind(j,1):reg_ind(j,2),:));
    session_contrib_no_eff(j,:) = sum(memebers_no_eff(reg_ind(j,1):reg_ind(j,2),:));
  end
  contrib = [contrib;session_contrib.'];
  contrib_no_eff = [contrib_no_eff;session_contrib_no_eff.'];
end
contrib = contrib(sum(contrib==0,2)~=3,:); % remove assemblies with no accepted members
contrib_no_eff = contrib_no_eff(sum(contrib_no_eff==0,2)~=3,:);
clear members memebers_eff memebers_no_eff session_activ_rate session_contrib session_contrib_no_eff % reg_ind weights

%% show that asmbs activate more after task
% pool all assemblies together
activ_rate = vertcat(activ_rate{:});
figure(Position=get(0,'ScreenSize'),Name='box',NumberTitle='off'); hold on
groups = [ones(size(activ_rate,1),1);repmat(2,size(activ_rate,1),1)];
boxplot([activ_rate(:,1);activ_rate(:,2)],groups);
title('Whole-brain-assemblies'' activation rate',FontSize=17,FontWeight='normal')
ylabel('activation rate (Hz)',FontSize=14)
set(gca,TickDir='out',Box='off')
[~,p] = ttest(activ_rate(:,1),activ_rate(:,2));
text(0.5,1.5,append('ttest p: ',string(p),' n: ',string(size(activ_rate,1))))
xticklabels(["pre-task sleep","post-task sleep"])

%% show the same with pair plot
figure(Position=get(0,'ScreenSize'),Name='pair',NumberTitle='off'); hold on
groups = repmat([1;2;NaN],size(activ_rate,1),1);
values = [activ_rate nan(size(activ_rate,1),1)].';
plot(groups,values(:),color=[0.5,0.5,0.5]);
plot([1,2],mean(activ_rate),color=[0,0,0]);
xlim([0.5,2.5]); ylim([-0.25,2]); set(gca,TickDir='out',Box='off')
clear groups values

%% show participations in assemblies which reactivate more after task
n_reu = sum(sum(contrib(:,1:2)==0,2)==2); % number of asmb having only reuniens members
n_reu_pfc = sum(contrib(:,2)==0 & sum(contrib(:,[1,3])~=0,2)==2);
n_reu_hpc = sum(contrib(:,3)==0 & sum(contrib(:,[1,2])~=0,2)==2);
n_reu_all = sum(sum(contrib~=0,2)==3);
n_other = sum(contrib(:,1)==0);
figure(Position=get(0,'ScreenSize'),Name='pie',NumberTitle='off');
piechart([n_reu,n_reu_pfc,n_reu_hpc,n_reu_all,n_other],["only REU","REU + mPFC","REU + HPC","REU + mPFC + HPC","no REU"])
title('Effect membership proportions')

n_reu = sum(sum(contrib_no_eff(:,1:2)==0,2)==2); % number of asmb having only reuniens members
n_reu_pfc = sum(contrib_no_eff(:,2)==0 & sum(contrib_no_eff(:,[1,3])~=0,2)==2);
n_reu_hpc = sum(contrib_no_eff(:,3)==0 & sum(contrib_no_eff(:,[1,2])~=0,2)==2);
n_reu_all = sum(sum(contrib_no_eff~=0,2)==3);
n_other = sum(contrib_no_eff(:,1)==0);
figure(Position=get(0,'ScreenSize'),Name='pie',NumberTitle='off');
piechart([n_reu,n_reu_pfc,n_reu_hpc,n_reu_all,n_other],["only REU","REU + mPFC","REU + HPC","REU + mPFC + HPC","no REU"])
title('No effect membership proportions')