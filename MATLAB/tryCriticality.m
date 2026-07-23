froot = fileparts(fileparts(fileparts(matlab.desktop.editor.getActiveFilename)));
saveFlag = false;

% set up batch
batch_file = fullfile(froot,'InfraSlowNRPaper/Data/IS_units.batch');

region = "nr";
states = ["all", "other", "rem", "sws", "sws"];
event_list = ["InfraSlowRhythm/slownr", "all", "all", "all", "^InfraSlowRhythm/slownr"];
window = 0; % WAS 0.01; % s
args = {region,states,event_list,window,'threshold',0,'labels',["ISR","wake","REM","nREM","nISR"]};

[alpha,a_range,a_D,beta,b_range,b_D,gam_th,gam_area,gam_shape,chi] = runBatch(batch_file,@CriticalExponents,args,'ignore_args',true,'sessions',10); % ,'sessions',10
alpha = reverseCellStruct(alpha);
a_range = reverseCellStruct(a_range);
a_D = reverseCellStruct(a_D);
beta = reverseCellStruct(beta);
b_range = reverseCellStruct(b_range);
b_D = reverseCellStruct(b_D);
gam_th = reverseCellStruct(gam_th);
gam_area = reverseCellStruct(gam_area);
gam_shape = reverseCellStruct(gam_shape);
chi = reverseCellStruct(chi);
diffgam = structFun(@minus, gam_area, gam_th);
stats_list = {alpha,a_range,a_D,beta,b_range,b_D,gam_th,gam_area,diffgam,gam_shape,chi};

%%
[fig,axs] = makeFigure('crit',upper(region),[2,6],'size',[30,10]);
stats_names = ["\alpha","range \alpha","D \alpha","\beta","range \beta","D \beta","\gamma_{th}","\gamma_{area}","DCC","\gamma_{shape}","\chi"];
for i = 1 : numel(stats_list)
  data = structCat(stats_list{i});
  [~,info] = boxPlot(structCat(stats_list{i},2),'label',fieldnames(stats_list{i}),'ax',axs(i));
  y_lim = [.9*min(info.lw),1.1*max(info.uw)];
  if ~any(isnan(y_lim))
    ylim(axs(i),y_lim)
  end
  title(axs(i),stats_names(i))
end

clearvars stats_names fig axs data group i b x_tick y_lim