function stats = stagerstat(obj,plt)

arguments
    obj (:,1) stager {mustBeA(obj,'stager')}
    plt (1,1) logical {mustBeMember(plt,[0,1])} = false;
end

if length(obj)==1
    stats = stgrSt(obj,plt);     
else
    z = double.empty;
    count = 0;
    stats = table;
    for h = 1:length(obj)
        count = count + 1;
        x = stgrSt(obj(h),plt);
        stats(h,1) = {x};
        [~,tagid] = fileparts(obj(h).tags.session);
        tagid = strsplit(tagid,'_');
        stats.Properties.RowNames{h} = tagid{2};     
    end
    stats.Properties.VariableNames = {'stats'};             
end


%% stager statistics

function [stats,y] = stgrSt(obj,plt)

[session,channel] = extags(obj);
[sws,rem,drowsy,rest,move] = extrbs(obj);
[~,~,evt] = getEvt(session);
id = string(fieldnames(evt));
stn = ["sws","rem","drowsy","rest","move"];
states = {sws,rem,drowsy,rest,move};
for i = 1:length(id)
    int = evt.(id(i));
    for k = 1:length(states)
        s = states{k};
        slim = SubtractIntervals(s,[-inf,int(1);int(2),inf]);
        y(i,k) = sum(diff(slim,1,2),"all");
    end  
end

[~,ID] = fileparts(session);

if plt
figure(Name=append(strrep(ID,'_',' '),'_channel_', string(channel)));
sgtitle(append(strrep(ID,'_',' '),' channel ', string(channel)));
cc = [0 0.4470 0.7410;0.6350 0.0780 0.1840;0.4940 0.1840 0.5560; ...
    0.3010 0.7450 0.9330;0.4660 0.6740 0.1880];
patch = bar(1:size(y,1),y,'stacked');
for i = 1:length(patch)
    patch(i).CData = cc;
    patch(i).FaceAlpha = 0.7;
end
xticklabels(id)
legend(stn,location='northeastoutside')
xlabel('event')
ylabel('duration (s)')
end

yy = [ y(:,1)+y(:,2),y];
yy = [yy, sum(yy,2)];
yy = [yy; sum(yy,1)];

yy = yy';
yy = round(yy,3);
stats = array2table(yy);

stats.Properties.RowNames = cellstr(["sleep",stn,"tot"]);
stats.Properties.VariableNames = cellstr([id;"tot"]);

end

end
