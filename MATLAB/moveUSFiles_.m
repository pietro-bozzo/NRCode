function moveUSFiles_(session)

[filepath,basename] = fileparts(session);
movefile(filepath+"/UltraSlowRythm",filepath+"/InfraSlowRhythm")

% R0 = regions(session,'legend',dataPath+"/bilateral.anat");
% neurons = R0.nNeurons;
% 
% [file_root,basename] = fileparts(session);
% labels = lower(R0.ids);
% file_root = fullfile(file_root,"InfraSlowRhythm",basename);
% extensions = join([repmat(".slow",size(labels)),labels],'');
% extensions_aval = join([repmat(".slowaval",size(labels)),labels],'');
% 
% for i = 1 : numel(labels)
%   try
%     intervals.(labels(i)) = readmatrix(file_root+extensions(i),FileType='text',CommentStyle='%');
%     slow_avals.(labels(i)) = readmatrix(file_root+extensions_aval(i),FileType='text',CommentStyle='%');
%     if neurons(i) < 30
%       saveMatrix([NaN,NaN],file_root+extensions(i),'beginnig','end of slow-rhythm interval');
%       saveMatrix([NaN,NaN,NaN],file_root+extensions_aval(i),'beginnig','end','size of slow-rhythm avalanches');
%     end
%   catch
%   end
% end