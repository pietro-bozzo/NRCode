function obj = text2stager(session)

arguments
    session (:,1) string {mustBeFile}
end

for j = 1:length(session)
    [directory, ID] = fileparts(session(j));
    state = ["sws","rem","drowsiness","rest","movement"];
    s = struct;
    for i = 1:length(state)
        filename = fullfile(directory,strcat(ID,'.',state(i)));        
        s.(state(i)) = dlmread(filename);
    end
    
    if isfile(fullfile(directory,strcat(ID,".speedpower")))
        speed = dlmread(fullfile(directory,strcat(ID,".speedpower")));
        obj(j) = stager(s.(state(1)),s.(state(2)),s.(state(3)),s.(state(4)),s.(state(5)),session(j),speed=speed);
    else
        obj(j) = stager(s.(state(1)),s.(state(2)),s.(state(3)),s.(state(4)),s.(state(5)),session(j));
    end
end