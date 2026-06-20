classdef stager    
    properties
        sws (:,2) double {mustBeNumeric}
        rem (:,2) double {mustBeNumeric}
        drowsy (:,2) double {mustBeNumeric}
        rest (:,2) double {mustBeNumeric}
        move (:,2) double {mustBeNumeric}
        tags (1,1) struct {mustBeA(tags,'struct')}
        data (1,1) struct {mustBeA(data,'struct')}
    end
    
    methods       
        function obj = stager(sws,rem,drowsy,rest,move,session,channel,dims,opt)
            arguments
                % brain states
                sws (:,2) double {mustBeNumeric}
                rem (:,2) double {mustBeNumeric}
                drowsy (:,2) double {mustBeNumeric}
                rest (:,2) double {mustBeNumeric}
                move (:,2) double {mustBeNumeric}
                % tags
                session (1,:) char {mustBeTextScalar}
                channel (:,1) double {mustBeNumeric} = NaN
                dims (1,1) double {mustBeNumeric} = NaN
                % datas
                opt.speed (:,2) double {mustBeNumeric} = double.empty
                opt.y (:,:) double {mustBeNumeric} = double.empty
                opt.t (:,1) double {mustBeVector} = double.empty
                opt.f (:,1) double {mustBeVector} = double.empty
            end
            
            if isempty(fileparts(session))
                error('session must be a directory')
            end

            % brain states
            obj.sws = sws;
            obj.rem = rem;
            obj.drowsy = drowsy;
            obj.rest = rest;
            obj.move = move;  
            % tags
            obj.tags.session = session;
            obj.tags.channel = channel;
            obj.tags.full = all([~isempty(opt.speed),~isempty(opt.y), ...
                ~isempty(opt.t),~isempty(opt.f)]);
            obj.tags.dims = dims;
            % datas
            obj.data.speed = opt.speed;
            obj.data.spect.y = opt.y;
            obj.data.spect.t = opt.t;
            obj.data.spect.f = opt.f;
        end
        
        function [sws,rem,drowsy,rest,move] = extrbs(obj)
            arguments
                obj (1,1) {mustBeA(obj,'stager')}
            end
            sws = obj.sws;
            rem = obj.rem;
            drowsy = obj.drowsy;
            rest = obj.rest;
            move = obj.move;
        end
        
        function [session,channel] = extags(obj)
            arguments
                obj (1,1) {mustBeA(obj,'stager')}
            end
            session = obj.tags.session;
            channel = obj.tags.channel;
        end
        
        function [speed,y,t,f] = exdata(obj)
            arguments
                obj (1,1) {mustBeA(obj,'stager')}
            end
            speed = obj.data.speed;
            y = obj.data.spect.y;
            t = obj.data.spect.t;
            f = obj.data.spect.f;
        end
        
        function couple = stager2couple(obj,add)
            arguments
                obj (1,:) stager
                add (:,2) cell {mustBeCouple(add,'stager','all')} = cell.empty
            end
            if isempty(obj)
                couple = {'',stager.empty};
            else
                couple = num2cell(zeros(length(obj),2));
                for i = 1:length(obj)
                    [~,id] = fileparts(obj(i).tags.session);
                    couple{i,1} = char(id);
                    couple{i,2} = obj(i);
                end
            end
            if ~isempty(add)
                couple = vertcat(add,couple);
            end
        end
        
        function cellstgr = stager2cell(obj)
            arguments
                obj (:,1) stager {mustBeA(obj,'stager')}
            end
            if isempty(obj)
                cellstgr = cell.empty;
                return
            end
            cellstgr = num2cell(zeros(length(obj),1));
            for i = 1:length(obj)
                cellstgr{i,1} = obj(i);
            end
        end
        
        function idx = tagfinder(obj,id)     
            arguments
                obj (:,1) stager {mustBeA(obj,'stager')}
                id (1,1) string {mustBeTextScalar}
            end
            mustBeNonempty(obj)
            idx = double.empty;
            for i = 1:length(obj)
                [~,tagid] = fileparts(obj(i).tags.session);
                if string(tagid) == id
                    idx = [idx i];
                end
            end
            idx = transpose(idx);
        end 

        function savestager(obj)
            arguments
                obj (:,1) stager {mustBeA(obj,'stager')}               
            end
            if length(obj) == 1
                svstgr(obj)
            elseif length(obj) > 1
                for ii = 1:length(obj)
                    svstgr(obj(ii))
                end
            end
            function svstgr(obj)
            [directory, ID] = fileparts(obj.tags.session);
            state = ["sws","rem","drowsy","rest","move"];
            statename = ["sws","rem","drowsiness","rest","movement"];           
            for i = 1:length(state)
                filename = fullfile(directory,strcat(ID,'.',statename(i)));              
                writematrix(obj.(state(i)),filename,FileType="text");
            end            
            end
        end        
    end
end