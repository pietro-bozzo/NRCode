function stgr = timer2stager(obj)

arguments
    obj (1,1) timer {mustBeA(obj,'timer')}
end

if isempty(obj)
    stgr = stager.empty;
    return
end

cellstgr = obj.UserData;

stgr = cell2stager(cellstgr);

end

