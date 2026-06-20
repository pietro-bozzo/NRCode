function obj = cell2stager(cell)
    
arguments
    cell (:,1) {mustBeA(cell,'cell')}
end

for i = 1:length(cell)
    mustBeA(cell{i},'stager')
end

obj = [cell{:}];

obj = transpose(obj);

end