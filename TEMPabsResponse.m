function [response,delay] = TEMPabsResponse(peth,t,limit)
  
  ind = find(t>limit,1);
  m = mean(peth(1:ind-1)); % mean over baseline
  s = std(peth(1:ind-1));
  [~,delay] = max(abs(peth(ind:end)-m));
  response = (peth(delay+ind-1) - m) / s; % standardize w.r.t. baseline
  delay = t(delay+ind-1); % convert to s

end

function AA = aaaa()
% old algo
figure, hold on
plot(t,peth)
plot(t,peth-m)
plot(t,abs(peth-m))
yline(max(abs(peth-m)),'r')
yline(0,'k')
title("delta: "+string(max(abs(peth-m)))+", r: "+string(response))

% new algo
figure, hold on
plot(t,peth)
plot(t,peth-m)
plot(t,abs(peth-m))
yline(max(abs(peth(ind:end)-m)),'r')
yline(0,'k')
title("delta: "+string(max(abs(peth-m)))+", r: "+string(response)+", d: "+delay+" s")
end