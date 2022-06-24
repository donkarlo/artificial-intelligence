function [norm, mindataNorm,maxdataNorm] = Normalization(data)
mindataNorm = min(data);
newdata = data - mindataNorm;
maxdataNorm = max(newdata);
norm = newdata/maxdataNorm;
end

