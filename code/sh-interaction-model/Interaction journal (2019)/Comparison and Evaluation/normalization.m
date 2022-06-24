function normAbnormal = normalization(Abnormal)
minDataNorm = min(Abnormal);
dataNorm = Abnormal - repmat(minDataNorm,size(Abnormal,1),1);
maxDataNorm = max(dataNorm);
normAbnormal   = dataNorm./repmat(maxDataNorm,size(Abnormal,1),1);
end

