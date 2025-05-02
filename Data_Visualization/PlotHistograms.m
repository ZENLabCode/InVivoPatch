edges = -650:100:-50;
centers = edges(1:end-1) + diff(edges)/2;
counts = histcounts(A,edges) / numel(A)*100;

histX = centers(:);
histY = counts(:);

xq = linspace(min(histX), max(histX), 100);

yq = spline(histX, histY, xq);