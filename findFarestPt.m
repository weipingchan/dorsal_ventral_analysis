function farestPt=findFarestPt(ptList,tarPt)
    %compute Euclidean distances:
    distances = sqrt(sum(bsxfun(@minus, ptList,tarPt).^2,2));
    %find the smallest distance and use that as an index into B:
    farestPt0=ptList(distances==max(distances),:);
    if size(farestPt0,1)>1
        farestPt=farestPt0(1,:);
    else
        farestPt=farestPt0;
    end
end