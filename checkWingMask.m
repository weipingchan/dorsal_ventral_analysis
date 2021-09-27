function outMask=checkWingMask(inMask, in_key, part)
    if nnz(inMask)<=10 %If one wing is missing, provide a fake one
        if part=='LF'
            keyPts=[in_key(7,:) ; in_key(2,:)];
            fakepolygon=[in_key(1,:); keyPts; in_key(10,:)];
        elseif part=='RF'
            keyPts=[in_key(9,:) ; in_key(4,:)];
            fakepolygon=[in_key(3,:); keyPts; in_key(11,:)];
        elseif part=='LH'
            keyPts=[in_key(7,:) ; in_key(8,:)];
            fakepolygon=[in_key(1,:); keyPts];
        elseif part=='RH'
            keyPts=[in_key(9,:) ; in_key(5,:)];
            fakepolygon=[in_key(4,:); keyPts];
        end
        outMask=roipoly(inMask,fakepolygon(:,1),fakepolygon(:,2));
    else
        outMask=inMask;
    end
end