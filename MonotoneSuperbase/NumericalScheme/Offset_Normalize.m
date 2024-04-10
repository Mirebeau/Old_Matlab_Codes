function offsets = Offset_Normalize(offsets)
% An offset is normalized iff its first non-zero coordinate is non-negative.
Dimension = size(offsets,1);
nPoints = size(offsets,2);
assert(all(size(offsets) == [Dimension,nPoints]));

done = zeros([1,nPoints]);
for k=1:Dimension
    active = offsets(k,:)<0 & ~done;
    offsets(:,active) = -offsets(:,active);
    done = offsets(k,:)>0 | done;
end

end