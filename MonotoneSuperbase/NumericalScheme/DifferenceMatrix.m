function data = DifferenceMatrix(data)
% Creates matrix and constant term for the first and second order differences in a grid.
% Options must contain : offsets, indices, multipliers,

% The product A1 x+b1 yields all first order differences of x.
% The product A2 x+b2 yields all second order differences of x.

assert(isfield(data.geometry,'points'));
assert(isfield(data.geometry,'offsets'));
assert(isfield(data.geometry,'neighIndices'));
assert(isfield(data.geometry,'multipliers'));
assert(isfield(data,'boundaryValue'));

Dimension   = size(data.geometry.points,1);
nPoints     = size(data.geometry.points,2);
nOffsets    = size(data.geometry.offsets,2); % Value changed below appropriately

assert(all(size(data.geometry.points)==[Dimension,nPoints]));
assert(all(size(data.geometry.offsets)==[Dimension,nOffsets]));
assert(all(size(data.geometry.neighIndices)==[1,nPoints,nOffsets]));
assert(all(size(data.geometry.multipliers)==[1,nPoints,nOffsets]));
% boundaryValue is expected to be a pure function

% ----- Construct b1 -----
nOffsets = size(data.offsets1,2);

if nOffsets>0 
    nDiff = nPoints*nOffsets;
    [loc,pos] = ismember(data.offsets1,data.geometry.offsets);
    assert(all(loc));

    indices = reshape(data.geometry.neighIndices(1,:,pos),[1,nDiff]);
    multipliers = data.geometry.multipliers(1,:,pos);

    % Construct points on the boundary.
    boundary = repmat(data.geometry.points,[1,1,nOffsets]) + ...
        repmat(multipliers,[Dimension,1,1]) .* ...
        repmat( reshape(data.offsets1,[Dimension,1,nOffsets]), [1,nPoints,1]);

    multipliers = reshape(multipliers,[1,nDiff]);
    boundary = reshape(boundary,[Dimension,nDiff]);

    pos = indices==0;

    b1 = zeros([nDiff,1]);
    b1(pos) = data.boundaryValue(boundary(pos)) ./ multipliers(pos);

    data.linear.b1 = b1;
% ----- Construct A1 -----

    Rows = [1:nDiff,1:nDiff];
    Cols = [indices, reshape(repmat(1:nPoints,[1,1,nOffsets]),[1,nDiff])];
    Vals = [1./multipliers, -1./multipliers];

    pos = Cols~=0;
    Rows = Rows(pos); Cols=Cols(pos); Vals=Vals(pos);
    data.linear.A1 = sparse(Rows,Cols,Vals,nDiff,nPoints);
else
    data.linear.A1 = zeros([0,nPoints]);
    data.linear.b1 = zeros([0,1]);
end

% ----- Construct b2 -----
nOffsets = size(data.offsets2,2);
if nOffsets>0
    nDiff = nPoints*nOffsets;
    
    % Positive offsets
    [loc,pos] = ismember( data.offsets2', data.geometry.offsets', 'rows');
    assert(all(loc));
    
    multP = data.geometry.multipliers(1,:,pos);    
    boundaryP = repmat(data.geometry.points,[1,1,nOffsets]) + ...
        repmat(multP,[Dimension,1,1]) .* ...
        repmat( reshape(data.offsets2,[Dimension,1,nOffsets]), [1,nPoints,1]);
    
    multP = reshape(multP,[1,nDiff]);
    boundaryP = reshape(boundaryP,[Dimension,nDiff]);
    indicesP = reshape(data.geometry.neighIndices(1,:,pos),[1,nDiff]);
    
    
    % Negative offsets
    [loc,pos] = ismember(-data.offsets2', data.geometry.offsets', 'rows');
    assert(all(loc));
    
    multN = data.geometry.multipliers(1,:,pos); 
    boundaryN = repmat(data.geometry.points,[1,1,nOffsets]) + ...
        repmat(multN,[Dimension,1,1]) .* ...
        repmat( reshape(-data.offsets2,[Dimension,1,nOffsets]), [1,nPoints,1]);
    
    multN = reshape(multN,[1,nDiff]);
    boundaryN = reshape(boundaryN,[Dimension,nDiff]);
    indicesN = reshape(data.geometry.neighIndices(1,:,pos),[1,nDiff]);

    
    % Mean of multipliers
    multAvg = (multP+multN)/2.;
    
    % Constant term in second order diffs ((up - u)/multP + (uN-u)/multN)/multAvg
    
    b2 = zeros([1,nDiff]);
    pos = indicesP==0;
    b2(pos) = b2(pos)+data.boundaryValue(boundaryP(:,pos)) ./ (multAvg(pos).*multP(pos));
    pos = indicesN==0;
    b2(pos) = b2(pos)+data.boundaryValue(boundaryN(:,pos)) ./ (multAvg(pos).*multN(pos));
    data.linear.b2 = b2';
    
    % ----- Construct A2 -----
    
    Rows = [1:nDiff,1:nDiff,1:nDiff];
    Cols = [indicesP,indicesN,reshape(repmat(1:nPoints,[1,1,nOffsets]),[1,nDiff])];
    Vals = [1./multP,1./multN, -(1./multP+1./multN)]./repmat(multAvg,[1,3]);
    
    pos = Cols~=0;
    Rows = Rows(pos); Cols=Cols(pos); Vals=Vals(pos);
    data.linear.A2 = sparse(Rows,Cols,Vals,nDiff,nPoints);
else
    data.linear.A2 = zeros([0,nPoints]);
    data.linear.b2 = zeros([0,1]);
end

end