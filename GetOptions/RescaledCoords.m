function [NewCoords] = RescaledCoords(OldCoords,Origin,Spacing,Direct)
    assert(nargin==3 || nargin==4);
    if nargin==3; Direct=true; end;
    assert(size(OldCoords,1)==size(Origin,1));
    assert(size(Origin,2)==1);
    assert(all(size(Origin)==size(Spacing)));
    
    len = size(OldCoords); len=len(2);
    Ones = ones(1,len);
    Origins = Origin*Ones;
    Spacings = Spacing*Ones;
    
    if Direct
        NewCoords = 1+(OldCoords-Origins)./Spacings;
    else
        NewCoords = Origins+(OldCoords-1).*Spacings;
    end
end