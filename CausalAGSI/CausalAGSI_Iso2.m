%Input : 
% - a Speed 2D array, defining the isotropic metric M = Id/Speed^2.
%  non-positive values of the Speed field act as walls.
% - seeds, given by their position and value [i,j,v]

%Copyright Jean-Marie Mirebeau, 2015 jm.mirebeau@gmail.com

function output = CausalAGSI_Iso2(input)
    s = size(input.Speed);
    
    output.Distance = zeros(s)+Inf;
    output.UpdateCount = zeros(s);
    improvementSinceLastAdded=zeros(s)+Inf;

    front=input.Seeds(1:2,:);
    output.Distance( sub2ind(s,front(1,:),front(2,:)) )=input.Seeds(3,:);
    
    while size(front)>0
        % Compute distances for front members
        val1=zeros([1,size(front,2)])+Inf; val2=val1;
        S1p = front(1,:)<s(1); val1(S1p) =          output.Distance(sub2ind(s,front(1,S1p)+1,front(2,S1p)));
        S1m=front(1,:)>1; val1(S1m) = min(val1(S1m),output.Distance(sub2ind(s,front(1,S1m)-1,front(2,S1m))));
        S2p=front(2,:)<s(2); val2(S2p) =            output.Distance(sub2ind(s,front(1,S2p),  front(2,S2p)+1));
        S2m=front(2,:)>1; val2(S2m) = min(val2(S2m),output.Distance(sub2ind(s,front(1,S2m),  front(2,S2m)-1)));
        val = val1; val1=min(val1,val2); val2=max(val,val2);
        
        indices = sub2ind(s,front(1,:),front(2,:));
        invSpeed = 1./input.Speed(indices);
        
        val = output.Distance(indices);
        val = min(val,val1+invSpeed);
        
        delta = 2*invSpeed - (val1-val2).^2;
        delta(delta<0)=Inf;
        delta=0.5*(val1+val2+sqrt(delta));
        delta(delta<val2)=Inf;
        val = min(val,delta);
        
        %Setup next front.
        threshold = min(val); % TO DO : mean or min, or smthg else ?
        improvement = output.Distance(indices)-val;
        improvementSinceLastAdded(indices)=improvementSinceLastAdded(indices)+improvement;
        
        output.Distance(indices)=val;
        output.UpdateCount(indices)=output.UpdateCount(indices)+1;
        
        kept = improvementSinceLastAdded(indices)>input.tol;
        added = kept & val<(threshold+invSpeed); %TO DO: add only neighbors with higher value
        % TO DO remove added from kept. 
        % TO DO maintain a silent kept, until threshold+speed fits.
        %added = kept; %basic AGSI
        if(all(~(kept|added))) break; end; %TO DO : increase instead factor in front invSpeed instead
        
        improvementSinceLastAdded(indices(added))=0;
        newFront = [front(1,kept), front(1,added & S1p)+1, front(1, added & S1m)-1, front(1,added & S2p), front(1,added & S2m); ...
                    front(2,kept), front(2,added & S1p),   front(2, added & S1m),   front(2,added & S2p)+1, front(2,added & S2m)-1];
        
        %Delete duplicates.
        newIndices = sub2ind(s,newFront(1,:),newFront(2,:));
        [~,newIndices]=cmunique(newIndices);
        indices=newIndices(:,1);
        [front1,front2]=ind2sub(s,indices);
        front=zeros([2,size(front1)]);
        front(1,:)=front1; front(2,:)=front2;
    end
end
