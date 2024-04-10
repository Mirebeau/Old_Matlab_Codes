% As such, speed up is not sufficient to justify for the added complexity


mex -Dchar16_t=uint16_T -output TripletMatrixProduct Toolbox/mex/TripletMatrixProduct.cpp

triplets = [1,2,1;2,2,1;1,2,3];
vector = [1,2;3,4];
times = [1,2];
[result,done]=TripletMatrixProduct(triplets,vector,times);
A=sparse(triplets(1,:),triplets(2,:),triplets(3,:),2,2);
result2=A*vector;

%        [result,dummy]=TripletMatrixProduct(triplets,smoothed,[1,2]);
%        smoothed = smoothed - dt*result;


%    A=sparse(triplets(1,:),triplets(2,:),triplets(3,:),prod(s),prod(s));
%  triplets = [row;column;coefficient];
