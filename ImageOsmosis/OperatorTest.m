addpath('../ParallelMatrix')

if false     % One dimensional sanity check
    options.dims = [10];
    options.gridScale=1;
    nPoints = prod(options.dims);
    m = ones(1,nPoints);
    w = ones(1,nPoints);

    [Rows,Cols,Vals]=ConvectionDiffusionSparseMatrix(m,w,options);
    A=sparse(Rows,Cols,Vals,nPoints,nPoints);

    [Rows,Cols,Vals]=ConvectionDiffusionSparseMatrix(m,[],options);
    B = sparse(Rows,Cols,Vals,nPoints,nPoints);
end


if false     % Reconstruction in one dimension
    x=linspace(0,1,100);
    options.dims = [numel(x)];
    options.gridScale = x(2)-x(1);
    options.upwind=false;

    %image = 1.+0.*x'; % OK
    image = exp(sin(2*pi*x))'; % OK
    %image = exp(0.+(x<=0.5))'; % Not OK due to discontinuity
    
    nPoints = prod(options.dims);
    m = ones(1,nPoints);
    w = gradient(log(image))'/options.gridScale;
    %w = gradient(image)./image; w = w'/options.gridScale; % Similar results

    [Rows,Cols,Vals]=ConvectionDiffusionSparseMatrix(m,w,options);
    A=sparse(Rows,Cols,Vals,nPoints,nPoints);
    
    % Check that the operator does have a kernel
    one = ones(nPoints,1);
    disp('Checking that transposed operator vanishes on constants.');
    disp(max(abs(one'*A)));
    
    sol = lsqlin(A,0*one,[],[],one',one'*image);
    disp('Relative error between the solution found and the original image.')
    disp(max(abs( (sol./image)-1)))
end

if true %  Reconstruction in two dimensions
    % Create axes
    %[x,y] = meshgrid(linspace(-1,1,21),linspace(0,1,11)); %OK
    [x,y] = meshgrid(linspace(-1,1,101),linspace(0,1,51));
    x=x'; y=y'; % Avoid Matlab's YXZ ordering
    
    % Create synthetic image
    %image = exp(x); % OK
    image = exp(-x.^2-y.^2); %OK
    
    % Set options
    options.upwind=true;
    options.dims=size(image);
    options.gridScale=x(2,1)-x(1,1);
    
    % Create tensor field (xx,xy,yy)
    m = zeros([3,options.dims]); 
    %m(1,:,:)=1; m(2,:,:) = 0; m(3,:,:)=1; %Identity matrix (xx,xy,yy)  = (1,0,1)
    m(1,:,:) = 1+x.*x; m(2,:,:)= x.*y; m(3,:,:)= 1+y.*y; % Mildly anisotropic
    
    % Compute image gradient
    [gradY,gradX] = gradient(log(image)/options.gridScale); % Again, YXZ ordering
    w=[shiftdim(gradX,-1); shiftdim(gradY,-1)];
    
    % Reshape for input
    nPoints = prod(options.dims);
    m = reshape(m,[3,nPoints]);
    w = reshape(w,[2,nPoints]);
    image = image(:);

    % Make system
    [Rows,Cols,Vals]=ConvectionDiffusionSparseMatrix(m,w,options);
    A=sparse(Rows,Cols,Vals,nPoints,nPoints);
    one = ones(nPoints,1);
    disp('Checking that transposed operator vanishes on constants.');
    disp(max(abs(one'*A)));
    
    % Solve
    sol = lsqlin(A,0*one,[],[],one',one'*image);
    disp('Relative error between the solution found and the original image.')
    disp(max(abs( (sol./image)-1)))
end

if false % Reconstruction in three dimensions
    %Create axes
    [x,y,z] = meshgrid(linspace(-1,1,21),linspace(0,1,11),linspace(-0.5,1,16)); %OK
    %[x,y,z] = meshgrid(linspace(-1,1,49),linspace(0,1,25),linspace(-0.5,1,37));   %OK, but longer solve
    x=permute(x,[2,1,3]); y=permute(y,[2,1,3]); z=permute(z,[2,1,3]);
    
    %Create synthetic image
    image = exp(-x.^2-y.^2-z.^2);
    
    % Set options
    options.upwind=true;
    options.dims=size(image);
    options.gridScale=x(2,1)-x(1,1);
    
    % Create tensor field (xx,xy,yy,xz,yz,zz)
    m = zeros([6,options.dims]);
    %m(1,:,:,:)=1; m(2,:,:,:)=0; m(3,:,:,:)=1; m(4,:,:,:)=0; m(5,:,:,:)=0; m(6,:,:,:)=1; % Identity matrix
    m(1,:,:,:)=1+x.*x; m(2,:,:,:)=x.*y; m(3,:,:,:)=1+y.*y; m(4,:,:,:)=x.*z; m(5,:,:,:)=y.*z; m(6,:,:,:)=1+z.*z; %Midly anisotropic
    
    % Compute image gradient
    [gradY,gradX,gradZ] = gradient(log(image)/options.gridScale); % Again, YXZ ordering
    w = [shiftdim(gradX,-1);shiftdim(gradY,-1);shiftdim(gradZ,-1)];

    % Reshape for input
    nPoints = prod(options.dims);
    m = reshape(m,[6,nPoints]);
    w = reshape(w,[3,nPoints]);
    image = image(:);

    % Make system
    [Rows,Cols,Vals]=ConvectionDiffusionSparseMatrix(m,w,options);
    A=sparse(Rows,Cols,Vals,nPoints,nPoints);
    one = ones(nPoints,1);
    disp('Checking that transposed operator vanishes on constants.');
    disp(max(abs(one'*A)));
    
    % Solve
    sol = lsqlin(A,0*one,[],[],one',one'*image);
    disp('Relative error between the solution found and the original image.')
    disp(max(abs( (sol./image)-1)))
end

