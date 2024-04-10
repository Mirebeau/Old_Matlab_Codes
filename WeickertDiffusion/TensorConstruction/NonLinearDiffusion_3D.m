% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function smoothed = NonLinearDiffusion_3D(image_to_denoise, options)
  time = 0.;
  smoothed = image_to_denoise;
  s = size(smoothed); s=s(1:3);
  layers=size(smoothed,4);

  for iter = 1:GetOptions(options,'max_diff_iter',100)  
    tensor = StructureTensor_3D(smoothed,options);
    tensor = WeickertTensor_3D(tensor,options);
    A = DiffusionSparseMatrix_3D(tensor);
    %maximum stable time step (matrix A has positive diagonal entries)
    time_delta = 1./full(max(A(:)));
    time_delta = time_delta*GetOptions(options,'time_delta_safety',0.7);
    
    smoothed = reshape(smoothed,[prod(s),layers]);
    for niter=1:GetOptions(options,'max_inner_iter',10)
        dt = min(time_delta, GetOptions(options,'final_time') - time);
        time = time+dt;
        fprintf('Time step : %f, leading to %f. Iteration %d\n',dt,time,iter);
        smoothed=smoothed - dt*A*smoothed;
    end
    smoothed = reshape(smoothed,[s,layers]);
    
    if dt < time_delta
      return;
    end
  end
end