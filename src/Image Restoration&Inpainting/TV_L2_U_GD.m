function [xest,itrs,r] = TV_L2_U_GD(b,b_k,lambda,eps,x0,tau,itrs_max,ch_th,debug_cw,xtr)

% OUTPUT:
  %
  %

% INPUT:
  % b: degradated image
  % b_k: kernel type
  % lambda: regularization parameter
  % eps: Parameter used for the differentiability of the regularization
  %       term
  % x0: initial condition
  % tau : descent step for the Gradient Descent
  % itrs_max: maximum number of iterations
  % ch_th: change treshold for the residual error
  % xtr: original noisy free image

% OUTPUT: 
  % xest : Reconstrructed image using Total Variation
  % itrs : Number of iterations used


%

% Let's create the blur matrix using the kernel b_k
  A_DFT = psf2otf(b_k,size(b));
  AT_DFT = conj(A_DFT);

  Ab_DFT = A_DFT .* fft2(b);

% Let's create the matrix of the second order finite differences
  Dh_DFT          = psf2otf([1,-1],size(b));
  Dv_DFT          = psf2otf([1;-1],size(b));
  DhT_DFT         = conj(Dh_DFT);
  DvT_DFT         = conj(Dv_DFT);


% First step of the gradient descent method fot the Total Variation
  xest = x0;
  r = norm ( xtr(:) - xest(:) ) / norm( xtr(:) );
  itrs = 1;

while( (r >= ch_th) && (itrs <= itrs_max) )
    %
    x0 = xest;

    TV_x = real( ifft2( Dh_DFT .* fft2(x0) ) );
    TV_y = real( ifft2( Dv_DFT .* fft2(x0) ) );
    TV_term = sqrt( TV_x.^2 + TV_y.^2 + eps^2 );

    F_1 = lambda * div( gradx(x0)./TV_term, grady(x0)./TV_term );
    %F_2 = real( ifft2( A_DFT .* fft2(x0) - fft2(b) ) );
    F_2 = x0 - b;

    grad_F = F_1 + F_2;

    xest = x0 + tau * grad_F;

    r = norm ( xtr(:) - xest(:) ) / norm( xtr(:) );
   
    itrs = itrs + 1;

end



