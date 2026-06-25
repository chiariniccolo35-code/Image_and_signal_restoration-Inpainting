function [rel_err] = compute_rel_err(xtr,x)

% INPUT:
  % xtr is the original noise and blur free image
  % x is the reconstructed image

rel_err = norm( xtr(:) - x(:) ) / norm( xtr(:) );

end