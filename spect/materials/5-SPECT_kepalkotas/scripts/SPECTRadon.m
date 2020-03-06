function [RT,xp] = SPECTRadon (I, theta)
%
% [RT, xp] = SPECTRadon(I, theta)
%               Calculates the 2D-Radon transform of the matrix I at angles given in
%               theta. To each element of theta corresponds a column in RT.
%               The variable xp represents the x-axis of the rotated coordinate.
%               If theta is not defined, then 0:179 is assumed.
%

  if (nargin == 0 || nargin > 2 || size(size(I),2)~=2)
    print_usage ();
    return;
  elseif (nargin == 1)
    theta = 0:179;
  end

  [m,n] = size (I);

  % center of image

  xc = floor ((m+1)/2);
  yc = floor ((n+1)/2);

  % divide each pixel into 2x2 subpixels

  d = reshape (I,[1 m 1 n]);
  d = d([1 1],:,[1 1],:);
  d = reshape (d,[2*m 2*n])/4;

  b = ceil (sqrt (sum (size (I).^2))/2 + 1);
  xp = [-b:b]';
  sz = size(xp);

  [X,Y] = ndgrid (0.75 - xc + [0:2*m-1]/2,0.75 - yc + [0:2*n-1]/2);

  X = X(:)';
  Y = Y(:)';
  d = d(:)';

  th = theta*pi/180;

  for l=1:length (theta)
    % project each pixel to vector (-sin(th),cos(th))
    Xp = -sin (th(l)) * X + cos (th(l)) * Y;
   
    ip = Xp + b + 1;

    k = floor (ip);
    frac = ip-k;

    RT(:,l) = accumarray (k',d .* (1-frac),sz) + accumarray (k'+1,d .* frac,sz);
  end

end



function print_usage()
fprintf( [ 
 '[RT, xp] = SPECTRadon(I, theta)\n' ...
 '              Calculates the 2D-Radon transform of the matrix I at angles given in\n'  ...
 '              theta. To each element of theta corresponds a column in RT.\n'  ...
 '              The variable xp represents the x-axis of the rotated coordinate.\n'  ...
 '              If theta is not defined, then 0:179 is assumed.\n\n'  ...
 ]);
end



%!test
%! A = radon(ones(2,2),30);
%! assert (A,[0 0 0.608253175473055 2.103325780167649 1.236538105676658 0.051882938682637 0]',1e-10) 
