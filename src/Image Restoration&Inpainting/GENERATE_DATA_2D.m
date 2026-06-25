%--------------------------------------------------------------------------
% Forward degradation model: given xtr (image) compute b = A*xtr+err
%--------------------------------------------------------------------------

% OUTPUT:
  %   xtr         original uncorrupted image
  %   Axtr A*xtr  blurred image
  %   b           blurred and noisy image

rng(1);  % Control on the random number generator
addpath('GENERATE_DATA_2D/input')

% Set Restoration
  SET_MATRIX_BLUR  = 0; 
  SHOW_CORRUPTIONS = 0;

% Set Inpainting
  SET_MATRIX_INPAINTING = 1;
  SHOW_INPAINTING = 1;

% -------------------------------------------------------------------------
% SELECT THE ORIGINAL IMAGE TO BE TESTED
im_id = 12; % index of the original image to be tested (see the switch-case below)
switch im_id    
    case 0 % cameraman
        im_file     = 'input/cameraman.bmp';
        im_name     = '00_CAMER';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 1 % satellite
        im_file     = 'input/satellite.png';
        im_name     = '01_SATEL';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 2 % checkboard
        im_file     = 'input/checkboard.png';
        im_name     = '02_CHECK';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 3 % checkboard fine
        im_file     = 'input/checkboard_fine.png';
        im_name     = '03_CHECKF';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 4 % barcode_1
        im_file     = 'input/qrcode_1.jpg';
        im_name     = '04_QRCODE1';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 5 % barcode_2
        im_file     = 'input/qrcode_2.png';
        im_name     = '05_QRCODE2';
        f           = 1;
        xtr         = double(imread(im_file))/f;
    case 6 % synthetic binary square
        im_name     = '06_SQUARE';
        f           = 255;
        xtr         = 20 * ones(200,200) / f;
        xtr(60:141,60:141) = 235 / f;    
    case 7 % synthetic binary rectangles
        im_name     = '07_RECTS';
        f           = 255;
        xtr         = 20 * ones(200,200) / f;
        xtr(21:180,18:58)   = 235 / f;
        xtr(41:160,73:103)  = 235 / f;
        xtr(61:140,118:138) = 235 / f;
        xtr(81:120,153:163) = 235 / f;
        xtr(91:110,178:183) = 235 / f;    
    case 8 % synthetic "eggs"
        im_name     = '08_EGGS';
        f       = 255;
        [X,Y]   = meshgrid(0:199,0:199);
        xtr     = (128 + 100 * sin(2*pi*X / 100) .* cos(2*pi*Y / 100) ) / f;
    case 9 % phantom for Xray reconstruction
        target = phantom('Modified Shepp-Logan',512);
        im_name     = '09_SHEPPLOGAN';
        % Choose sparse measurement angles (given in degrees, not radians)
        Nang    = 20;
        angle0  = -90;
        measang = angle0 + [0:(Nang-1)]/Nang*180;
        % Construct measurement (projection or sinogram). 
        f           = 1;
        xtr = radon(target,measang)./f;
        % Construct noisy data
        noiselevel = 0.1; % Choose relative noise level in simulated noisy data
        b = xtr + noiselevel*max(abs(xtr(:)))*randn(size(xtr));
    case 10 % van gogh
        im_file     = 'input/2.png';
        im_name     = '10_VANGOGH';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 11 % parrot
        im_file     = 'input/parrot.png';
        im_name     = '11_PARROT';
        f           = 255;
        xtr         = double(imread(im_file))/f;
    case 12 % inpainted
        im_file     = 'input/reachup1.jpg';
        im_name     = '12_SQUIRRELINP';
        f           = 255;
        xtr         = double(imread(im_file))/f;
        xtr         = rgb2gray(xtr);
end

% extract image dimensions
[h,w] = size(xtr);
N = h * w;

if SET_MATRIX_BLUR
% -------------------------------------------------------------------------
% SET THE BLUR PARAMETERS, THEN GENERATE AND STORE THE BLUR PSF (KERNEL) b_k, 
% FINALLY BLUR THE ORIGINAL IMAGE xtr --> K xtr

    % set the blur parameters then generate and store the blur kernel, b_k
    b_type_id = 3; % index of the blur type (see the switch-case below)
    switch b_type_id
        case 0 % no blur
            b_k             = 1; % kernel
            b_type_descr    = 'NO BLUR';            
        case 1 % average
            b_r             = 5; % radius (pixels)
            b_bc_type       = 0; % boundary conditions type (0->periodic;1->adiabatic)
            b_k             = fspecial('average',(1 + 2 * b_r)); % kernel 
            b_type_descr    = 'AVER';            
        case 2 % disk
            b_r             = 7; % radius (pixels)
            b_bc_type       = 0; % boundary conditions type (0->periodic;1->adiabatic)
            b_k             = fspecial('disk',b_r); % kernel
            b_type_descr    = 'DISK';
        case 3 % Gaussian
            b_r             = 3; % radius (pixels) ...the band is 2 * b_r + 1
            b_s             = 1.5; % standard deviation
            b_bc_type       = 0; % boundary conditions type (0->periodic;1->adiabatic)
            b_k             = fspecial('gaussian',(1 + 2 * b_r),b_s); % kernel
            b_type_descr    = 'GAUSS';
        case 4 % motion
            b_l             = 10; % length (pixels) 
            b_t             = 45; % angle
            b_bc_type       = 0; % boundary conditions type (0->periodic;1->adiabatic)
            b_k             = fspecial('motion',b_l,b_t); % kernel
            b_type_descr    = 'MOTION';
    end
    b_k = b_k / sum(b_k(:)); % normalize the kernel (actually, it is done by fspecial!)
    
    % compute and store the blurred image
    % OTF = psf2otf(PSF) computes the fast Fourier transform (FFT) of the point-spread function (PSF) array
    % blur
        switch b_bc_type
            case 0 % periodic
                K_DFT = psf2otf(b_k,size(xtr));
                Kxtr  = real( ifft2( K_DFT .* fft2(xtr) ) );
            case 1 % adiabatic (Neumann homogeneous)
                % ...              
        end
    else % no blur
        Kxtr  = xtr;
end

if SET_MATRIX_INPAINTING 
% -------------------------------------------------------------------------
% SET THE (SYNTHETIC) INPAINTING MASK PARAMETERS, THEN GENERATE AND STORE  
% THE INPAINTING_MASK IMAGE, FINALLY MASK THE ORIGINAL IMAGE u_true --> Au_true

    % set the masking parameters then generate and store the mask image, M
    mask_type_id = 1; % index of the mask type (see the switch-case below)
    switch mask_type_id
        case 0 % no masking
            M                  = ones(h,w); % mask image
            mask_type_descr    = 'NO'; % type description            
        case 1 % random points
            P_mask             = 0.2; % probability  for a pixel to be masked
            M                  = rand(h,w);
            M(M<=P_mask)       = 0;
            M(M>P_mask)        = 1;
            mask_type_descr    = 'RPS';            
        case 2 % random horizontal strips
            strips_n           = 15; % number of strips
            strips_ht          = 1; % half-thickness of strips (in pixels)
            js                 = round(1 + (h-1) * rand(1,strips_n));
            M                  = ones(h,w);
            for strip_i = 1:strips_n
                j     = js(strip_i);
                j_min = max(j-strips_ht,1);
                j_max = min(j+strips_ht,h);
                M(j_min:j_max,:) = 0;
            end
            mask_type_descr    = 'RHSS'; % type description
        case 3 % random vertical strips
            strips_n           = 15; % number of strips
            strips_ht          = 1; % half-thickness of strips (in pixels)
            is                 = round(1 + (w-1) * rand(1,strips_n));
            M                  = ones(h,w);
            for strip_i = 1:strips_n
                i     = is(strip_i);
                i_min = max(i-strips_ht,1);
                i_max = min(i+strips_ht,w);
                M(:,i_min:i_max,:) = 0;
            end
            mask_type_descr    = 'RVSS'; % type description
        case 4  % load mask from file 
            mask        = double(imread('input/mask_parrot.png'))/255;
            M           = ones(h,w); % mask image
            M           = mask./max(max(mask));
            mask_type_descr    = 'PARROT'; % type description
        case 5 % load mask from file
            mask        = imread('input/reachupmask.gif');
            mask        = double(1 - mask);
            M           = mask./max(max(mask));
            mask_type_descr    = 'SQUIRREL'; % type description
    end
    % mask the original image
    Kxtr = M .* xtr;
end

% -------------------------------------------------------------------------
% SET THE NOISE PARAMETERS, THEN ADD NOISE TO THE BLURRED IMAGE K xtr --> b
n_type_id = 1; % index of the noise type (see the switch-case below)
switch n_type_id
    case 0 % no noise
        n_mean          = 0.0;  % mean (intended for images in [0,255])
        n_sigma         = 0.0;  % stdv (intended for images in [0,255])
        n_type_descr    = 'NO';
        b               = Kxtr;
    case 1 % additive white Gaussian noise (AWGN)
        n_mean          = 0.0;  % mean (intended for images in [0,255])
        n_sigma         = 200;   % stdv (intended for images in [0,255])
        n_type_descr    = 'AWGN';
        n_realiz        = n_mean + n_sigma * randn(h,w);
        b               = Kxtr + n_realiz/255;
    case 2 % additive white uniform noise (AWUN) 
           % (given n_mean and n_sigma -> uniform pdf in [n_mean-sqrt(3)n_sigma,n_mean+sqrt(3)n_sigma])
        n_type_descr    = 'AWUN';
        n_mean          = 0.0;  % mean (intended for images in [0,255])
        n_sigma         = 50;   % stdv (intended for images in [0,255])
        n_hs            = n_sigma * sqrt(3); % half-size of the support        
        n_realiz        = n_mean + n_hs * ( 2 * (rand(h,w) - 0.5) );
        b               = Kxtr + n_realiz/255;
    case 3 % additive white Laplace noise (AWLN)
        n_type_descr    = 'AWLN';
        n_mean          = 0.0;  % mean (intended for images in [0,255])
        n_sigma         = 10;   % stdv (intended for images in [0,255])
        n_sp            = n_sigma / sqrt(2); % scale parameter
        U               = rand(h,w) - 0.5;
        n_realiz        = n_mean - n_sp * ( sign(U) .* log(1 - 2 * abs(U)) );
        b               = Kxtr + n_realiz/255;
    case 4 % multiplicative white Gaussian noise (MWGN)
        n_mean          = 0.0;  % mean (intended for images in [0,255])
        n_sigma         = 20;   % stdv for gray level 127.5 (intended for images in [0,255])
        n_type_descr    = 'MWGN';
        n_realiz        = 1 + ( (n_sigma/127.5) * randn(h,w) );
        b               = Kxtr .* n_realiz;    
    case 7 % impulsive: salt & pepper
        n_type_descr    = 'ISPN';
        n_sigma = 11;
        p   = 0.1;
        ph  = p/2;
        P   = rand(h,w);
        b   = Kxtr;
        b(P < ph)       = 0;
        b(P > (1-ph))   = 1;
    case 8 % impulsive: random-valued
        n_type_descr    = 'IRVN';
        p   = 0.1;
        P   = rand(h,w);
        P1  = rand(h,w);
        b   = Kxtr;
        b(P < p) = P1(P < p);
end

if SHOW_CORRUPTIONS
    figure(100)
    set(gcf,'Position',get(0,'ScreenSize'));

    subplot(2,3,[1,4])
    imshow(uint8(255*xtr));
    title(sprintf('\\textbf{ORIGINAL} ( %d x %d )',h,w),'Interpreter','latex')
    subplot(2,3,2)
    imshow(uint8(255*Kxtr));
    if (SET_MATRIX_BLUR == 0)
        title('\textbf{CORRUPTED} - \textbf{NO BLUR}','Interpreter','latex');
    else
        title(sprintf('\\textbf{CORRUPTED} - %s BLUR',b_type_descr),'interpreter','latex');
    end
    subplot(2,3,3)
    imshow(uint8(255*b));
    if (SET_MATRIX_BLUR == 0)
        if (n_type_id == 0)
            title('\textbf{CORRUPTED} - \textbf{NO BLUR and NO NOISE}','Interpreter','latex');
        else
            title(sprintf('\\textbf{CORRUPTED} - NO BLUR and %s NOISE',n_type_descr),'Interpreter','latex');
        end
    else
        if (n_type_id == 0)
            title(sprintf('\\textbf{CORRUPTED} - %s BLUR and NO NOISE',b_type_descr),'Interpreter','latex');
        else
            title(sprintf('\\textbf{CORRUPTED} - %s BLUR and %s NOISE',b_type_descr,n_type_descr),'Interpreter','latex');
        end
        subplot(2,3,5)
        [b_k_h,b_k_w] = size(b_k);
        enl     = 1;
        b_k_2   = zeros(b_k_h+2*enl,b_k_w+2*enl);
        b_k_2((1+enl):(end-enl),(1+enl):(end-enl)) = b_k;
        image(255*b_k_2/max(b_k_2(:)));
        axis equal;
        axis tight;
        title(sprintf('\\textbf{BLUR KERNEL (PSF)}: %s %dx%d',b_type_descr,size(b_k,1),size(b_k,2)),'Interpreter','latex');
    end
    subplot(2,3,6)
    if (n_type_id == 0)
        imshow(uint8(127.5*ones(size(xtr))));
    elseif (n_type_id < 4)
        imshow(uint8(127.5+n_realiz));
    elseif (n_type_id < 7)
        imshow(uint8(127.5+255*(b-xtr)));
    elseif (n_type_id == 7)
        Omega_0             = zeros(size(xtr));
        Omega_0(P < ph)     = 1;
        Omega_0(P > (1-ph)) = 1;
        tmp = 127.5*ones(size(xtr));
        tmp(Omega_0>0) = 255*b(Omega_0>0);
        imshow(uint8(tmp));        
    elseif (n_type_id == 8)
        Omega_0         = zeros(size(xtr));
        Omega_0(P < p)  = 1;
        tmp = 127.5*ones(size(xtr));
        tmp(Omega_0>0) = 255*b(Omega_0>0);
        imshow(uint8(tmp)); 
    end
    if (n_type_id == 0)
        title('\textbf{NO NOISE}','interpreter','latex');
    else
        title(sprintf('%s NOISE',n_type_descr),'Interpreter','latex');
    end
end
    
if SHOW_INPAINTING
% reset to zero the inpainting mask pixels
    b(M==0) = 1;
% choose the color for visualization of inpainting mask
    M_COL = [255;0;0];    %not used
    figure(100)
    set(gcf,'Position',get(0,'ScreenSize'));
    subplot(2,2,1)
    imshow(uint8(255*xtr));
    title(sprintf('\\textbf{ORIGINAL} ( %d x %d )',h,w),'Interpreter','latex');
    subplot(2,2,2)
    imshow(uint8(255*b));
    title(sprintf('\\textbf{CORRUPTED by} %s \\textbf{MASK} and %s \\textbf{NOISE}',mask_type_descr,n_type_descr),'Interpreter','latex');
    subplot(2,2,3)
    imshow(uint8(255*M));
    title(sprintf('%s \\textbf{BINARY INPAINTING MASK}',mask_type_descr),'interpreter','latex');

    subplot(2,2,4)
    if (n_type_id == 0)
        imshow(uint8(127.5*ones(size(xtr))));
    elseif (n_type_id < 4)
        imshow(uint8(127.5+n_realiz));
    elseif (n_type_id < 7)
        imshow(uint8(127.5+255*(b-xtr)));
    elseif (n_type_id == 7)
        Omega_0             = zeros(size(xtr));
        Omega_0(P < ph)     = 1;
        Omega_0(P > (1-ph)) = 1;
        tmp = 127.5*ones(size(xtr));
        tmp(Omega_0>0) = 255*b(Omega_0>0);
        imshow(uint8(tmp));        
    elseif (n_type_id == 8)
        Omega_0         = zeros(size(xtr));
        Omega_0(P < p)  = 1;
        tmp = 127.5*ones(size(xtr));
        tmp(Omega_0>0) = 255*b(Omega_0>0);
        imshow(uint8(tmp)); 
    end
    if (n_type_id == 0)
        title('\textbf{NO NOISE}','interpreter','latex');
    else
        title(sprintf('%s NOISE',n_type_descr),'Interpreter','latex');
    end
 end



