% fccAlignment.m
% 
% Perform alignment of image sequence to reference image based on 
%  maximizing image crosscorrelation in Fourier space.
% Same code as fourierCrossCorrelAlignment.m (provided by Marion Silies),
% except only the alignment part.
% 
% INPUT:
%   imageArray - 3D array, corresponding to height x width x nImages
%       sequence of images to align
%   refFrame - single 2D array corresponding to reference to align image to
%   filetype - 'lsm' or 'xml', corresponds to format of imageArray
%
% OUTPUT:
%   registeredImages - 3D array corresponding to aligned image sequence
%
% Updated: 10/3/18 - HHY (int16 output for ScanImage)

function [registeredImages] = fccAlignment(imageArray, refFrame, filetype)
                                                           
    imageArray = double(imageArray);
    refFrame = double(refFrame);
    
    if nargin < 3 || isempty(filetype)
        filetype = 'xml';
    end
    dimension = size(imageArray);
    height = dimension(1);
    width = dimension(2);
    nImages = dimension(3);
    
    % Initialize image arrays:
%     registeredImages = zeros(height, width, nImages, 'uint16'); % aligned
%     unregisteredImages = zeros(height, width, nImages, 'uint16'); % unaligned

    registeredImages = zeros(height, width, nImages,'int16'); % aligned
%     unregisteredImages = zeros(height, width, nImages); % unaligned
    
    % Fourier transform of reference image, DC in (1,1)   [DO NOT FFTSHIFT].
    refImageFT = fft2(refFrame);
    % Upsampling factor (integer). Images will be registered to within 1/usfac 
    % of a pixel. For example usfac = 20 means the images will be registered 
    % within 1/20 of a pixel. (default = 1)
    upSamplingFactor = 1;

    for jImage = 1: nImages
        % Fourier transform of image to register, DC in (1,1) [DO NOT FFTSHIFT].
        switch filetype
            case 'lsm'
                jImageToRegister = imageArray(jImage).data;
                jImageFT = fft2(jImageToRegister);
            case 'xml'
                jImageToRegister = imageArray(:, :, jImage);
                jImageFT = fft2(jImageToRegister);
            otherwise
                error(['Unrecognized file type: ' filetype])        
        end

        [~, registeredImageFT] = dftregistration(refImageFT, jImageFT, ...
                                                 upSamplingFactor); 

        % Convert image to 16 bit format. Check this for use with different bit
        % depths.
%         registeredImages(:, :, jImage) = uint16(ifft2(registeredImageFT));
        registeredImages(:, :, jImage) = int16(ifft2(registeredImageFT));
%         unregisteredImages(:, :, jImage) = jImageToRegister;


        % sometimes, ifft2 inverts image sign 181212 HHY (why?)
        meanValOrig = mean2(jImageToRegister);
        meanValReg = mean2(registeredImages(:,:,jImage));
        
        if (sign(meanValOrig) ~= sign(meanValReg))
            registeredImages(:,:,jImage) = registeredImages(:,:,jImage) * -1;
        end
                                
    end 
end 