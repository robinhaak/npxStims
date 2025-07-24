function gratingFrames = generateOrientedGrating(angle_deg, patchSize, pixelsPerCycle, temporalFreq_Hz, fps, nFrames, background, white)
% generateOrientedGrating creates a sequence of drifting grating frames
% inside a fixed-size square patch with specified orientation.

[X, Y] = meshgrid(1:patchSize, 1:patchSize);
X = X - mean(X(:));  % center coordinates
Y = Y - mean(Y(:));

theta = deg2rad(angle_deg);
Xrot = cos(theta) * X + sin(theta) * Y;

spatialFreq = 2 * pi / pixelsPerCycle;
phaseShiftPerFrame = 2 * pi * temporalFreq_Hz / fps;

gratingFrames = zeros(patchSize, patchSize, nFrames);

for f = 1:nFrames
    phase = (f - 1) * phaseShiftPerFrame;
    frame = cos(spatialFreq * Xrot + phase);
    frame = background + ((white - background)/2) * frame;
    gratingFrames(:,:,f) = frame;
end
end
