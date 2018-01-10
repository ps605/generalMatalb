filenameIn='C:\Users\ps605\Documents\PhD\OpenSim\Media\ FYP_Bushing.avi';  %Movie to be read
filenameOut='C:\Users\ps605\Documents\PhD\OpenSim\Media\FYP_Bushing.avi';  %New filename
desiredDur=10;   %Desire Duration in Seconds


mInfo=VideoReader(filenameIn);
fps=mInfo.NumberOfFrames/desiredDur;

vid = VideoReader(filenameIn);
writerObj = VideoWriter(filenameOut);
writerObj.FrameRate = fps; 
open(writerObj);
nFrames = vid.NumberOfFrames;
vidHeight = vid.Height;
vidWidth = vid.Width;
mov(1:nFrames) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),'colormap', []);
for k = 1 : nFrames
    mov(k).cdata = read(vid, k);
    writeVideo(writerObj, mov(k).cdata);
end
close(writerObj);