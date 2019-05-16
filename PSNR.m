function s = psnr(noisydata, original)

noisydata=double(noisydata);
original=double(original);

[m,n] = size(noisydata);

peak=255*255*m*n;

noise  = noisydata - original;
nostotal = sum(sum(noise.*noise));

if nostotal == 0
    s = 999.99; %% INF. clean image
else
    s = 10 * log10(peak./nostotal);
end

return