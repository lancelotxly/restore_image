function s = snr(noisydata, original)

noisydata   =   double(noisydata);
original    =   double(original);


tmp           = original;
var_original  = sum(sum(tmp.*tmp));

noise      = noisydata - original;

tmp        = noise ;
var_noise  = sum(sum(tmp.*tmp));

if var_noise == 0
    s = 999.99; %% INF. clean image
else
    s = 10 * log10(var_original / var_noise);
end
return