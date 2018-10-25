function y=convert_back_to_time(fyo, Nsamples, freq_indtest)

% Nsamples=109;yo=demean(randn(Nsamples,1));
% NumUniquePts=ceil((Nsamples+1)/2); freq_indtest=2:25; fyo=fft(yo); fyo=fyo(freq_indtest); 
% y=convert_back_to_time(fyo,Nsamples, freq_indtest);
% sfigure;plot(y);ho;plot(yo,'r--')

NumUniquePts=ceil((Nsamples+1)/2);

if(freq_indtest==-1),
    freq_indtest=1:NumUniquePts;
end;

full_fyo=zeros(Nsamples,1);

full_fyo(freq_indtest)=fyo;  


tmp=zeros(1,length(full_fyo));

tmp(1:length(full_fyo))=full_fyo; 
tmp(length(NumUniquePts)+1:end)=tmp(length(NumUniquePts)+1:end)+full_fyo(end:-1:2)';
 
y=real(ifft(tmp,'symmetric'));
