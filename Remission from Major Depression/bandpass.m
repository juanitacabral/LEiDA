function x=bandpass(y,freqrange,fres,do_plot)

% function x=bandpass(y,lp,hp)
%
% does bandpass filtering via FFT
% 0<lp<1 proportion of high freq to cut out
% 0<hp<1 proportion of low freq to cut out

if(nargin<4)
  do_plot=0;
end;
 
Nsamples=length(y);

Nunique_points=ceil((Nsamples+1)/2);

fHz = (0:Nunique_points-1)*fres/Nsamples;

freq_ind=intersect(find(fHz>=freqrange(1)),find(fHz<=freqrange(2)));

fy=fft(y);
fyo=fy(freq_ind);

x=convert_back_to_time(fyo, Nsamples, freq_ind);

if(do_plot),   
  fy2=zeros(size(fy));
  fy2(freq_ind)=fyo;
  figure;plot(fHz,abs(fy2(1:length(fHz))));
  figure;plot(x);ho;plot(y,'r--');
end;