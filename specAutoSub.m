function [specSub,kSol] = specAutoSub(spcIn,spcRef)
%Autosubstract thereference spectrum from the input spectrum to minimise
%spectrum lenght
%the spectrum length minimizing (SLM) reported by Kozlov from Novosibirsk (Kozlov, D.; Besov, A.; Applied spectroscopy 2011, 65, 918-923.). 

if nargin ~=2
    disp('Insuficient Arguments')
    disp('[specSub,kSol] = specAutoSub(spcIn,spcRef)')
    return
end

if isvector(spcRef)==0 || isvector(spcIn)==0 
   disp('No valid spectra selected')
   disp('[specSub,kSol] = specAutoSub(spcIn,spcRef)')
   return
end

if size(spcIn,2)>size(spcIn,1)
    spcIn = transpose(spcIn);
end

if size(spcRef,2)>size(spcRef,1)
    spcRef = transpose(spcRef);
end

if length(spcRef) ~= length(spcIn)
   disp('No comaptible spectra')
end

kMax = max([1/(max(spcIn)/ max(spcRef)), -1/(min(spcIn)/ min(spcRef))]);
kMin = min([1/(max(spcIn)/ max(spcRef)), -1/(min(spcIn)/ min(spcRef))]);

%tic %Start timer
optfun = @(k) sum(abs((spcIn(1:end-1)-(k*spcRef(1:end-1)))-(spcIn(2:end)-(k*spcRef(2:end)))));

kSol = fminbnd(optfun,kMin,kMax);
if isempty(kSol)
      disp('Optmisation failed');
      return
end %if
specSub  = spcIn - (kSol.*spcRef);%Calculates the resulting spectrum
%disp(['Caculated in ',num2str(toc*1000),' ms']);
end %function