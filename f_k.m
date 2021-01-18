clear
close
clc 
%% 1// Organizing Data 

Nsigma=3; % Nombre de sigma qu'on va prendre en compte comme signal 
load('data_bin_line_selected.mat') % all stations !
% load('data_bin_line_I.mat')
% load('all_2-60sec.mat')
[dist, iA,~] = unique(dist);

data = data(:,iA);

[~,coldata]=find(data~=0);
coldata=unique(coldata);

data1=data(:,coldata);
x=dist(coldata);
U=data1/max(max(data1));
freqlimits=[0 0.5];

fs=1/mean(diff(time)); % sampling frequency Hz
min_x = min(x); % min offset [m] of reciever spread
max_x = max(x); % max offset [m] of reciever spread

d_x = mean(diff(x)); % distance between receivers [m]
x = min_x:d_x:max_x; % spatial sampling vector (receiver position)
t=1/fs:1/fs:length(U(:,1))/fs;  % time vector in seconds

L=length(t); % number os samples
T=max(t); % max time
if ~mod(L,2)
    f = (1/T)*(0:L/2-1); % n is even
else
    f = (1/T)*(0:(L-1)/2); % n is odd
end

%% 2// On applique la SVD

Ns = 7; % Arbitrary chosen 
aNr_U2=size(U,2)-11;
dtemp=zeros(Ns,aNr_U2);
for iNs = 1 : Ns
    dtemp(iNs,:) = randperm(size(U,2),aNr_U2);
end
U2 = zeros(Ns,size(U,1),aNr_U2);
x2 = zeros(Ns,aNr_U2);
for iNs = 1:size(dtemp,1)
    dtemp(iNs,:) = sort(dtemp(iNs,:));
    U2(iNs,:,:) = U(:,dtemp(iNs,:));
    x2(iNs,:) = x(dtemp(iNs,:));
end
aNr_U2=round(aNr_U2*1.5);
x2_int=zeros(Ns,aNr_U2);
for iNs=1:Ns
maxi = max(x2(iNs,:));
mini = min(x2(iNs,:));
x2_int(iNs,:) = linspace(mini,maxi,aNr_U2);
end
U2_int = zeros(Ns,size(U,1),aNr_U2);
for iNs = 1 : Ns
    U2_int(iNs,:,:) = interp2(time',x2(iNs,:),squeeze(U2(iNs,:,:))',time',x2_int(iNs,:))';
end
U2_int_fft=zeros(Ns,size(U,1),aNr_U2); % Rq: N=length(time)
for iNs=1:Ns
U2_int_fft(iNs,:,:)=fft(squeeze(U2_int(iNs,:,:)));
end
U2_int_fft=U2_int_fft(:,1:length(f),:);
U2_int_fft=permute(U2_int_fft,[3,1,2]);

U_singul = zeros(aNr_U2,Ns,length(f));
V_singul = zeros(Ns,Ns,length(f));
S = zeros(Ns,Ns,length(f));
sigma = zeros(Ns,length(f));

for Nf = 1 : length(f)
    [U_singul(:,:,Nf),S(:,:,Nf),V_singul(:,:,Nf)] = svd(squeeze(U2_int_fft(:,:,Nf)),'econ');
    sigma(:,Nf) = diag( squeeze(S(:,:,Nf)) );
end
figure(2);clf;
plot(f,20*log10(sigma/max(max(abs(sigma)))))
xlabel('frequency (Hz)'); ylabel('Singular values (in dB of peak value)')
xlim([0 max(f)])
[fcut,sigmacutsup] = ginput(14);
fcut = [f(1),fcut.',f(end)];
sigmacutsup = [sigmacutsup(1),sigmacutsup',sigmacutsup(end)];
sigmacutsup = interp1(fcut,sigmacutsup,f);
hold on; plot(f,sigmacutsup,'k--'); drawnow;

for Ns = 1:Ns
    titi = 20*log10(sigma(Ns,:)/max(max(abs(sigma))));
    [~,idfx]= find(titi < sigmacutsup);
    U_singul(:,Ns,idfx) = 0;
    V_singul(Ns,:,idfx) = 0;
    sigma(Ns,idfx) = 0;
end
U2_int_fft_filtered = zeros(size(U2_int_fft));
V = zeros(Ns,Ns);
U = zeros(aNr_U2,aNr_U2);
S = zeros(aNr_U2,Ns);
for Nf = 1 : length(f)
    U(1:aNr_U2,1:Ns) = U_singul(:,:,Nf);
    V(1:Ns,1:Ns) = V_singul(:,:,Nf);
         for iNs = 1:Ns
             S(iNs,iNs) = sigma(iNs,Nf);
         end
    U2_int_fft_filtered(:,:,Nf) = U*S*V';
end

%% // (Vecteur k)
% nombre d'onde 

D=size(x2_int,2);

if ~mod(D,2)
    k=2*pi*(1/max(max(x2_int)))*(0:D/2-1);
else
    k=2*pi*(1/max(max(x2_int)))*(0:(D-1)/2);
end


%% K 2
% k=zeros(Ns,round(size(x2_int,2)/2));
% xe=zeros(Ns,1);
% dk=zeros(Ns,1);
% ke=zeros(Ns,1);
% kmax=zeros(Ns,1);
% kNyq=zeros(Ns,1);
% D=zeros(Ns,1);  % idk1=zeros(Ns,1);
% for ff=1:Ns  
% xe(ff,1) = mean(diff( x2_int(ff,:)));
% dk(ff,1)=2*pi/(max( x2_int(ff,:)));
% ke(ff,1) = 2*pi/xe(ff,1);
% kmax(ff,1)=2*pi/xe(ff,1);
% kNyq(ff,1)=ke(ff,1)/2;
% D(ff,1)=length( x2_int(ff,:));
% k(ff,:)=linspace(0,(ke(ff,1)/2),round(size(x2_int,2)/2)); % k = 0:dk:dk*(D-1);
% end

%% 3// Calculate the F_K image 

c=k; % j'ai juste remplacé c dans la boucle par k 
     % w./c a été remplacé par c=k(wavenumber)
[~,fcol]=find(squeeze(U_singul(:,1,:))~=0);
fcol=(unique(fcol))';
f=f(fcol);
Norm=zeros(Nsigma,length(c),length(fcol));
                        %c=k
Uu=U_singul(:,:,fcol);
for iNs=1:Nsigma
% aicounter=iNs;
U=squeeze(Uu(:,iNs,:)); 
[~,fcol1]=find(U~=0);
fcol1=unique(fcol1);
if ~isempty(fcol1)
Rnorm=U(:,fcol1)./abs(U(:,fcol1))./sqrt(size(x2_int,2));
w=f(fcol1)*2*pi; % convert the freq vector to rads/sec
AbsSumTestc=zeros(length(c),length(w));
                        %c=k
for iw=1:length(w)
           % exp(j*x*k)           % c=k                           % c=k
        test=exp(1j*x2_int(iNs,:)'.*c)./norm(exp(1j*x2_int(iNs,:)'.*c));
        AbsSumTestc(:,iw)=(abs(sum(real(test).*real(Rnorm(:,iw)./norm((Rnorm(:,iw))))+1i*(imag(test).*imag(Rnorm(:,iw)./norm((Rnorm(:,iw)))))))).^2;        
end
end
Norm(iNs,:,fcol1)=AbsSumTestc;
end

Norm_sum=zeros(size(Norm,2),size(Norm,3));
for iNs=1:Nsigma
    Norm_sum=squeeze(Norm(iNs,:,:))+Norm_sum;
end

a=jet ;   
colormap(a)
colormapeditor
figure(3);clf;
          %c=k
h=pcolor(f,c,Norm_sum./max(Norm_sum)); 
set(h, 'EdgeColor', 'none');
xlim(freqlimits)
ylim([min(c) max(c)])
colormap (a);
caxis('auto');
xlabel('Freq 0-0.5(Hz)');
ylabel('k : wave number  (1/km)');
gca;
title('F-K')
set(gca,'FontSize',15,'fontweight','bold') 
drawnow
