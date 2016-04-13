clear;
clc;
close all;

fpi=fopen('overlap_dnb');
%hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
%field=textscan(hline{1}{1},'%s');
clear format;
format='%s';
for i=2:70
    format=[format,' %f'];
end
lines =textscan(fpi, format,1000000,'delimiter', '\t');
egf_ipi=lines{1};
profile = [];
for i = 2 :70
    profile = [profile, lines{i}];
end
fclose(fpi);

psize=size(profile);
control=profile(:,1);
EGF=profile(:,2:35);
HRG=profile(:,36:69);
egf_control=profile(:,2:5);
hrg_control=profile(:,[1,36:37]);

for i=1:psize(1)
    EGF(i,:)=(EGF(i,:)-mean(egf_control(i,:)))/std(egf_control(i,:));
    HRG(i,:)=(HRG(i,:)-mean(hrg_control(i,:)))/std(hrg_control(i,:));
end

for i=1:16
    egf_case(:,i,:)=EGF(:,2*i-1:2*i+2);
    hrg_case(:,i,:)=HRG(:,2*i-1:2*i+2);
end
%***********************************************************
fpi=fopen('sample_network');
hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
field=textscan(hline{1}{1},'%s');
clear format;
format='%s';
for i=2:70
    format=[format,' %f'];
end
lines =textscan(fpi, format,1000000,'delimiter', '\t');
egf_ipi=lines{1};
profile = [];
for i = 2 :70
    profile = [profile, lines{i}];
end
fclose(fpi);

rand_psize=size(profile);
rand_EGF=profile(:,2:35);
rand_HRG=profile(:,36:69);
rand_egf_control=profile(:,2:5);
rand_hrg_control=profile(:,[1,36:37]);

for i=1:rand_psize(1)
    rand_EGF(i,:)=(rand_EGF(i,:)-mean(rand_egf_control(i,:)))/std(rand_egf_control(i,:));
    rand_HRG(i,:)=(rand_HRG(i,:)-mean(rand_hrg_control(i,:)))/std(rand_hrg_control(i,:));
end

for i=1:16
    rand_egf_case(:,i,:)=rand_EGF(:,2*i-1:2*i+2);
    rand_hrg_case(:,i,:)=rand_HRG(:,2*i-1:2*i+2);
end

stdRNA=zeros(1,16);
pccRNA=zeros(1,16);
outpcc=zeros(1,16);
for j=1:16
    for i=1:psize(1)
        stdRNA(j)=stdRNA(j)+std(egf_case(i,j,:));
    end
    tmp_case=reshape(egf_case(:,j,:),psize(1),4);
    pccall=corr(tmp_case');
    pccRNA(j)=(sum(sum(abs(pccall)))-psize(1))/(psize(1)*(psize(1)-1));
end
stdRNA=stdRNA/psize(1);

for i=1:psize(1)
    for j=1:rand_psize(1)
        for k=1:16
            tmp1=reshape(egf_case(i,k,:),1,4);
            tmp2=reshape(rand_egf_case(j,k,:),1,4);
            outpcc(k)=outpcc(k)+abs(corr(tmp1',tmp2'))/psize(1)/rand_psize(1);
        end
    end
end

egf_compositindex=zeros(16,1);
for s=1:16
    egf_compositindex(s)=stdRNA(s)*pccRNA(s)/outpcc(s);
end
%**************************************************************************
stdRNA=zeros(1,16);
pccRNA=zeros(1,16);
outpcc=zeros(1,16);
for j=1:16
    for i=1:psize(1)
        stdRNA(j)=stdRNA(j)+std(hrg_case(i,j,:));
    end
    tmp_case=reshape(hrg_case(:,j,:),psize(1),4);
    pccall=corr(tmp_case');
    pccRNA(j)=(sum(sum(abs(pccall)))-psize(1))/(psize(1)*(psize(1)-1));
end
stdRNA=stdRNA/psize(1);

for i=1:psize(1)
    for j=1:rand_psize(1)
        for k=1:16
            tmp1=reshape(hrg_case(i,k,:),1,4);
            tmp2=reshape(rand_hrg_case(j,k,:),1,4);
            outpcc(k)=outpcc(k)+abs(corr(tmp1',tmp2'))/psize(1)/rand_psize(1);
        end
    end
end

hrg_compositindex=zeros(16,1);
for s=1:16
    hrg_compositindex(s)=stdRNA(s)*pccRNA(s)/outpcc(s);
end

t=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
plot(t,egf_compositindex);
hold on;
plot(t,hrg_compositindex,'r');
