fn=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\mesh.txt','r');
node=fscanf(fn,'%f',[1,inf]);%Nf,Ng,Ns,Ny,M,N,Nt
fclose(fn);
Nf=node(1);Ng=node(2);Ns=node(3);Ny=node(4);
M=Ny+1;
N=Nf+Ng+Ns+1;

fid=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\temperature.txt','r');
%T=textread('E:\documents\doctors degree\software\tansistant\output\temperature.txt');
T=fscanf(fid,'%f',[M-1,N]);
fclose(fid);

fctime=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\ctime.txt','r');
%T=textread('E:\documents\doctors degree\software\tansistant\output\temperature.txt');
ctime=fscanf(fctime,'%f',[1,inf]);
fclose(fctime);

ftout=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tout.txt','r');
%T=textread('E:\documents\doctors degree\software\tansistant\output\temperature.txt');
tout=fscanf(ftout,'%f',[1,inf]);
fclose(ftout);

ftpow=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tpow.txt','r');
%T=textread('E:\documents\doctors degree\software\tansistant\output\temperature.txt');
tpow=fscanf(ftpow,'%f',[1,inf]);
fclose(ftpow);

ftuin=fopen('E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tuin.txt','r');
%T=textread('E:\documents\doctors degree\software\tansistant\output\temperature.txt');
tuin=fscanf(ftuin,'%f',[1,inf]);
fclose(ftuin);
Nt=length(ctime)
y=1:1:Nt
subplot(2,2,1);plot(y,tpow),title('功率随时间分布'),xlabel('时间'),ylabel('功率')
subplot(2,2,2);plot(y,tuin),title('入口流速随时间分布'),xlabel('时间'),ylabel('入口流速')
subplot(2,2,3);plot(y,tout),title('堆芯出口温度随时间分布'),xlabel('时间'),ylabel('堆芯出口温度')
