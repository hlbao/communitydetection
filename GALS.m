%This code is designed by Honglin Bao @ Banzhaf Lab, CSE Michigan State University and NSF BEACON Center for the Study of Evolution in Action%
%Version 1, April 2019%
%Community detection problem in heterogenous and evolving networks through clustering and genetic algorithm.%
%Contact: baohongl@msu.edu%
%https://cse.msu.edu/person/3%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function GALS %%%%%%%%%%%%%%%%%%%%%%%% function Q as fitness evaluation%%%%
clear,clc;
close all;
global str;
str='football';%%%% Test data set%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global labels;
global Q;
global tt;
[labels,Q,tt]=GACD(str);




function [labels,Q,tt]=GACD(str)
%%%%%%%%%%Global Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global W;
global degrees;
global m;
global jin;
global pop;
global pop_Q;
global pop2;
global pop2_Q;
global p1;
global pop2_size;
%%%%%%%%%%Global Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%Variable Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop2_size=0;
pop_size=80;%%% 80
T=200;%%% 200
p1=0.8;
W= load(str);
W = W.ijv;
count=W(size(W,1),2);
pop=zeros(pop_size,count);
pop_Q=zeros(1,pop_size);
pop2=zeros(pop_size,count);
pop2_Q=zeros(1,pop_size);
temp=diff(W(:,2)');
jin=find([1,temp,1]);
degrees=zeros(1,count);
for i=1:length(jin)-1
    v=W(jin(i):jin(i+1)-1,3);
    degrees(i)=sum(v);
end
m=sum(degrees)/2;
%%%%%%%%%%Variable Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%Population Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels=zeros(1,count);
for i=1:pop_size
    for j=1:count
        v=W(jin(j):jin(j+1)-1,1);
        ix=ran(length(v));
        labels(j)=v(ix);
    end
    pop(i,:)=labels(:);
    labels2=convert_labels(labels);
    pop_Q(i)=compute_Q2(labels2);
end
%%%%%%%%%%Population Initialization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%进化
x=[1;T];
y=zeros(1,T);
tic;
for t=1:T
    Crossover;
    Mutation;
    Select;
    y(t)=pop_Q(1);
end
tt=toc;
plot(x,y);
%%%%%%%%%%%%进化
labels=pop(1,:);
labels=convert_labels(labels);
Q=pop_Q(1);


function Crossover
global pop;
global pop2;
global p1;
global pop2_size;
global t2;
n=size(pop,1);
n2=n/2;
m=size(pop,2);
count=0;
t1=rand(1,n);
[t2,ix]=sort(t1);
for i=1:n2
    pp=rand();
    if pp<p1
        r1=randi(1,m,[0,1]);
        r2=-r1+1;
        i1=i*2-1;
        i2=i*2;
        v1=pop(ix(i1),:);
        v2=pop(ix(i2),:);
        count=count+1;
        pop2(count,:)=v1.*r1+v2.*r2;
        count=count+1;
        pop2(count,:)=v1.*r2+v2.*r1;
    end
end
pop2_size=count;

function Mutation
global pop2_size;
for i=1:pop2_size
    one_scan(i);
end

function Select
global pop;
global pop_Q;
global pop2;
global pop2_Q;
global t1;
t3=[pop_Q,pop2_Q];
[t1,ix]=sort(t3,'descend');
t2=[pop;pop2];
pop=t2(ix(1:size(pop,1)),:);
pop_Q=t3(ix(1:length(pop_Q)));
pop2(:,:)=0;
pop2_Q(:)=0;

function Q=compute_Q2(labels)
global W;
global degrees;
global m;
global jin;
Q=0;
[indexs,di]=find_communities(labels);
V=zeros(1,length(labels));
for i=1:length(di)-1
    index=indexs(di(i):di(i+1)-1);
    for j=1:length(index)
        id=index(j);
        nn2=W(jin(id):jin(id+1)-1,[1,3])';
        V(nn2(1,:))=nn2(2,:);
        left=sum(V(index));
        right=sum(degrees(id)*degrees(index))/(2*m);
        val=left-right;
        Q=Q+val;
        V(nn2(1,:))=0;
    end
end
Q=Q/(2*m);


function one_scan(num) 
global W;
global degrees;
global m;
global jin;
global pop2;
global pop2_Q;
global temp;
global co;
o_labels=pop2(num,:);
vs=count_vs(o_labels);
labels=convert_labels(o_labels);
count2=length(labels);
V=zeros(1,count2);
[indexs,di]=find_communities(labels);
new_labels=zeros(1,count2);
for i=1:length(di)-1
    new_labels(indexs(di(i):di(i+1)-1))=i;
end
labels=new_labels;
ff=rand(1,count2);
[temp,gg]=sort(ff);
clear ff temp;
for rr=1:length(labels)
    i=gg(rr);
    if vs(i)==0
        %%%%%
%         tp=0.5;
%         if rand()>tp
%             continue;
%         end
        %%%%%
        nn2=W(jin(i):jin(i+1)-1,[1,3])';
        nn=nn2(1,:);
        colors=unique2(labels(nn));
        ma=-realmax;
        co=-1;
        du=degrees(i);
        V(nn2(1,:))=nn2(2,:);
        for j=1:length(colors)
            color=colors(j);
            index2=indexs(di(color):di(color+1)-1);
            if color~=labels(i)
                index=[i,index2];
            else
                index=index2;
            end
            left=sum(V(index));
            sub_d=degrees(index);
            right=sum(du*sub_d)/(2*m);
            val=left-right;
            if val==ma
                co=color;
            end
            if val>ma
                ma=val;
                co=color;
            end
        end
        V(nn2(1,:))=0;
        if length(co)==1
            labels(i)=co;
        else
            l=length(co);
            id=ran(l);
            labels(i)=co(id);
        end
        ixs=find(labels(nn)==labels(i));
        id=ran(length(ixs));
        ix=ixs(id);
        vs(o_labels(i))=vs(o_labels(i))-1;
        o_labels(i)=nn(ix);
        vs(o_labels(i))=vs(o_labels(i))+1;
    end
end
pop2_Q(num)=compute_Q2(labels);
pop2(num,:)=o_labels(:);

function [indexs,di]=find_communities(labels)
[v,indexs]=sort(labels);
temp=diff(v);
di=find([1,temp,1]);

function id=ran(l)
id = ceil(l.*rand());

function y=unique2(x)
x = sort(x);  
difference = diff([x,NaN]);  
y = x(difference~=0);

function labels=convert_labels(g)
global sp;
n=length(g);
x=[1;n];
W=sparse([x,g],[g,x],ones(1,2*n));
[r,c]=find(W);
r=r';
c=c';
temp=diff(c);
jin=find([1,temp,1]);

count=1;
labels=zeros(1,n);
v=zeros(1,n);

s=zeros(1,n);
sp=1;

for i=1:n
    if v(i)==1
        continue;
    end
    sp=1;
    s(sp)=i;
    sp=sp+1;
    while sp~=1
        sp=sp-1;
        t=s(sp);
        labels(t)=count;
        v(t)=1;
        for j=jin(t):jin(t+1)-1
            k=r(j);
            if v(k)~=1
                s(sp)=k;
                sp=sp+1;
            end
        end
    end
    count=count+1;
end

function vs=count_vs(g)
n=length(g);
vs=zeros(1,n);
for i=1:n
    j=g(i);
    if i~=j
        vs(j)=vs(j)+1;
    end
end
