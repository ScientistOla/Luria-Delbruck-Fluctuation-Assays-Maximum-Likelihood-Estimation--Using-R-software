data=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 3 3 3 3 3 4 4 6 7 7 7 7 8 8 9 10 11 11 12 13 13 13 15 16 21 27 31 35 81 94 102 145];
m_d=findMLmTwoParam(data);
m=m_d(1,1);
d=m_d(1,2);
disp("m = ");
disp(m);
disp("d = ");
disp(d);
maxValue=max(data);
SSDScoreLDResult=SSDScoreLD(data);
disp(SSDScoreLDResult);
SSDScoreTwoParamResult=SSDScoreTwoParam(data);
dataLength=length(data);
SSDScoreLDResultForDraw=zeros(1,dataLength);
for c = 1:dataLength
SSDScoreLDResultForDraw(1,c)=SSDScoreLDResult(1,(data(1,c))+1);
end
SSDScoreTwoParamResultForDraw=zeros(1,dataLength);
for c = 1:dataLength
SSDScoreTwoParamResultForDraw(1,c)=SSDScoreTwoParamResult(1,(data(1,c))+1);
end
ln =plot(data,SSDScoreTwoParamResultForDraw,"ro");
ln.LineWidth=2.5;
ln.Color=[0.4940 0.1840 0.5560];
hold on
ln3=owerCdfplot(data);
ln3.Marker=".";
ln2= plot(data,SSDScoreLDResultForDraw);
ln2.LineWidth=2;
ln2.Color = [0.4660 0.6740 0.1880];
ln2.MarkerEdgeColor = 'b';
xlabel('Mutants per culture ')
ylabel('Cumulative distribution')
title('Clone 1 ( + Tet ) ')
legend('Two-parameter Model ','Data','One-parameterModel',"location","southeast");
hold off