(* ::Package:: *)

BeginPackage["W`Matcher`",{"StatisticalPlots`"}];
VF::usage="VF[true_group, predicted_group] returns Error rate, FAR and FRR";
RocPlot::usage="";
eucDist::usage="";
eucDistV::usage="Vector of subtract, eucDist mod";
dDist::usage="distance dist. dDist[data,distFunc]";
EER::usage="";
ScatterPlot::usage="scatter plot";
DistList::usage="";


Begin["`PP`"];


(*matcher*)
(*EulerDist[d1_,d2_,th_]:=Module[
{},
If[Norm[d1-d2]>th,-1,Norm[d1-d2]]
];*)



classify[data_,templ_,group_,th_]:=Module[
{dg,dl,g,tp,v,tgg},
(*data and templ should be flattened first*)
(*distance and neighbour*)
dl=Length[data];
g=DeleteDuplicates[group];

tp=Table[Extract[templ,Position[group,g[[i]]]],{i,Length[g]}];

dg=Table[0,{dl}];
For[i=1,i<=dl,i++,
v=Table[0,{Length[tp]}];
For[j=1,j<=Length[tp],j++,
tgg=Select[Table[EulerDist[da[[i]],tp[[j]][[k]],th],(*19*){k,Length[tp[[j]]]}],#>0&];
If[Length[tgg]>0,v[[j]]=Mean[tgg],v[[j]]=4096];
];
(*find the nearest*)
dg[[i]]=Ordering[v,1][[1]];
];
dg
]



(* for multi-class classification*)
(*tg:true group,pg:predicted group*)
Errorf[tg_,pg_]:=Module[
{l,er,g,gl},
If[Length[tg]!=Length[pg],Print["error input"];Return[]];
l=Length[pg];
er=N[HammingDistance[pg,tg]/l];
g=DeleteDuplicates[tg];
gl=Length[g];
(* FRR*)

For[i=1;nfr=0,i<=gl,i++,
md=Extract[dg,Position[group,g[[i]]]];
nfr=nfr+Length[Select[md,#!=g[[i]]&]];
];
(*FAR*)
For[i=1;nfa=0,i<=gl,i++,
md=Extract[group,Position[dg,g[[i]]]];
nfa=nfa+Length[Select[md,#!= g[[i]]&]];
];]


(* for binary-classification only*)
(* for distance threshold classify only *)


(*  =1: positive, =0: negative
tg=true group, pg=predicted group, 
    input should be flattened*)
VF[tg_,pg_]:=Module[
{ttg,tpg,Far,Frr,err},
(*lets flatten it first*)
{ttg,tpg}={Flatten[tg],Flatten[pg]};
If[Length[ttg]!=Length[tpg],Print["error input"];Return[]];
(*error rate:*)
err=HammingDistance[ttg,tpg]/Length[tpg]//N;
(*find the result which should be negative first, and how many are positive*)
Far=Length[Select[Extract[tpg,Position[ttg,0]],#==1&]]/Length[Select[ttg,#==0&]]//N;
(*find the truth of the negative result, ...*)
Frr=Length[Select[Extract[tpg,Position[ttg,1]],#==0&]]/Length[Select[ttg,#==1&]]//N;
{err,Far,Frr}];


EER[sdata_,templ_,tg_,fc_,step_]:=Module[
{rr,de,trr,tde=2,i},
For[i=step[[1]],i<step[[2]],i=i+step[[3]],
rr=VF[tg,fc[sdata,templ,i]][[{2,3}]];
de=Abs[rr[[1]]-rr[[2]]];
If[de<=tde,tde=de;trr=rr,Break[]]];
{i,trr}
];


(* fc is a classifier which returs err,far and frr
st=0 , return raw plot data points*) 
RocPlot[sdata_,templ_,tg_,fc_,step_,st_]:=Module[
{rr},
rr=Reap[For[i=step[[1]],i<step[[2]],i=i+step[[3]],
Sow[VF[tg,fc[sdata,templ,i]][[{2,3}]]];
]][[2,1]];
(*optimize according to:
 FAWCETT T. An introduction to roc analysis. Pattern Recognition Letters. 2006, 27 (8): 861-874.*)
(*to do*)
(*add zero*)
If[Total[rr[[1]]]!=0,rr=Prepend[rr,{0,0}]];
If[Total[rr[[-1]]]!=0,rr=Append[rr,{0,0}]];

If[st==0,rr*100,
ListPlot[rr*100,Joined->True,Frame->True,PlotRange->{{0,100},{0,100}},
FrameLabel->{"FAR (%)","FRR (%)"},AxesOrigin->{0,0},PlotStyle->Dashed]]
]


(* each row of pm is a class,each element in a row is a sample*)
eucDist[pm_]:=Module[
{dsame,ddiff,i,j,k},

ddiff=Flatten[Reap[For[i=1,i<=Length[pm]-1,i++,
For[j=i+1,j<=Length[pm],j++,
Sow[Flatten[Outer[EuclideanDistance,pm[[i]],pm[[j]],1]]]
]]
][[2,1]]];

dsame=Reap[For[i=1,i<=Length[pm],i++,
For[j=1,j<=Length[pm[[i]]]-1,j++,
For[k=j+1,k<=Length[pm[[i]]],k++,
Sow[EuclideanDistance[pm[[i]][[j]],pm[[i]][[k]]]];
(*Sow[-Mdistance[fdata[[i]][[j]],fdata[[i]][[k]]]];*)
]]]
][[2,1]];
{ddiff,dsame}
];

eucDistV[pm_]:=Module[
{dsame,ddiff,i,j,k},

ddiff=Flatten[Reap[For[i=1,i<=Length[pm]-1,i++,
For[j=i+1,j<=Length[pm],j++,
Sow[Flatten[Outer[Subtract,pm[[i]],pm[[j]],1],1]]
]]
][[2,1]],1];

dsame=Reap[For[i=1,i<=Length[pm],i++,
For[j=1,j<=Length[pm[[i]]]-1,j++,
For[k=j+1,k<=Length[pm[[i]]],k++,
Sow[Subtract[pm[[i]][[j]],pm[[i]][[k]]]];
(*Sow[-Mdistance[fdata[[i]][[j]],fdata[[i]][[k]]]];*)
]]]
][[2,1]];
{ddiff,dsame}
];

(* func distance distribution*)
dDist[pm_,func_]:=Module[
{dsame,ddiff,i,j,k},

ddiff=Flatten[Reap[For[i=1,i<=Length[pm]-1,i++,
For[j=i+1,j<=Length[pm],j++,
Sow[Flatten[Outer[func,pm[[i]],pm[[j]],1]]]
]]
][[2,1]]];

dsame=Reap[For[i=1,i<=Length[pm],i++,
For[j=1,j<=Length[pm[[i]]]-1,j++,
For[k=j+1,k<=Length[pm[[i]]],k++,
Sow[func[pm[[i]][[j]],pm[[i]][[k]]]]
(*Sow[-Mdistance[fdata[[i]][[j]],fdata[[i]][[k]]]];*)
]]]
][[2,1]];
{ddiff,dsame}
];
(*multi dist for plot,each row of data is a dist*) 
DistList[data_,range_]:=Module[{dd,distx,pdist},
(*preprocess distribution list for plot,
range is {minx,maxx,dx}*)

dd=Map[BinCounts[#,range]/Length[#]&,data];
distx=Table[i*range[[3]],{i,0,Length[dd[[1]]]-1}];
pdist=Map[Transpose[{distx,#}]&,dd]
];


(*Draw something*)
(* pairwise scatterplot*)
ScatterPlot[data_,drange_]:=Module[{color,pp},
(* data must be form of m*n*p, m classes(usually 2),n features and p points each,
there must be an equal p for all classes *)
color=ColorData[1,"ColorList"];
If[Length[data]>Length[color],Print["too many classes"];Return[]];
If[Length[Select[data,Dimensions[#]<2&]]>0,Print["data structure error"];Return[]];
pp=Table[PairwiseScatterPlot[data[[i]],
         DataTicks->True,DataRanges->drange,
         PlotStyle->{color[[i]],PointSize[Small]}],{i,Length[data]}];
Show[pp]
];



End[];
EndPackage[];
