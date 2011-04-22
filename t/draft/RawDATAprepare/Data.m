(* ::Package:: *)

BeginPackage["W`Data`"];
mData::usage="measured data";
ColeData::usage="Cole fitted data";
ccData::usage="complex data";
xyData::usage="transform one group of data to x,y";
Freq23::usage="measured frequency points, 23 p";
Freq32::usage="measured frequency points, 32 p";
nxrData::usage="new R and X data";
Filter3d::usage="3\[Sigma] princple filter";



Begin["`PP`"];


ColeData[]:=Module[{ff,dd,i},
ff=FileNames["Ca_*.txt"];
dd=Reap[
For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]],"Table"]];
]];
dd[[2,1]]];

mData[]:=Module[{ff,dd,i},
ff=FileNames["Calib_*.dat"];
dd=Reap[For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]],"Table"]];
]];
dd[[2,1]]];
(*return only count groups of data*)
mData[count_]:=Module[{ff,dd,i},
ff=FileNames["Calib_*.dat"];
dd=Reap[For[i=1,i<=Min[Length[ff],count],i++,
Sow[Import[ff[[i]],"Table"]];
]];
dd[[2,1]]];


ccData[]:=Module[{ff,dd,i},
ff=FileNames["cc_*.mat"];
dd=Reap[For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]],"Mat"]];
]];
dd[[2,1]]];
ccData[count_]:=Module[{ff,dd,i},
ff=FileNames["cc_*.mat"];
dd=Reap[For[i=1,i<=Min[Length[ff],count],i++,
Sow[Import[ff[[i]],"Mat"]];
]];
dd[[2,1]]];



(*nxrData[]:=Module[{ff,dd,i},
ff=FileNames["nXR_*.mx"];
dd=Reap[For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]]]];
]];
dd[[2,1]]];*)
nxrData[]:=Module[{ff,dd,i,fm,sm},
(*first time*)
ff=FileNames["nxr_*1.txt"];
dd=Reap[For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]],"Table"]]
]
][[2,1]];
ff=Length[dd[[1,1]]];
fm=Table[Transpose[{dd[[1,All,i]],-dd[[2,All,i]]}],{i,ff}];
(*second time*)
ff=FileNames["nxr_*2.txt"];
dd=Reap[For[i=1,i<=Length[ff],i++,
Sow[Import[ff[[i]],"Table"]]
]
][[2,1]];
ff=Length[dd[[1,1]]];
sm=Table[Transpose[{dd[[1,All,i]],-dd[[2,All,i]]}],{i,ff}];
Table[{fm[[i]],sm[[i]]},{i,Length[fm]}]];


(*nxrData[count_]:=Module[{ff,dd,i},
ff=FileNames["nXR_*.mx"];
dd=Reap[For[i=1,i<=Min[Length[ff],count],i++,
Sow[Import[ff[[i]]]];
]];
dd[[2,1]]];*)
(* return all std data if count=-1 *)
nxrData[count_]:=Module[{ff,dd,i},
If[count==-1,
	ff=FileNames["Cache_stdnxy.mx"];
	If[Length[ff]!=0,Import[ff[[1]]],{}],
 If[count>0,dd=nxrData[];dd[[1;;count]],{}]
]
]


(*change data to (x,y) form*)
xyData[data_]:=Module[{ff,dd,i,X,Y,Z},
X=Map[ Function[x,x[[1]]Cos[-x[[2]]\[Pi]/(180)]],data];
Y=Map[ Function[x,x[[1]]Sin[-x[[2]]\[Pi]/(180)]],data];
Z=Table[{X[[i]],Y[[i]]},{i,Length[X]}]
];

(*Freq[]:=Module[{},
Flatten[Import["freq32.dat","Table"]]
];*)
Freq23[]:=Module[{},
Flatten[Import["freq23.txt","Table"]]
];
Freq32[]:=Module[{},
Flatten[Import["freq32.dat","Table"]]
];


(* dd must be a 2D matrix*)
Filter3d[dd_]:=Module[
{t,m,i,td},
td=dd;
t=StandardDeviation[td];
m=Mean[td];
For[i=1,i<=Length[td],i,
(*3\[Sigma] princple here*)
If[Length[Select[Abs[td[[i]]-m]-3t,#>0&]]!=0, td=Drop[td,{i}];Continue[],i++];
];
td
];


End[];
EndPackage[];
