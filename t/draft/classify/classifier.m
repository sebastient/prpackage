(* ::Package:: *)

BeginPackage["W`Classifier`","KNN`"];(* load knn*)



classifyDist1::usage="euclidean distance classfier";
classifyDist2::usage="Cosine distance classfier";


Begin["`Private`"];


(*the group of tmp and classifier result should be (number of sample) row 
and (number of class) column*)  


(* template tmp has one sample for each class, indata has unknown numbers of samples*)
(*indata should first flattened*)
classifyDist1[indata_,tmp_,th_]:=Module[
{group,i,j},
group=ConstantArray[ConstantArray[0,Length[tmp]],Length[indata]];
For[i=1,i<=Length[indata],i++,
For[j=1,j<=Length[tmp],j++,
If[EuclideanDistance[indata[[i]],tmp[[j]]]<th,
group[[i,j]]=1];
]];
group
]
(* using cos dist*)
classifyDist2[indata_,tmp_,th_]:=Module[{},
classifyDist[indata,tmp,th,CosineDistance]
]
(* give your own distance*)
classifyDist[indata_,tmp_,th_,dfunc_]:=Module[
{group,i,j},
group=ConstantArray[ConstantArray[0,Length[tmp]],Length[indata]];
For[i=1,i<=Length[indata],i++,
For[j=1,j<=Length[tmp],j++,
If[dfunc[indata[[i]],tmp[[j]]]<th,
group[[i,j]]=1];
]];
group
]


End[];

EndPackage[];
