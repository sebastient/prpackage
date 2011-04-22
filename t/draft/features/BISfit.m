(* ::Package:: *)

(* ::Input:: *)
(**)


BeginPackage["BISfit`"];


BISLSfit::usage="using least squares fit measured data, reture circle paras";
ColeP::usage="get the resistances(r0,r\[Infinity],\[Alpha])";
ColePxy::usage="get the resistances(r0,r\[Infinity],\[Alpha]) from x,y data";
Filterd::usage="filter bad Cole data";
BISffit::usage="fit fast using analytic solution";
BISffitxy::usage="fit fast using analytic solution and x y data";
BISLSfitxy::usage="using least squares fit measured x,y data, reture circle paras";
LSfit::usage="raw fitted solution";



Maxd::usage="upper boundary";\
Mind::usage="lower boundary";\


Begin["`Private`"];



LSfit[x_,y_]:=
Module[{F},
F[m_,n_,r_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[x]\)]
\*SuperscriptBox[\((
\*SuperscriptBox[\((x[\([k]\)] - m)\), \(2\)] + 
\*SuperscriptBox[\((y[\([k]\)] - n)\), \(2\)] - 
\*SuperscriptBox[\(r\), \(2\)])\), \(2\)]\);
NSolve[{\!\(
\*SubscriptBox[\(\[PartialD]\), \(m\)]\(F[m, n, r]\)\)==0&&\!\(
\*SubscriptBox[\(\[PartialD]\), \(n\)]\(F[m, n, r]\)\)==0&&\!\(
\*SubscriptBox[\(\[PartialD]\), \(r\)]\(F[m, n, r]\)\)==0},{m,n,r}]]


ColeP[data_,ff_]:=
Module[{m,n,r,t,a,R,X,r0,r1,fc},
{m,n,r}=BISffit[data];
(*ff=*)
(*Print[{m,n,r}]/.sol[[i]];*)
R=Map[Function[x,x[[1]]Cos[-x[[2]]\[Pi]/(180)]],data];
X=Map[Function[x,x[[1]]Sin[-x[[2]]\[Pi]/(180)]],data];
t=Sqrt[r^2-n^2];
a=2/\[Pi] ArcCos[-n/r];
r0=m+t;
r1=m-t;
n=Length[ff];
(*freq is a little harder*)
fc=-Re[1/n \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(ff[\([i]\)] Im[\(-
\*RadicalBox[
FractionBox[\(X[\([i]\)] I + R[\([i]\)] - r1\), \(r0 - X[\([i]\)] I - R[\([i]\)]\)], \(a\)]\)]\)\)];
{r0,r1,a,fc}
]
(*ff is the freqs,(m,n) is center of circle, r the radius*)
ColeP[m_,n_,r_,data_,ff_]:=
Module[{t,a,R,X,r0,r1,fc,n},
R=Map[Function[x,x[[1]]Cos[-x[[2]]\[Pi]/(180)]],data];
X=Map[Function[x,x[[1]]Sin[-x[[2]]\[Pi]/(180)]],data];
t=Sqrt[r^2-n^2];
a=2/\[Pi] ArcCos[-n/r];
r0=m+t;
r1=m-t;
n=Length[ff];
fc=-Re[1/n \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(ff[\([i]\)] Im[\(-
\*RadicalBox[
FractionBox[\(X[\([i]\)] I + R[\([i]\)] - r1\), \(r0 - X[\([i]\)] I - R[\([i]\)]\)], \(a\)]\)]\)\)];
{r0,r1,a,fc}
];

ColePxy[data_,ff_]:=
Module[{m,n,r,t,a,R,X,r0,r1,fc},
{m,n,r}=BISffitxy[data];
(*ff=*)
(*Print[{m,n,r}]/.sol[[i]];*)
R=data[[All,1]];
X=data[[All,2]];
t=Sqrt[r^2-n^2];
a=2/\[Pi] ArcCos[-n/r];
r0=m+t;
r1=m-t;
n=Length[ff];
(*freq is a little harder*)
fc=-Re[1/n \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(ff[\([i]\)] Im[\(-
\*RadicalBox[
FractionBox[\(X[\([i]\)] I + R[\([i]\)] - r1\), \(r0 - X[\([i]\)] I - R[\([i]\)]\)], \(a\)]\)]\)\)];
{r0,r1,a,fc}
];
ColePxy[m_,n_,r_,data_,ff_]:=
Module[{t,a,R,X,r0,r1,fc},
R=data[[All,1]];
X=data[[All,2]];
t=Sqrt[r^2-n^2];
a=2/\[Pi] ArcCos[-n/r];
r0=m+t;
r1=m-t;
n=Length[ff];
fc=-Re[1/n \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(ff[\([i]\)] Im[\(-
\*RadicalBox[
FractionBox[\(X[\([i]\)] I + R[\([i]\)] - r1\), \(r0 - X[\([i]\)] I - R[\([i]\)]\)], \(a\)]\)]\)\)];
{r0,r1,a,fc}
];


BISLSfit[data_]:=
Module[{X,Y,Z,i=0,sol,t},
X=Map[Function[x,x[[1]]Cos[-x[[2]]\[Pi]/(180)]],data];
Y=Map[Function[x,x[[1]]Sin[-x[[2]]\[Pi]/(180)]],data];
(*Z=Table[{X[[i]],Y[[i]]},{i,Length[X]}];*)
(*While[Z[[2,2]]<Z[[1,2]] &&  i<Length[Z]/3,Z=Drop[Z,1];i++];*)
sol=LSfit[X,Y];
i=1;
While[i<=Length[sol],
If[r<=0 /.sol[[i]],i++,
If[Im[m]>0 || Im[n]>0/.sol[[i]],i++,Break[]]]
];
{m,n,r}/.sol[[i]]
]

BISLSfitxy[data_]:=
Module[{X,Y,Z,i=0,sol,t},
X=data[[All,1]];
Y=data[[All,2]];
sol=LSfit[X,Y];
i=1;
While[i<=Length[sol],
If[r<=0 /.sol[[i]],i++,
If[Im[m]>0 || Im[n]>0/.sol[[i]],i++,Break[]]]
];
{m,n,r}/.sol[[i]]
]


BISffit[data_]:=Module[
{m,n,r,T,RR,XX},
RR=Map[Function[x,x[[1]]Cos[-x[[2]]\[Pi]/(180)]],data];
XX=Map[Function[x,x[[1]]Sin[-x[[2]]\[Pi]/(180)]],data];
T[a_,b_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[a]\)]\(a[\([i]\)]*b[\([i]\)]\)\)-1/Length[a] \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[a]\)]\(a[\([i]\)]\ \(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[b]\)]b[\([i]\)]\)\)\);
m=(T[RR,RR*RR+XX*XX]T[XX,XX]-T[XX,RR*RR+XX*XX]T[RR,XX])/(2(T[RR,RR]T[XX,XX]-(T[RR,XX])^2));
n=(T[XX,RR*RR+XX*XX]T[RR,RR]-T[RR,RR*RR+XX*XX]T[RR,XX])/(2(T[RR,RR]T[XX,XX]-(T[RR,XX])^2));
r=Sqrt[1/Length[RR] \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[RR]\)]\((
\*SuperscriptBox[\((RR[\([i]\)] - m)\), \(2\)] + 
\*SuperscriptBox[\((XX[\([i]\)] - n)\), \(2\)])\)\)];
{m,n,r}]


BISffitxy[data_]:=Module[
{m,n,r,T,RR,XX},
RR=data[[All,1]];
XX=data[[All,2]];
T[a_,b_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[a]\)]\(a[\([i]\)]*b[\([i]\)]\)\)-1/Length[a] \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[a]\)]\(a[\([i]\)]\ \(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[b]\)]b[\([i]\)]\)\)\);
m=(T[RR,RR*RR+XX*XX]T[XX,XX]-T[XX,RR*RR+XX*XX]T[RR,XX])/(2(T[RR,RR]T[XX,XX]-(T[RR,XX])^2));
n=(T[XX,RR*RR+XX*XX]T[RR,RR]-T[RR,RR*RR+XX*XX]T[RR,XX])/(2(T[RR,RR]T[XX,XX]-(T[RR,XX])^2));
r=Sqrt[1/Length[RR] \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(Length[RR]\)]\((
\*SuperscriptBox[\((RR[\([i]\)] - m)\), \(2\)] + 
\*SuperscriptBox[\((XX[\([i]\)] - n)\), \(2\)])\)\)];
{m,n,r}]



(*Options[Filterd]={Maxd->{136,76},Mind->{80,34}};*)
(* we do not really need this option*)
(*maxd and mind is r0, r1, fc *)(* max freq we used is 10^6*)
Options[Filterd]={Maxd->{200,100,1000000},Mind->{30,20,1}};
Filterd[d_,opts___]:=
Module[(*filter one person's multi-time Cole data*)
{td=d,t=0,m=0,mx1,mn1},
(*data is r0,r1,a,fc*)
For[i=1,i<=Length[td],i,
(*delete ith element if not real or <0 or a>1*)(*  fc must be real and >0*)
If[Length[Select[td[[i]],(NumberQ[#]==False||#<0)&]]>0||td[[i]][[3]]>1,td=Drop[td,{i}];Continue[],i++];
];


(*t=StandardDeviation[td];
m=Mean[td];
For[i=1,i<=Length[td],i,
(*3\[Sigma] princple and m,n,r>0 here*)
If[Length[Select[Abs[td[[i]]-m]-3t,#>0&]]!=0 ||Length[Select[td[[i]],#<0&]]!=0,td=Drop[td,{i}];Continue[],i++];
];*)
mx1=Maxd/.{opts}/.Options[Filterd];
mn1=Mind/.{opts}/.Options[Filterd];
For[i=1,i<=Length[td],i,
(*Another value boundary*)
If[td[[i]][[1]]>mx1[[1]]\[Or]td[[i]][[2]]>mx1[[2]]\[Or]td[[i]][[4]]>mx1[[3]]
   ||td[[i]][[1]]<mn1[[1]]\[Or]td[[i]][[2]]<mn1[[2]]\[Or]td[[i]][[4]]<mn1[[3]],
td=Drop[td,{i}];Continue[],i++];
];
td
];



(* ::Input:: *)
(**)


End[];


EndPackage[];
