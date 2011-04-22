(* ::Package:: *)

(** User Mathematica initialization file **)


oldPath=$Path;
MainDir="E:\\work\\thesis\\\:8bba\:6587\\nb";
SetDirectory[MainDir];

dirs=FileNames[WordCharacter..];
DataDir="data";

myPath=Map[FileNameJoin[{MainDir,#}]&,dirs];


(* add all sub dirs in work dir and set wd to data*)
$Path=Join[oldPath, myPath];
SetDirectory[FileNameJoin[{MainDir,DataDir}]];


(*add work context*)
$ContextPath=Join[{"W`"},$ContextPath];

