(* : Title: BootStrapPackage *)

(* : Author: Enis Siniksaran *) 

(* : Mathematica Version: 4.0 *)

(* : Date: June,2001 *)

(* : Summary : This package is designed to produce some summary statistics 
and confidence intervals based on bootstrap percentiles and BCa percentiles
for the mean, simple linear regression and correlation coefficient. *)

BeginPackage["BootStrapPackage`"]

BootMean::usage = "BootMean[list,B] returns a table which shows some summary 
statistics including bootstrap estimate of mean, standard error and bias, confidence intervals 
based on bootstrap percentiles and BCa percentiles. B which is the number of 
bootstrap samples must be an integer.It also gives a histogram of B bootstrap replications of 
the mean. Dashing line indicates the mean of the original sample";

BootRegPair::usage = "BootRegPair[data,B] finds a bootstrap fit for simple 
linear regression model.The method is based on the bootstrapping pairs.The data must have the 
form {{x1, y1 },{x2, y2}, ... ,{xk , yk}}. B which is the number of bootstrap 
samples must be an integer.BootRegPair[list,B] returns a table which shows summary 
statistics for model parameters beta0 and beta1.It also gives histograms of B bootstrap 
replications of the estimates of the model parameters. Dashing lines indicate the 
estimates of the beta0 and beta1 for the original data";

BootRegRes::usage = "BootRegRes[data,B] finds a bootstrap fit for simple 
linear regression model.The method is based on the bootstrapping residuals.The data must have 
the form {{x1, y1 },{x2, y2}, ... ,{xk , yk}}. B which is the number of bootstrap 
samples must be an integer.BootRegRes[list,B] returns a table which shows summary 
statistics for model parameters beta0 and beta1.It also gives histograms of B bootstrap 
replications of the estimates of the model parameters. Dashing lines indicate 
estimates of the beta0 and beta1 for the original data";


BootCorr::usage ="BootCorr[data,B] returns a table which shows some summary 
statistics including bootstrap estimate of correlation coefficient, standard error and bias, 
confidence intervals based on bootstrap percentiles and BCa percentiles.The data must have the 
form {{x1, y1 },{x2, y2}, ... ,{xk , yk}}.B which is the number of bootstrap 
samples must be an integer.It also gives a histogram of B bootstrap replications of 
the correlation coefficient.Dashing line indicates the estimate of the correlation 
coefficient for the original data";

Begin["`Private`"]

Needs["DiscreteMath`Combinatorica`"]
Needs["Statistics`MultiDescriptiveStatistics`"]
Needs["Statistics`NormalDistribution`"]
Needs["Statistics`LinearRegression`"]
Needs["Graphics`Graphics`"]

            (* Bootstrapping The Mean *)

BootMean[list_List, B_Integer?Positive] := 
  Module[{eboos, eboot, i, means, meanboot, k, bias, seobserved, seboot, 
      meanobserved, zetzero, alphahat, p, ones, ejack, ejacks, mjack, 
      meanjacks, percentiles, BCas, hist}, 
      k = Length[list];
    Do[eboot[i] = Table[RandomKSubset[list, 1], {k}], {i, 1, B}];
    eboos = Table[eboot[i], {i, 1, B}];
    means = Table[Mean[Flatten[eboos[[i]], 1]], {i, 1, B}];
    meanboot = N[Mean[means]];
    seboot = N[StandardDeviation[means]];
    meanobserved = N[Mean[list]];
    biasmean = meanobserved - meanboot;
    seobserved = (StandardDeviation[list]/
              Sqrt[Length[list]])*(Sqrt[Length[list] - 1]/
              Sqrt[Length[list]]) // N;
    p = {.005, .025, .05, .95, .975, .995};
    percentiles = Table[N[Quantile[means, p[[i]]]], {i, 1, 6}];
    Do[ejack[i] = Drop[list, {i}], {i, 1, k}];
    ejacks = Table[ejack[i], {i, 1, B}];
    meanjacks = N[Table[Mean[ejacks[[i]]], {i, 1, k}]];
    mjack = N[Mean[meanjacks]];
    ones = Table[1, {i, 1, k}];
    alphahat = ((mjack - 
                  meanjacks)^3).ones/(6*(((mjack - meanjacks)^2).ones)^(3/2));
    zetzero = 
      N[Quantile[NormalDistribution[0, 1], 
          N[Length[Select[means, # < Mean[list] &]]/B]]];
    BCas = 
      Table[N[Quantile[means, 
            N[CDF[NormalDistribution[0, 1], 
                zetzero + ((zetzero + 
                          N[Quantile[NormalDistribution[0, 1], p[[i]]]])/
                          (1 - alphahat*(zetzero + 
                                N[Quantile[NormalDistribution[0, 1], 
                                    p[[i]]]])))]]]], {i, 1, 6}];
    Print["*", "*", "Bootstrap Results:", B "Replications", "*", "*"];   
    Print["Summary Statistics"];
    Print[
      TableForm[{{"meanobserved", "meanboot", "bias", "seobserved", 
            "seboot"}, {meanobserved, meanboot, biasmean, seobserved, 
            seboot}}]];
    Print["Bootstrap Percentiles"];
    Print[
      TableForm[{{"0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"}, 
          percentiles}]];
    Print["BCa Percentiles"];
    Print[
      TableForm[{{"0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"}, 
          BCas}]];          
     hist = Histogram[means, BarStyle -> {RGBColor[1, 1, 1]}, 
        Epilog -> {Dashing[{.01}], Line[{{ meanobserved, 0}, { meanobserved, B}}]}, 
        AxesLabel -> {"Bootstrap Means", "Frequencies"}, 
        ImageSize -> {500, 200}];
    ];
    

             (* Bootstrapping Regression:Pairs *)

BootRegPair[data_, B_Integer?Positive] := 
  Module[{regress, bestfitparameters, beta0, beta1, popboot, eboot, i, eboos,

      regressboot, parameterestimates, pairbeta0, pairbeta1, regressjack, 
      parameterestimatesboot, beta0boots, beta1boots, beta0boot, beta1boot, 
      sebootbeta0, sebootbeta1, ejack, ejacks, parameterestimatesjack, 
      beta0jack, beta1jack, beta0jacks, beta1jacks, k, p, ones, alphahat, 
      zetzero, BCas, histbeta0, histbeta1}, 
    regress = 
        Regress[data, {1, x}, x, RegressionReport -> {BestFitParameters}];;
    bestfitparameters = BestFitParameters /. regress;
    beta0 = bestfitparameters[[1]];
    beta1 = bestfitparameters[[2]];
    popboot = Flatten[Join[Table[data, {i, 1, Length[data]}]], 1];
    Do[eboot[i] = RandomKSubset[popboot, Length[data]], {i, 1, B}];
    eboos = Table[eboot[i], {i, 1, B}];
    regressboot = 
      Table[Regress[eboos[[i]], {1, x}, x, 
          RegressionReport -> {BestFitParameters}], {i, 1, B}];
    parameterestimates = (bestfitparameters = 
          BestFitParameters /. regressboot);
    beta0boots = Transpose[parameterestimates][[1]];
    beta1boots = Transpose[parameterestimates][[2]];
    beta0boot = Mean[beta0boots];
    beta1boot = Mean[beta1boots];
    sebootbeta0 = StandardDeviation[beta0boots];
    sebootbeta1 = StandardDeviation[beta1boots];
    biasbeta0 = beta0boot - beta0;
    biasbeta1 = beta1boot - beta1;
    k = Length[data];
    Do[ejack[i] = Drop[data, {i}], {i, 1, k}];
    ejacks = Table[ejack[i], {i, 1, k}];
    regressjack = 
      Table[Regress[ejacks[[i]], {1, x}, x, 
          RegressionReport -> {BestFitParameters}], {i, 1, k}];
    parameterestimatesjack = (bestfitparameters = 
          BestFitParameters /. regressjack);
    beta0jacks = Transpose[parameterestimatesjack][[1]];
    beta1jacks = Transpose[parameterestimatesjack][[2]];
    beta0jack = Mean[beta0jacks];
    beta1jack = Mean[beta1jacks];
    p = {.005, .025, .05, .95, .975, .995};
    percentilesbeta0 = Table[N[Quantile[beta0boots, p[[i]]]], {i, 1, 6}];
    percentilesbeta1 = Table[N[Quantile[beta1boots, p[[i]]]], {i, 1, 6}];
    ones = Table[1, {i, 1, k}];
    alphahatbeta0 = ((beta0jack - 
                  beta0jacks)^3).ones/(6*(((beta0jack - 
                          beta0jacks)^2).ones)^(3/2));
    zetzerobeta0 = 
      N[Quantile[NormalDistribution[0, 1], 
          N[Length[Select[beta0boots, # < beta0 &]]/B]]];
    BCasbeta0 = 
      Table[N[Quantile[beta0boots, 
            N[CDF[NormalDistribution[0, 1], 
                zetzerobeta0 + ((zetzerobeta0 + 
                          N[Quantile[NormalDistribution[0, 1], p[[i]]]])/
                          (1 - alphahatbeta0*(zetzerobeta0 + 
                                N[Quantile[NormalDistribution[0, 1], 
                                    p[[i]]]])))]]]], {i, 1, 6}];
    alphahatbeta1 = ((beta1jack - 
                  beta1jacks)^3).ones/(6*(((beta1jack - 
                          beta1jacks)^2).ones)^(3/2));
    zetzerobeta1 = 
      N[Quantile[NormalDistribution[0, 1], 
          N[Length[Select[beta1boots, # < beta1 &]]/B]]];
    BCasbeta1 = 
      Table[N[Quantile[beta1boots, 
            N[CDF[NormalDistribution[0, 1], 
                zetzerobeta1 + ((zetzerobeta1 + 
                          N[Quantile[NormalDistribution[0, 1], p[[i]]]])/
                          (1 - alphahatbeta1*(zetzerobeta1 + 
                                N[Quantile[NormalDistribution[0, 1], 
                                    p[[i]]]])))]]]], {i, 1, 6}];
    Print["*", "*", "Bootstrap Results for Paired Regression:", 
      B "Replications", "*", "*"];   Print["Summary Statistics"];
    Print[
      TableForm[{{"", "observed", "bootstrap", "bias", "seboot"}, {"beta0", 
            beta0, beta0boot, biasbeta0, sebootbeta0}, {"beta1", beta1, 
            beta1boot, biasbeta1, sebootbeta1}}]];
    Print["Bootstrap Percentiles"];
    Print[
      TableForm[{{"","0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"}, 
         Prepend[percentilesbeta0,"beta0"],
         Prepend[percentilesbeta1,"beta1"]}]];
    Print["BCa Percentiles"];
    Print[
      TableForm[{{"","0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"}, 
      Prepend[BCasbeta0,"beta0"],Prepend[BCasbeta1,"beta1"]}]];          
    histbeta0 = 
      Histogram[beta0boots, BarStyle -> {RGBColor[1, 1, 1]}, 
        Epilog -> {Dashing[{.01}], Line[{{beta0, 0}, {beta0, B}}]}, 
        AxesLabel -> {"Bootstrap Beta0", "Frequencies"}, 
        ImageSize -> {500, 200}];
    histbeta1 = 
      Histogram[beta1boots, BarStyle -> {RGBColor[1, 1, 1]}, 
        Epilog -> {Dashing[{.01}], Line[{{beta1, 0}, {beta1, B}}]}, 
        AxesLabel -> {"Bootstrap Beta1", "Frequencies"}, 
        ImageSize -> {500, 200}];
        ];
        

                    (* Bootstrapping Regression:Residuals *)

  
    BootRegRes[data_, B_Integer?Positive] :=Module[{listx, listy, ones, 
datamatrix, regress, residuals , bestfitparameters, varyhat, beta0, beta1, k, i, eboot, eboos, yy, 
booty,g, beta0boots, beta1boots, beta0boot, beta1boot, sebootbeta0, 
      sebootbeta1, biasbeta0, biasbeta1, regressjack, ejack, ejacks, 
      parameterestimatesjack, beta0jack, beta1jack, beta0jacks, beta1jacks, 
p,alphahat, zetzero, BCas, histbeta0, histbeta1},
 
        
    listx = Transpose[data][[1]];
    listy = Transpose[data][[2]];
    ones = Table[1, {i, 1, Length[listx]}];
    datamatrix = Transpose[{ones, listx}];



regress = Regress[data, {1, x}, x,
RegressionReport ->
 {FitResiduals, SinglePredictionCITable,
BestFitParameters
 }];

residuals = FitResiduals /. regress;
bestfitparameters = BestFitParameters /. regress;
varyhat = Transpose[(SinglePredictionCITable /. regress)[[1]]][[2]];
    
beta0 = bestfitparameters[[1]];
beta1 = bestfitparameters[[2]];

k = Length[data];

    
Do[eboot[i] = Table[RandomKSubset[residuals, 1], {k}], {i, 1, B}];
eboos = Table[eboot[i], {i, 1, B}];
    
yy = Table[varyhat, {B}];
booty = yy + eboos;
g = Inverse[
          Transpose[datamatrix].datamatrix].(Transpose[datamatrix].
          Transpose[booty]);

     beta0boots = Flatten[g[[1]]];
     beta1boots = Flatten[g[[2]]];
     beta0boot = Mean[beta0boots];
     beta1boot = Mean[beta1boots];
   
    sebootbeta0 = StandardDeviation[beta0boots]; 
    sebootbeta1 = StandardDeviation[beta1boots];
    
    biasbeta0 = beta0boot - beta0;
    biasbeta1 = beta1boot - beta1;
        
    Do[ejack[i] = Drop[data, {i}], {i, 1, k}];
    ejacks = Table[ejack[i], {i, 1, k}];
    
    regressjack = 
      Table[Regress[ejacks[[i]], {1, x}, x, 
          RegressionReport -> {BestFitParameters}], {i, 1, k}];
    parameterestimatesjack = (bestfitparameters = 
          BestFitParameters /. regressjack);
    beta0jacks = Transpose[parameterestimatesjack][[1]];
    beta1jacks = Transpose[parameterestimatesjack][[2]];
    beta0jack = Mean[beta0jacks];
    beta1jack = Mean[beta1jacks];
    p = {.005, .025, .05, .95, .975, .995};
    percentilesbeta0 = Table[N[Quantile[beta0boots, p[[i]]]], {i, 1, 6}];
    percentilesbeta1 = Table[N[Quantile[beta1boots, p[[i]]]], {i, 1, 6}];
    ones = Table[1, {i, 1, k}];
    alphahatbeta0 = ((beta0jack - 
                  beta0jacks)^3).ones/(6*(((beta0jack - 
                          beta0jacks)^2).ones)^(3/2));
    zetzerobeta0 = 
      N[Quantile[NormalDistribution[0, 1], 
          N[Length[Select[beta0boots, # < beta0 &]]/B]]];
    BCasbeta0 = 
      Table[N[Quantile[beta0boots, 
            N[CDF[NormalDistribution[0, 1], 
                zetzerobeta0 + ((zetzerobeta0 + 
                          N[Quantile[NormalDistribution[0, 1], p[[i]]]])/
                          (1 - alphahatbeta0*(zetzerobeta0 + 
                                N[Quantile[NormalDistribution[0, 1], 
                                    p[[i]]]])))]]]], {i, 1, 6}];
    alphahatbeta1 = ((beta1jack - 
                  beta1jacks)^3).ones/(6*(((beta1jack - 
                          beta1jacks)^2).ones)^(3/2));
    zetzerobeta1 = 
      N[Quantile[NormalDistribution[0, 1], 
          N[Length[Select[beta1boots, # < beta1 &]]/B]]];
    BCasbeta1 = 
      Table[N[Quantile[beta1boots, 
            N[CDF[NormalDistribution[0, 1], 
                zetzerobeta1 + ((zetzerobeta1 + 
                          N[Quantile[NormalDistribution[0, 1], p[[i]]]])/
                          (1 - alphahatbeta1*(zetzerobeta1 + 
                                N[Quantile[NormalDistribution[0, 1], 
                                    p[[i]]]])))]]]], {i, 1, 6}];
    Print["*", "*", "Bootstrap Results for Residual Regression:", 
      B "Replications", "*", "*"];   Print["Summary Statistics"];
    Print[
      TableForm[{{"", "observed", "bootstrap", "bias", "seboot"}, {"beta0", 
            beta0, beta0boot, biasbeta0, sebootbeta0}, {"beta1", beta1, 
            beta1boot, biasbeta1, sebootbeta1}}]];
    Print["Bootstrap Percentiles"];
    Print[
      TableForm[{{"","0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"}, 
          Prepend[percentilesbeta0,"beta0"],Prepend[percentilesbeta1,"beta1"]}]];
    Print["BCa Percentiles"];
    Print[
      TableForm[{{"","0.5%", "2.5%", "5%", "95%", "97.5%", "99.5%"},  
      Prepend[BCasbeta0,"beta0"],Prepend[BCasbeta1,"beta1"]}]];         
    histbeta0 = 
      Histogram[beta0boots, BarStyle -> {RGBColor[1, 1, 1]}, 
        Epilog -> {Dashing[{.01}], Line[{{beta0, 0}, {beta0, B}}]}, 
        AxesLabel -> {"Bootstrap Beta0", "Frequencies"}, 
        ImageSize -> {500, 200}];
    histbeta1 = 
      Histogram[beta1boots, BarStyle -> {RGBColor[1, 1, 1]}, 
        Epilog -> {Dashing[{.01}], Line[{{beta1, 0}, {beta1, B}}]}, 
        AxesLabel -> {"Bootstrap Beta1", "Frequencies"}, 
        ImageSize -> {500, 200}];
        ];
    
    
                      (* Bootstrapping The Correlation Coefficient *)
    
     BootCorr[data_,B_Integer?Positive]:=

          Module[{corr,k,i,eboos,bootmatrix,corrsboot,corrboot,sebootcorr,biascorr,
      ejack,ejacks,  cojacks,corrjacks,corrjack,p,percentilescorr,ones,
      alphahat,zetzero,BCas,hist},
    
    
    corr=N[Correlation[Transpose[data][[1]],Transpose[data][[2]]]];
        k=Length[data];
    Do[eboot[i]=Table[RandomKSubset[data,1],{k}],{i,1,B}];
    eboos=Table[Flatten[Table[eboot[i],{i,1,B}][[i]],1],{i,1,B}];
    
    bootmatrix=Table[CorrelationMatrix[eboos[[i]]],{i,1,B}];
    
         corrboots=N[Transpose[Transpose[bootmatrix][[1]]][[2]]];
         corrboot=N[Mean[corrboots]];
    sebootcorr=N[StandardDeviation[corrboots]];
    
    biascorr=corrboot-corr;
    
    Do[ejack[i]=Drop[data,{i}],{i,1,k}];
    ejacks=Table[ejack[i],{i,1,k}];
    cojacks=Table[CorrelationMatrix[ejacks[[i]]],{i,1,k}];
    corrjacks=Transpose[Transpose[cojacks][[1]]][[2]];
    corrjack=N[Mean[corrjacks]];
    
    
    p={.005,.025,.05,.95,.975,.995};
    
    percentiles=Table[N[Quantile[corrboots,p[[i]]]],{i,1,6}];
    
    ones=Table[1,{i,1,k}];
    alphahat=((corrjack-
                  corrjacks)^3).ones/(6*(((corrjack-corrjacks)^2).
                  ones)^(3/2));
    zetzero=
      N[Quantile[NormalDistribution[0,1],
          N[Length[Select[corrboots,#<corrboot&]]/B]]];
    BCas=Table[
        N[Quantile[corrboots,
            N[CDF[NormalDistribution[0,1],
                zetzero+((zetzero+
                          N[Quantile[NormalDistribution[0,1],p[[i]]]])/
                          (1- alphahat*(zetzero+
                                N[Quantile[NormalDistribution[0,1],
                                    p[[i]]]])))]]]],{i,1,6}];
    
    Print["*","*","Bootstrap Results for Correlation:",B "Replications","*",
      "*"];
    Print["Summary Statistics"];
    Print[
      TableForm[{{"observed","boot","bias","seboot"},{corr,corrjack,biascorr,
            sebootcorr}}]];
    Print["Bootstrap Percentiles"];
    Print[
      TableForm[{{"0.5%","2.5%","5%","95%","97.5%","99.5%"},percentiles}]];
    Print["BCa Percentiles"];
    Print[TableForm[{{"0.5%","2.5%","5%","95%","97.5%","99.5%"},BCas}]];
    
    
    hist=Histogram[corrboots,BarStyle->{RGBColor[1,1,1]},
        Epilog->{Dashing[{.01}],Line[{{corr,0},{corr,B}}]},
        AxesLabel->{"Bootstrap Correlations","Frequencies"},
        ImageSize->{500,200}];
    ]
         
        
        
End[]
EndPackage[]
