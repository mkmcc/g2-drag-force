(* ============================================================================== *)
(* given a numerical solution to the orbit, calculate a chi^2 against the data *)
chisqG2[sol_, tbeg_, tend_] :=
    Module[{chisqpos, chisqvel,
            myposdata = Select[g2posdata, First[#] < tend &],
            myveldata = Select[g2veldata, First[#] < tend &]},

           chisqpos = Map[Apply[With[{t = #1, ra = #2, raerr = #3, dec = #4, decerr = #5},
                                     (  (r[t].{1,0,0} - ra)^2  / raerr^2
                                      + (r[t].{0,1,0} - dec)^2 / decerr^2) /. sol] &, #] &,
                          myposdata];

           chisqvel = Map[Apply[With[{t = #1, vlos = #2, verr = #3},
                                     ((v[t].{0,0,1} /. sol) - vlos)^2 / verr^2] &, #] &,
                          myveldata];


           (Total[chisqpos] + Total[chisqvel])];


chisqG1[sol_, tbeg_, tend_, dtg1_?NumericQ] :=
    Module[{chisqpos, chisqvel,
            myposdata = Select[g1posdata, First[#] < tend - dtg1 &],
            myveldata = Select[g1veldata, First[#] < tend - dtg1 &]},

           chisqpos = Map[Apply[With[{t = #1 + dtg1, ra = #2, raerr = #3, dec = #4, decerr = #5},
                                     (  (r[t].{1,0,0} /. sol) - ra)^2  / raerr^2
                                     + ((r[t].{0,1,0} /. sol) - dec)^2 / decerr^2] &, #] &,
                          myposdata];

           chisqvel = Map[Apply[With[{t = #1 + dtg1, vlos = #2, verr = #3},
                                     ((v[t].{0,0,1} /. sol) - vlos)^2 / verr^2] &, #] &,
                          myveldata];

           (Total[chisqpos] + Total[chisqvel])];


chisqG1[sol_, tbeg_?NumericQ, tend_?NumericQ] :=
    NMinimize[{chisqG1[sol, tbeg, tend, dtg1], dtg1 < tend - 2010.5},
              {dtg1, 0, Max[1.0, tend - 2010.5]},
              Method -> "SimulatedAnnealing"];


chisqerr[modelParams_?(VectorQ[#[[All,2]], NumericQ]&)] :=
    Module[{sol, tbeg, tend,
            chig1, chig2, chimin},

           {sol, tbeg, tend} = mkorbit[modelParams, 2040];

           chig2 = chisqG2[sol, tbeg, tend];
           {chig1, tg1rule} = chisqG1[sol, tbeg, tend];

           chimin = chig2 + chig1;

           (* use a global variable to store orbit and fit parameters *)
           (* This is pretty hacky... it requires $params to be defined before calling this function *)
           PrependTo[$params,
                     N[Flatten[{chisq -> chimin,
                                modelParams,
                                tg1rule}]]];

           chimin];


Print["loaded chi-sq.m"];

(* Local Variables: *)
(* mode: mathematica *)
(* coding: utf-8-unix *)
(* End: *)
