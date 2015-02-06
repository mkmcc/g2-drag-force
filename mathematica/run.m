(* separate out code into separate files *)
Get["preliminaries.m"];
Get["initial-cond.m"];
Get["mk-orbit.m"]
Get["chi-sq.m"];

dobs = g2posdata[[All, {3, 5}]] // Max;

$params = {};
c = 0;
Print[Timing[NMinimize[{chisqerr[{alpha     -> mAlpha,
                                  beta      -> mBeta,
                                  fkep      -> mFkep,
                                  theta     -> mTheta,
                                  phi       -> mPhi,
                                  gf        -> mGF,
                                  sigmag2   -> mSigmaG2,
                                  raShift   -> mRaShift,
                                  decShift  -> mDecShift,
                                  vraShift  -> mVRaShift,
                                  vdecShift -> mVDecShift}],
                        (*  *)
                        (* parameter constraints (or "priors") *)
                        0.0       < mAlpha     < 1.2,
                        0.0       < mFkep      < 1.0,
                        0.05      < mBeta      < 100.0,
                        0.05      < mTheta     < Pi,
                        0.05      < mPhi       < 2*Pi,
                        0.0       < mGF        < 100.0,
                        1 * 10^5  < mSigmaG2   < 1 * 10^7,
                        -3 * dobs < mRaShift   < 3 * dobs,
                        -3 * dobs < mDecShift  < 3 * dobs,
                        -3 * dobs < mVRaShift  < 3 * dobs,
                        -3 * dobs < mVDecShift < 3 * dobs},
                       (*  *)
                       (* initial search space *)
                       {{mAlpha,        0.55,      0.65}, 
                        {mSigmaG2,      0.8*10^6,  1.2*10^6},
                        {mGF,          34.0,      36.0},
                        {mBeta,         0.7,       0.9}, 
                        {mTheta,        2.3,       2.4}, 
                        {mPhi,          0.1,       0.2}, 
                        {mFkep,         0.55,      0.65},
                        (*  *)
                        {mRaShift,   -0.1*dobs, 0.1*dobs},
                        {mDecShift,  -0.1*dobs, 0.1*dobs},
                        {mVRaShift,  -0.1*dobs, 0.1*dobs},
                        {mVDecShift, -0.1*dobs, 0.1*dobs}},

                       Method -> {"SimulatedAnnealing", "SearchPoints" -> 100},

                       EvaluationMonitor :> Block[{}, c++; Print["step " <> ToString[c]]]]]];


(* -- check below!  may need to manually delete the file first *)
outfile = "../data/obserr-big.symb";
Save[outfile, $icfit];
Save[outfile, $params];

(* Local Variables: *)
(* mode: mathematica *)
(* coding: utf-8-unix *)
(* End: *)
