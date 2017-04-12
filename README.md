# LangmuirProbe_Matlab
Calculates plasma parameters from swept Langmuir probe I-V measurements.

Performs simple linear fits to the data.

## Inputs
The input data for the code are two files: single column ASCII data for probe voltage, and single column ASCII data for the measured probe current.

## Example for Matlab Command Window: 
%       (First fitting with default parameters)
%           >> [pAns,pFit,rFit] = analyze_lp();
%       (Repeat fit and adjust parameters)
%           >> [pAns,pFit,rFit] = analyze_lp(pAns,pFit,rFit);
%       (Repeat until satisfied with fit)
%
%       (Another data set fitting but using default fit parameters)
%           >> [pAns,pFit,rFit] = analyze_lp(pAns);
%       (Only the initial prompt is filled out with the last input)

