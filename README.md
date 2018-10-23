# EFT constraints


`SL_OS_SS_combined_blind.json` contains the upper limit on tttt cross section.
Edit this file if you want to get different constrains.
## 2d plots
```bash
python EFT.py -c conf2d_cff.py
```

## independent and marginal limits tables
result is stored in build folder
```bash
python EFT.py -c conf_cff.py 
```

## API for EFT cross sections
1. Select 13 TeV of 14 TeV predictions in **mg_calcuations.py**
```python
sig_SM=sig_SM_13TeV
sig_SM=sig_SM_14TeV
# and
MG_SM=MG_SM_13TeV
MG_SM=MG_SM_14TeV
```

2. To calculate tttt cross sections as function of EFT parameters

```python
from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM, sig_SM

eft = EftPredictions(wilson_coefficients, MG_SM, sig_SM)

tttt_xs = eft.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])
# or
tttt_xs = eft.vgen_eft_xs(C_OR,C_OL1,C_OL8,C_B1,CB8) # optimzed for numpy
```
alternatively one can use 13 and 14 TeV constants directly
```python
from eft_coefficients import EftPredictions
from mg_calculations import wilson_coefficients, MG_SM_13TeV, sig_SM_13TeV, MG_SM_14TeV, sig_SM_14TeV

eft13 = EftPredictions(wilson_coefficients, MG_SM_13TeV, sig_SM_13TeV)
eft14 = EftPredictions(wilson_coefficients, MG_SM_14TeV, sig_SM_14TeV)
tttt_xs_14_to_13_ratio = eft13.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])/eft14.gen_eft_xs([C_OR,C_OL1,C_OL8,C_B1,CB8])

# test
eft13.gen_eft_xs([0.,0.,0.,0.,0.]) # should give 9.201
eft14.gen_eft_xs([0.,0.,0.,0.,0.]) # should give 11.31723
```

-------------------------------------
# Example individual plots for 1d and 2d EFT constraints

```bash
python matrix_coef.py 
# python matrix_coef_old.py  # outdated version is kept for historical purposes only. Use the script above!
```

`python matrix_coef.py` is an improved version of `python matrix_coef_old.py`. The logic of limits calculation and plotting is separated into several modules:
* `mg_calculations.py` --- Contains predictions such as tttt cross section values in a given point of the EFT parameter space. The predictions are obtained from MG, substituting a vector of Wilson coefficients in the steering cards.
* `eft_coefficients.py` --- Contains logic for the determination of analytic parametrisation of the tttt cross section as function of Wilson coefficients

**IMPORTANT**
> Results were compared with Mathias Mancini. Full agreement in the predicted values was obtained!

