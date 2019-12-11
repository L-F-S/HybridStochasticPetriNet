# HybridStochasticPetriNet
Modelling budding yeast cell cycle using a deterministc model and a "hybrid" stochastic petri net model.
![Hybrid Stochastic Petri Net of all variables.](/img/SPN_completo.png)

### Repo structure:
```bash
├── deterministic
│   ├── Main_Tyson_Novak_det.m
│   ├── Tyson_Novak_det.m
│   ├── Tyson_Novak_det_variables_checkpoints.m
├── hybrid_stochastic_petri_net
│   ├── HSPN.py
├── img
│   ├── various image files
├── README.md
```

# Determinisitic model

The deterministic model is written in Matlab and should be run from inside Matlab from the Main_tyson_novak_det.m

# Hybrid stochastic model
The hybrid stochastic petri net is written in Python 3.7. 
To run:
```bash
python HSPN.py <n_steps>
```
and insert the number of steps you wish to simulate

# References:

* Tyson, John J., and Bela Novak. "Regulation of the eukaryotic cell cycle: molecular antagonism, hysteresis, and irreversible transitions." Journal of theoretical biology 210.2 (2001): 249-263.
* Mura, Ivan, and Attila Csikász-Nagy. "Stochastic Petri Net extension of a yeast cell cycle model." Journal of theoretical biology 254.4 (2008): 850-860.
* Wilkinson, Darren J. Stochastic modelling for systems biology. Chapman and Hall/CRC, 2006.

# Credits:
This project was written by me (Lorenzo Federico Signorini), Marco Giulini, and Nicole Salvatori for our Mathematical Modelling course in the MSc in Quantitative and Computational Biology at the University of Trento, Italy, in January 2018.
