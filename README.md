# Robust-and-Resilient-Control
This repository host the code developed during my master thesis: "Toward Robust and Resilient Cyber-Physical Systems: State awareness and Control under Physical faults and Cyber-attacks" at the IMS laboratory.

This repository is also host the code used to generate the illustrative exemple of the conference paper: *Bob Aubouin--Pairault, Arthur Perodou, Christophe Combastel, Ali Zolghadri. Resilient tube-based MPC for Cyber-Physical Systems Under DoS Attacks. 11th IFAC Symposium on Fault Detection, Supervision and Safety for Technical Processes - SAFEPROCESS 2022, Jun 2022, Pafos, Cyprus. pp.278-284*, [doi](10.1016/j.ifacol.2022.07.142). Please cite this paper if your are using the code.


# Folders Organization

* **Functions MPC:** Contains all the MPC fnction (tube based - imporved tube-based - classic - classic with zonotopes ..).
* **Resilient mpc SUN:** Contains a reproduction of the numerical example from the article of *SUN et al.* 2020 "Resilient MPC" and a correction.
* **Robust-Resilient mpc order 2:** Contains my numerical example of order 2 for the Robust and Resilient MPC proposed in my report.
* **Tube-based MPC order 2:** Contains an implementation of the Robuste tube-based MPC (both imrpoved and simple version) for an LTI system of order 2.
* **Smart building simu_1 room:** Contains simulation on the one room model (simple mpc, robuste, and robust and resilient MPC).
* **Smart building simu_8 room:** Contains simulation on the height rooms model (simple mpc and robuste MPC).
* **calcul set:** Contains various functions to calcul specific set (like terminal set, RPI set etc..)
* **Bench:** Use the BRCM toolbox to generate the one room and height rooms model. Also includes the reduction order of the height rooms model.
* **Bench:** Code to generate the figures present in our article.

# Requirement

* MPC toolbox from MATLAB
* CORA 2021 : https://tumcps.github.io/CORA/
* MPT3 toolbox : https://www.mpt3.org
* BRCM toolbox :https://brcm.ethz.ch/doku.php
