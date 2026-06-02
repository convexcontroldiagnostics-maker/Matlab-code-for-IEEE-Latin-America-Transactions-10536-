# A Quasi-LPV Dynamic Output Feedback Stabilizer for Nonlinear Descriptor Systems via Convex

**Manuscript ID:**  IEEE Latin America Transactions Submission ID: 10536 

**Authors:**  
- Arturo Alvarado, Department of Mechatronics of the Universidad Politécnica de Pachuca, C.P. 43830, Zempoala, Hidalgo, Mexico
- Tonatiuh Hernández Cortés, Universidad Autónoma del Estado de Hidalgo, Instituto de Ciencias Básicas e Ingeniería, Pachuca, Hidalgo, Mexico.  
- Jaime González Sierra, Unidad Profesional Interdisciplinaria de Ingeniería Campus Hidalgo, Instituto Politécnico Nacional, Carretera Pachuca—Actopan Kilómetro 1+500, Distrito de Educación, Salud, Ciencia, Tecnología e Innovación, C.P. 42162, San Agustín Tlaxiaca, Hidalgo, Mexico 
- Víctor Estrada Manzo, Department of Mechatronics of the Universidad Politécnica de Pachuca, C.P. 43830, Zempoala, Hidalgo, Mexico  

---

## 📁 Included Scripts

This repository contains all scripts required to reproduce the simulation and numerical results presented in the article.

| Script | Related Figure(s) | Description |
|--------|-------------------|-------------|
| `Exa_1b.m`| None | Compute LMI conditions in Theorems 1 and 2 for Example 1. |
| `Exa_1a.m`, `Exa_1.slx` | Figs. 2 and 3 | Implement the controler (23) for Example 1. |
| `Exa1a_Fang.m`,`Exa1_Fang.slx` | Fig. 4 | Implement the controller in [26]. |
| `Exa_2grid.m`, `Controller_Observer_Th1_Th2.m` | Fig. 5 | Compute the feasible set solution for LMI conditions in Theorems 1 and 2 while compares with previous results [27] . |
|`data_Ricardo.mat` | Fig. 5 | Contains the data from the approach in [27] . |
| `Exa_3RIP_control.m`, `Exa_3RIP_obs.m` | None | Compute the LMI conditions in Theorems 1 and 2 for the RIP. |
| `Experimental_Data.mat` | Figs 7,8,9,10 | Contains the data from the real-time implementation of the RIP system. |

---


## 💻 Requirements

- MATLAB R2015a or later.

---

## ✉️ Contact

For questions or replication of results:  
victor_estrada@upp.edu.mx
