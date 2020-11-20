# MAPbI3 single crystal growth simulation by MCMC sampling
This project is done for an assignment of the computational physics course.
First, we dissolve the PbI2 to supersaturation state. With the temperature of solution increasing, PbI2 crystal will separate out and precipitate on the SiO2 substrate. Then we take out the substrate with PbI2 sediment, and put it on the Methylamine gas. Thus we can get PAPbI3 crystal.
Observed in the optical microscope, crystal of regular hexagon can be seen. So, why hexagon? As we can see the crystal face of PbI2 is parallelogram, naturally the singel crystal of PbI2 should be parallelogram. So we use simulation to prove that something has gone wrong in our experiment. Unfortunately, the result of simulation shows the hexagon our the crystal growth, fitting the experiment perfectly. I give a explanation for this phenomenon in my manuscript.
This simulation uses the Monte Carlo Markov Chain (MCMC) sampling, as MCMC sampling shows a much higher efficiency compared with the traditional sampling method.  Considering the limitted memories of PC, we can only simulate 30*30 lattice sites
