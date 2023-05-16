

 implementation greatly inspired from: https://github.com/geodynamics/seismic_cpml/blob/master

 # Credits 
 Author: Pawan Bharadwaj 
         (bharadwaj.pawan@gmail.com)
 
 * original code in FORTRAN90: March 2013
 * modified: 11 Sept 2013
 * major update: 25 July 2014
 * code optimization with help from Jan Thorbecke: Dec 2015
 * rewritten in Julia: June 2017
 * added parrallelization over supersources in Julia: July 2017
 * efficient parrallelization using distributed arrays: Sept 2017
 * optimized memory allocation: Oct 2017



 
As forward modeling method, the 
finite-difference method is employed. 
It uses a discrete version of the two-dimensional isotropic acoustic wave equation.

```math
\pp[\tzero] - \pp[\tmo] = \dt \mB \left({\partial_x\vx}[\tmh]
 + \partial_z \vz[\tmh]  + \dt\sum_{0}^{\tmo}\sfo\right)
 ```
 ```math
\pp[\tpo] - \pp[\tzero] = \dt \mB \left(\partial_x \vx[\tph]
 + {\partial_z \vz}[\tph]  + \dt\sum_{0}^{\tzero}\sfo\right)
 ```

#attenuation related
#Rmemory::StackedArray2DVector{Array{Float64,3}} # memory variable for attenuation, see Carcione et al (1988), GJ
#Rmemoryp::StackedArray2DVector{Array{Float64,3}} # at previous time step