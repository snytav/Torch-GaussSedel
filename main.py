import torch
import numpy as np

# solves Poisson equation using Gauss-Seidel*/
def GaussSeidel(phi,rho,dh,max_solver_it):
    EPS_0 = 8.85418782e-12  	# C/(V*m), vacuum permittivity
    idx2 = 1.0/(dh[0]*dh[0])
    idy2 = 1.0/(dh[1]*dh[1])
    idz2 = 1.0/(dh[2]*dh[2])

    L2        = 0.0			#norm
    converged = False
    ni,nj,nk = rho.shape

    for it in range(max_solver_it):
        for i in range(1,ni):
             for j in range(1,nj):
                 for k in range(1,nk):
					#standard internal open node
                    phi_new = (rho[i][j][k]/EPS_0 +
							idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
							idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
							idz2*(phi[i][j][k-1]+phi[i][j][k+1]))/(2*idx2+2*idy2+2*idz2)

					# SOR
                    phi[i][j][k] = phi[i][j][k] + 1.4*(phi_new-phi[i][j][k])


        # check for convergence
        if it%25==0:
           sum = 0
           for i in range(1, ni):
               for j in range(1,nj):
                   for k in range(1,nk):
                       R = (-phi[i][j][k]*(2*idx2+2*idy2+2*idz2) +
										rho[i][j][k]/EPS_0 +
										idx2*(phi[i-1][j][k] + phi[i+1][j][k]) +
										idy2*(phi[i][j-1][k]+phi[i][j+1][k]) +
										idz2*(phi[i][j][k-1]+phi[i][j][k+1]))

                       sum += R*R


           L2 = np.sqrt(sum/(ni*nj*nk))
           if L2<tolerance:
              converged=True
              break



    if not converged:
       print("GS failed to converge, L2= " + L2)

	# write_3D(world, rho, "rho_after_solve", world.getTs(), 0);
	# write_3D(world, world.phi, "phi_after_solve", world.getTs(), -2);
    return converged

