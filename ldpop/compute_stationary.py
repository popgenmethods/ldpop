
import numpy, sys, logging, math, time
from numpy import zeros, array
import scipy.sparse as sp
from scipy.linalg import norm

def assertValidProbs(x):
    assert numpy.all(x >= 0) and abs(1.0 - numpy.sum(x)) < 1e-10, numpy.sum(x)

''' Q must be tridiagonal and have no absorbing states'''
def stationary1d_tridiagonal(Q):
#     logging.info("Starting Tridiagonal Stationary")
    assert Q.shape[0] == Q.shape[1]
    n = Q.shape[0] - 1
    # check Q is tridiagonal
    for i in range(n+1):
        for j in range(n+1):
            if (abs(i-j) <= 1 and Q[i,j] == 0.0) or (abs(i-j) > 1 and Q[i,j] != 0.0):
                raise Exception("at Q[%d,%d]: Q must be strictly tridiagonal" % (i,j))
    
    ret = [0.0]
    for i in range(n):
        # pi[i] * Q[i,i+1] = pi[i+1] * Q[i+1,i]
        ret.append(ret[i] + math.log(Q[i,i+1]) - math.log(Q[i+1,i]))
    ret = numpy.array(ret)
    ret = ret - numpy.max(ret)
    ret = numpy.exp(ret)
    ret = ret / numpy.sum(ret)
#     logging.info("Done!")
    return ret

def stationary(Q, init=None, norm_order=1, epsilon=1e-8):
#     logging.info( "Starting Stationary")
    start = time.time()
    
    size = Q.shape[0]
    assert size == Q.shape[1]
    # Compute a suitable stochastic matrix by means of uniformization
#     l = min(Q.values())*1.001  # avoid periodicity, see the book of Bolch et al.
    l = Q.min() * 1.001
    P = sp.eye(size, size) - Q/l
    # compute Pi
    P =  P.tocsr()
#     pi = zeros(size);  pi1 = zeros(size)
#     pi[0] = 1;
#     n = norm(pi - pi1,1)
    if init is not None:
        assertValidProbs(init)
        pi = init
    else:
        pi = array([1 / float(size) for i in range(size)])
    n = float("inf")
    
    i = 0;
    #maxIts = 1e5
    while n > epsilon:# and i < maxIts:
#     while n > 1e-10 and i < 1e6:
#     while i < 1e6:
        pi1 = pi*P
        pi = pi1*P   # avoid copying pi1 to pi
        
        not_zero = numpy.logical_and(pi != 0., pi1 != 0.)
        if numpy.all((pi == 0.) == (pi1 == 0.)):
#             print str(pi)
#             print str(pi1)
            n = norm(numpy.log(pi[not_zero] / pi1[not_zero]),norm_order); i += 1
        else:
            n = float("inf")
#         n = norm(pi - pi1,1); i += 1
    end = time.time()
    logging.info( "Computed stationary with L-%f stopping, epsilon=%g, in %d iterations, %f seconds " % (norm_order, epsilon, i, end - start))
#     logging.info( "Done with Stationary" )
    return pi


'''
#I stole this from SO: http://stackoverflow.com/questions/21308848/markov-chain-stationary-distributions-with-scipy-sparse
def stationary(Q, direct = False, tol=1e-35):
    Qcsr = Q.tocsr()
    logging.info("Computing Stationary.")
    n = Qcsr.shape[0]
    Qt = Qcsr.T[1:,:]
    A = (sp.vstack([numpy.ones(n), Qt])).tocsr()
    rhs = numpy.zeros((n,))
    rhs[0] = 1

    if direct:
        # Requires that the solution is unique
        toReturn = array(spsolve(A,rhs))
        return toReturn
    else:
        # GMRES does not care whether the solution is unique or not, it
        # will pick the first one it finds in the Krylov subspace
        #print "Starting GMRES"
        #if you want to use the following preconditioner, determine a suitably
        #easier to solve matrix Aapprox
        logging.info("Starting Incomplete LU Decomp")
        M_init = spilu(A.tocsc())
        M_x = lambda x: M_init.solve(x)
        M = LinearOperator((n, n), M_x)
        logging.info("Done!")
        logging.info("Starting GMRES")
        p, info = gmres(A, rhs, M=M, tol=tol)
        logging.info("Done!")
        if info != 0:
            logging.info("GMRES did not converge!")
            exit(-1)
        #print "DONE!"
        return array(p)


'''
'''
def stationary(Q):
    eig1 = sp.linalg.eigs(Q, 1)
    return eig1
''' 
'''
def stationary(Q, epsilon=1e-10):
    t1 = 20.0
    t2 = 25.0
    assert t1 < t2
    size = Q.shape[0]
    initial = array([1 / float(size) for i in range(size)])
    ends = expm_multiply(Q.transpose(), initial, start=t1, stop=t2, num=2)
    end1 = ends[0,:]; end2 = ends[1,:]
    n = norm(end1 - end2, 1)
    assert n < epsilon
    return end2
'''
