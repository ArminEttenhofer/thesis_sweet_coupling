#!/usr/bin/env python3
#
# Author: Boriskova E.Y. <boriskova.eyu@gmail.com>
# Date 2023-12-29
# 
import numpy as np

from mule.sdc.nodes import NodesGenerator
from mule.SWEETFileDict import SWEETFileDict

try:
    from .fd_weights_explicit import get_fd_stencil
except:
    from fd_weights_explicit import get_fd_stencil

def getSetup(
    nNodes:int=3, nIter:int=-1, nodeDistr:str='LEGENDRE'
    )-> SWEETFileDict:
    """
    Generate SWEETFileDict for one given ETD SDC setup

    Parameters
    ----------
    nNodes : int, optional
        Number of nodes (time substeps!).
    nIter : int, optional
        Number of iterations (sweeps).
    nodeDistr : str, optional
        Node distribution ('EQUID', 'LEGENDRE', 'CHEBY-1', 'CHEBY-2', 'CHEBY-3', 'CHEBY-4').

    Returns
    -------
    out : SWEETFileDict
        The SDC parameters and coefficients, with the following keys :
            
        - tauNodes : np.ndarray(nNodes+2)
            Time substeps between (0, 1).
        - tau : np.ndarray(nNodes+2)
            Time substeps including 0 and 1.
        - deltaTau : np.ndarray(nNodes+1)
            Time substeps, node-to-node distance.
        - nIters : str
            Number of iterations (sweeps).
        - A : np.ndarray(nNodes+1, nNodes+2, nNodes+2)
            Finite difference coefficients a(i)jl.
        - idString : str
            Unique string identifier.
    """
    assert nNodes >= 0, "nNodes needs to be 0 or positive"
    assert nNodes < 20, "Only up to 20 substeps are supported"

    # Get inner nodes and aijl
    nodes, tau, deltaTau, a = genWeights(nNodes, nodeDistr)
    assert a.shape[-2] == nNodes + 2, "a(i)jl size inconsistent with node number"

    # Save M as Aijl dimention
    M = nNodes + 2

    if nIter == -1:
        # default case: optimal sweep amount
        # for ord = min(M, nIter + 1) [= min(N, M+1) Buvoli et al]
        nIter = M - 1
    else:
        assert nIter >= 0, "nIter needs to be 0 or positive"
    
    # Save nNodes + 2 for total amount, including 0. and 1.
    idString = f'M{nIter}_{nodeDistr}_N{nNodes+2}'

    out = SWEETFileDict(initDict={
        # Necessary parameters
        "tauNodes": nodes,
        "tau": tau,
        "deltaTau": deltaTau,
        "A": a,
        "nIter": nIter,
        "postSweep": "LASTNODE",
        "preSweep": "ZEROS",
        # SDC description variables
        "idString": idString,
    })
    return out
    

def genWeights(nNodes, distr):
    """
    Generate the nodes and aijl matrices

    Parameters
    ----------
    nNodes : int
        Number of nodes.
    distr : str
        Node distribution. Can be selected from :
    - LEGENDRE : nodes from the Legendre polynomials
    - EQUID : equidistant nodes distribution
    - CHEBY-{1,2,3,4} : nodes from the Chebyshev polynomial (1st to 4th kind)
    
    Returns
    -------
    nodes : np.ndarray(nNodes+2)
        Time substeps between (0, 1) following distr.
    tau : np.ndarray(nNodes+2)
        Time substeps including 0 and 1.
    deltaTau : np.ndarray(nNodes+1)
        Time substeps, node-to-node distance.
    a : np.ndarray(nNodes+1, nNodes+2, nNodes+2)
        Finite difference coefficients a(i)jl.
    """

    # Generate nodes [0, 1]
    nodes = NodesGenerator(node_type=distr).getNodes(nNodes + 2)
    nodes += 1
    nodes /= 2
    np.round(nodes, 14, out=nodes)

    tau = nodes

    deltaTau = np.zeros(nNodes+2)
    deltaTau[0] = nodes[0]
    deltaTau[1:] = nodes[1:] - nodes[:-1]

    M = nNodes + 2

    # Generate 3D matrix a(i)jl
    a = np.zeros((M - 1, M, M))

    for j in range(1, M):
        q = (tau - tau[j-1]) / deltaTau[j]
        for l in range(0, M):
            a[j-1, :, l] = get_fd_stencil(l, 0, q)

    return nodes, tau, deltaTau, a

    
def getSetupFromString(idString:str):
    try:
        params = idString.split('_')
        nNodes, nIter, nodeDistr = params
        
        assert nNodes.startswith('M')
        # idString containes nNodes + 2
        nNodes = int(nNodes[1:]) - 2

        assert nIter.startswith('N')
        nIter = int(nIter[1:])

        return getSetup(
            nNodes=nNodes, nodeDistr=nodeDistr, nIter=nIter)

    except Exception:
        raise ValueError(f'{idString} is not a valid idString for ETDSDC, eg must be M4_CHEBY-1_N5')

