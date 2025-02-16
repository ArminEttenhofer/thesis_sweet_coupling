#! /usr/bin/env python3

"""
This is a prototyping version of a banded diagonal solver which can be used in SWEET

A refers to the full dense matrix
a refers to the banded-only storage
"""

import numpy as np
import scipy as sp
import sys



def solver_Aa(A, a, b, mat_size, halo_size, debug=False):
    A = A.copy()
    a = a.copy()
    b = b.copy()

    if debug:
        print(f"A: {A}")
        print(f"a: {a}")
        print(f"b: {b}")

    # Eliminate lower diagonal matrix
    for ipiv in range(mat_size-1):
        # Pivot element given by A[ipiv,ipiv]

        # Set col below A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(ipiv+1, min(ipiv+halo_size+1, mat_size)):

            if np.isclose(A[ipiv,ipiv], 0):
                print(f"A: {A}")
                raise Exception("A[ipiv,ipiv] == 0")

            if np.isclose(a[ipiv,halo_size+1], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,ipiv] == 0")

            pivA = A[jrow,ipiv]/A[ipiv,ipiv]
            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]
            assert np.isclose(pivA, piva)

            # Iterate over columns
            for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
                A[jrow,jcol] -= A[ipiv,jcol]*pivA
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva

                assert np.isclose(A[jrow,jcol], a[jrow,jcol-jrow+halo_size])

            b[jrow] -= b[ipiv]*pivA



    # Eliminate upper diagonal matrix
    for ipiv in range(mat_size-1, -1, -1):
        # Pivot element given by A[ipiv,ipiv]

        # Set col above A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(max(ipiv-halo_size, 0), ipiv):

            if np.isclose(A[ipiv,ipiv], 0):
                print(f"A: {A}")
                raise Exception("A[ipiv,ipiv] == 0")

            if np.isclose(a[ipiv,halo_size], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,halo_size] == 0")

            pivA = A[jrow,ipiv]/A[ipiv,ipiv]
            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]
            assert np.isclose(pivA, piva)

            # Iterate over columns
            for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
                A[jrow,jcol] -= A[ipiv,jcol]*pivA
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva
                assert np.isclose(A[jrow,jcol], a[jrow,jcol-jrow+halo_size])

            b[jrow] -= b[ipiv]*pivA

        b[ipiv] /= A[ipiv,ipiv]
        A[ipiv,ipiv] = 1.0

    return b


def solver_A(A, b, mat_size, halo_size, debug=False):
    A = A.copy()
    b = b.copy()

    if debug:
        print(f"A: {A}")
        print(f"b: {b}")

    # Eliminate lower diagonal matrix
    for ipiv in range(mat_size-1):
        # Pivot element given by A[ipiv,ipiv]

        # Set col below A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(ipiv+1,min(ipiv+halo_size+1, mat_size)):

            pivA = A[jrow,ipiv]/A[ipiv,ipiv]

            # Iterate over columns
            for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
                A[jrow,jcol] -= A[ipiv,jcol]*pivA

            b[jrow] -= b[ipiv]*pivA


    # Eliminate upper diagonal matrix
    for ipiv in range(mat_size-1, -1, -1):
        # Pivot element given by A[ipiv,ipiv]

        # Set col above A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(max(ipiv-halo_size, 0), ipiv):

            if np.isclose(A[ipiv,ipiv], 0):
                print(f"A: {A}")
                raise Exception("A[ipiv,ipiv] == 0")

            pivA = A[jrow,ipiv]/A[ipiv,ipiv]

            # Iterate over columns
            for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
                A[jrow,jcol] -= A[ipiv,jcol]*pivA

            b[jrow] -= b[ipiv]*pivA

        b[ipiv] /= A[ipiv,ipiv]
        A[ipiv,ipiv] = 1.0

    return b



def solver_a(a, b, mat_size, halo_size, debug=False):
    a = a.copy()
    b = b.copy()

    if debug:
        print(f"a: {a}")
        print(f"b: {b}")

    # Eliminate lower diagonal matrix
    for ipiv in range(mat_size-1):
        # Pivot element given by A[ipiv,ipiv]

        # Set col below A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(ipiv+1,min(ipiv+halo_size+1, mat_size)):
            if np.isclose(a[ipiv,halo_size+1], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,ipiv] == 0")

            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]

            # Iterate over columns
            for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva

            b[jrow] -= b[ipiv]*piva



    # Eliminate upper diagonal matrix
    for ipiv in range(mat_size-1, -1, -1):
        # Pivot element given by a[ipiv,ipiv]

        # Set col above A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(max(ipiv-halo_size, 0), ipiv):

            if np.isclose(a[ipiv,halo_size], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,halo_size] == 0")

            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]

            # Iterate over columns
            for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva

            b[jrow] -= b[ipiv]*piva

        b[ipiv] /= a[ipiv,halo_size]
        a[ipiv,halo_size] = 1.0

    return b


"""
Optimized version with some implicit assumptions.
The resulting banded storage "a" doesn't store the diagonal in the end
"""
def solver_a_opti(a, b, mat_size, halo_size, debug=False):
    a = a.copy()
    b = b.copy()

    if debug:
        print(f"a: {a}")
        print(f"b: {b}")

    # Eliminate lower diagonal matrix
    for ipiv in range(mat_size-1):
        # Pivot element given by A[ipiv,ipiv]

        # Process rows below below A[i,i] to zero
        for jrow in range(ipiv+1,min(ipiv+halo_size+1, mat_size)):
            if np.isclose(a[ipiv,halo_size+1], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,ipiv] == 0")

            assert ipiv-jrow+halo_size >= 0
            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]

            # Iterate over columns
            # Opti: Removed ipiv starting index
            #for jcol in range(ipiv+1, min(ipiv+halo_size+1, mat_size)):
            for jcol in range(ipiv, min(ipiv+halo_size+1, mat_size)):
                assert jcol-jrow+halo_size >= 0
                assert jcol-ipiv+halo_size >= 0
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva

            # Opti: Processed here, but then eliminated
            #a[jrow,ipiv-jrow+halo_size] = 0

            b[jrow] -= b[ipiv]*piva



    # Eliminate upper diagonal matrix
    for ipiv in range(mat_size-1, -1, -1):
        # Pivot element given by a[ipiv,ipiv]

        # Set col above A[i,i] to zero
        # => Iterate over rows below
        for jrow in range(max(ipiv-halo_size, 0), ipiv):

            if np.isclose(a[ipiv,halo_size], 0):
                print(f"a: {a}")
                raise Exception("a[ipiv,halo_size] == 0")

            piva = a[jrow,ipiv-jrow+halo_size]/a[ipiv,halo_size]

            # Iterate over columns
            for jcol in range(max(ipiv-halo_size+1, 0), min(ipiv+1, mat_size)):
                assert jcol-jrow+halo_size >= 0
                assert jcol-ipiv+halo_size >= 0
                a[jrow,jcol-jrow+halo_size] -= a[ipiv,jcol-ipiv+halo_size]*piva

            b[jrow] -= b[ipiv]*piva

        b[ipiv] /= a[ipiv,halo_size]
        # No required
        #a[ipiv,halo_size] = 1.0

    return b



np.set_printoptions(precision=4)

def test():

    #debug = True
    debug = False

    for mat_size in range(6, 16):
        for halo_size in range(0, mat_size):


            print("*"*80)
            print("*"*80)
            print("*"*80)
            print(f"mat_size: {mat_size}")
            print(f"halo_size: {halo_size}")
            print("")

            num_diagonals = halo_size*2+1

            # Matrix
            A = np.zeros((mat_size, mat_size))

            # Banded storage
            a = np.zeros((mat_size, num_diagonals))

            b = np.zeros(mat_size)

            c = 1
            for i in range(mat_size):
                k = i
                for j in range(num_diagonals):
                    l = j-halo_size+k

                    c+=1

                    k_ = np.pi*k/mat_size
                    l_ = np.pi*l/mat_size
                    value = np.cos(np.pi*(k+0.5)*l/mat_size)

                    if 1:
                        value += mat_size*(1.0/(1+np.abs(k-l)))**2


                    ja = j+i-halo_size

                    if ja < 0:
                        value = np.inf
                    if ja >= mat_size:
                        value = np.inf

                    a[i][j] = value

                    if np.isinf(value):
                        continue

                    A[i][ja] = value

                b[i] = np.cos(k/mat_size)


            if debug:
                print(f"a: {a}")
                print(f"A: {A}")
                print(f"b: {b}")

            a_copy = a.copy()
            A_copy = A.copy()
            b_copy = b.copy()


            # Solve using numpy
            x = np.linalg.solve(A_copy, b_copy)
            print(f"x: {x}")

            for i in range(4):
                print("*"*80)

                if i == 0:
                    print("solver_A(...)")
                    x_my = solver_A(A, b, mat_size, halo_size, debug)

                elif i == 1:
                    print("solver_a(...)")
                    x_my = solver_a(a, b, mat_size, halo_size, debug)

                elif i == 2:
                    print("solver_Aa(...)")
                    x_my = solver_Aa(A, a, b, mat_size, halo_size, debug)

                elif i == 3:
                    print("solver_a_opti(...)")
                    x_my = solver_a_opti(a, b, mat_size, halo_size, debug)


                res = np.max(np.abs(b_copy-np.dot(A_copy,x_my)))
                print(f"res: {res}")

                error = np.max(np.abs(x-x_my))
                print(f"error: {error}")

                assert error < 1e-12


test()


print("*"*80)
print("* FIN")
print("*"*80)
