from __future__ import print_function
"""
GaloisPy
A Python library for computations involving finite fields
"""
import copy
from functools import reduce

__author__ = "enjmiah / jerryyin.me"


"""
A finite field.  You can create a field with optional argument verbose as True, 
which will have some methods print to the console the steps it took to get to
the answer.  This is useful if you wish to know *how* to get to an answer just
as much as what the answer is.  You can change the verbose member any time
after initialization as well.
    
Usage:
>>> GF2 = GF(2)
>>> GF4 = GF(4, verbose=True)
>>> GF7 = GF(7)
>>> GF7.verbose = True
"""
class GF:
    size = 0
    verbose = False
    elements = []
    
    def __init__(self, size, verbose=False):
        self.size = size
        self.verbose = verbose
        self._modular = _is_prime(size)
        if self._modular == True:
            self.elements = list(range(0, self.size))
        elif self.size == 4:
            self.elements = [0, 1, "a", "b"]
        else:
            raise NotImplementedError()
        
    def identity(self, x):
        """Returns an equivalent scalar or vector in the field, if possible"""
        if isinstance(x, list):
            return [self.identity(a) for a in x]
        else:
            if self._modular:
                return x % self.size
            elif self.size == 4:
                if x not in self.elements:
                    raise ValueError("Passed value not in GF(4).")
                return x
            else:
                raise NotImplementedError();   
    
    def mult_scalar(self, x, y):
        """Multiply two scalars x and y and return the result."""
        if self._modular:
            return (x * y) % self.size
        elif self.size == 4:
            if x == 0 or y == 0:
                return 0
            elif x == 1:
                return y
            elif y == 1:
                return x

            if x == "a" and y == "a":
                return "b"
            elif (x == "a" and y == "b") or (x == "b" and y == "a"):
                return 1
            elif x == "b" and y == "b":
                return "a"

            raise ValueError("Passed values %s, %s not in GF(4)." % (x, y))
        else:
            raise NotImplementedError()
    
    def add(self, x, y):
        """Add two scalars or vectors x and y and return the result."""
        if isinstance(x, list) and isinstance(y, list):
            if len(x) != len(y):
                raise ValueError("Tried to add vectors of different dimensions")
            return self.add_vec(x, y)
        else:
            return self.add_scalar(x, y)
    
    def add_scalar(self, x, y):
        """Add two scalars x and y and return the result."""
        if self._modular:
            return (x + y) % self.size
        elif self.size == 4:
            if x == y:
                return 0

            if x == 0:
                return y
            elif y == 0:
                return x
            elif (x == "a" and y == "b") or (x == "b" and y == "a"):
                return 1
            elif (x == 1 and y == "b") or (x == "b" and y == 1):
                return "a"
            elif (x == 1 and y == "a") or (x == "a" and y == 1):
                return "b"
        else:
            raise NotImplementedError()

    def exp_scalar(self, a, n):
        """Returns a**n over the appropriate finite field."""
        if self._modular:
            return (a**n) % self.size
        elif self.size == 4:
            if n == 0:
                return 1
            elif n == 1:
                return a
            else:
                return self.mult_scalar(a, self.exp_scalar(a, n-1))
        else:
            raise NotImplementedError()

    def add_inverse(self, x):
        """Returns the additive inverse of scalar or vector x"""
        if isinstance(x, list):
            return [self.add_inverse(a) for a in x]
        else:
            if self._modular:
                return self.size - x
            elif self.size == 4:
                return x
            else:
                raise NotImplementedError()
    
    def negative(self, x):
        """
        Returns the additive inverse of scalar or vector x.
        See add_inverse()
        """
        return self.add_inverse(x)

    def mult_inverse(self, a, verbose=None):
        """Returns the multiplicative inverse of scalar a"""
        if a == 0:
            raise ZeroDivisionError()
        elif a == 1:
            return 1
        elif self._modular:
            return self._prime_field_mult_inverse(a % self.size, verbose)
        elif self.size == 4:
            if a == "a":
                return "b"
            elif a == "b":
                return "a"
            else:
                raise TypeError("Value %s was not in GF(4)." % a)
        else:
            raise NotImplementedError();
    
    def _prime_field_mult_inverse(self, a, verbose=None):
        """
        Calculate the inverse of a in the field GF(n) where n is prime, using
        the Euclidean Algorithm in reverse.
        """
        if a == 0:
            raise ZeroDivisionError("0 is not invertible.")
        
        self._v_print("", verbose)
        t = 0
        r = self.size
        newt = 1
        newr = a
        MSG1 = "%s = %s * %s + %s"
        MSG2 = "   ==>   %s = %s - %s * %s"
        MSG3 = "   ==>   t = %s - %s * %s = %s"
        while newr != 0:
            quotient = r // newr
            remain = (r % newr if newr != 0 else 1)
            self._v_print(MSG1 % (r, r // newr, newr, remain), verbose, end="")
            self._v_print(MSG2 % (remain, r, r // newr, newr), verbose, end="")
            if (r - quotient * newr) != 0:
                self._v_print(MSG3 % (t, quotient, newt, t - quotient * newt),
                              verbose)
            else:
                self._v_print("", verbose)
            (t, newt) = (newt, t - quotient * newt)
            (r, newr) = (newr, r - quotient * newr)

        if t < 0:
            t += self.size
            MSG = u"%s mod %s = %s"
            self._v_print(MSG % (t - self.size, self.size, t), verbose)
        self._v_print("Multiplicative inverse is %s"%(t), verbose)
        
        return t

    def add_vec(self, u, v):
        """Add two vectors u and v and returns the result."""
        return [self.add_scalar(a, b) for a,b in zip(u, v)]

    def scale_vec(self, a, v):
        """Multiplies vector v by scalar a."""
        return [self.mult_scalar(a, b) for b in v]

    def is_lin_indep(self, S):
        """Determine whether a set is linearly independent."""
        return self.rank(S) == len(S)

    def dot_vec(self, u, v):
        """Return the inner (dot) product of vectors u and v"""
        if len(u) != len(v):
            raise ValueError("Vectors must be same length.")
        return reduce(self.add_scalar, self._mult_vec(u, v), 0);
        
    def _mult_vec(self, u, v):
        return [self.mult_scalar(a, b) for a,b in zip(u, v)]

    def is_generator_matrix(self, M):
        """
        Returns true iff M is a generator matrix (has linearly independent rows)
        """
        return self.is_lin_indep(M)
        
    def is_pc_matrix(self, A, B, verbose=None):
        """
        Returns true iff B and A are parity-check matrices of each other.
        """
        for u in A:
            for v in B:
                if self.dot_vec(u, v) != 0:
                    MSG = "%s and %s are not orthogonal."
                    self._v_print(MSG % (u, v), verbose)
                    return False
        return True
    
    def create_pc_matrix(self, G, verbose=None):
        """Returns a parity-check matrix for generator matrix G"""
        rows = len(G)
        cols = len(G[0])
        for row in G:
            if len(row) != cols:
                raise ValueError("Matrix not valid, check row lengths.")
        if not self.is_standard_form(G, 'g'):
            G = self.rref(G)
        if not self.is_generator_matrix(G):
            raise ValueError("Passed matrix was not a generator matrix.")
        
        H = self._tranpose(G)
        H = self.negative(H[rows:])
        offset = rows # rows == len(H[0])
        for i in range(len(H)):
            for j in range(len(H)):
                H[i].insert(j+offset, 1 if i == j else 0)
        return H
    
    def is_standard_form(self, M, type="generator"):
        """
        Using 'generator' or 'g' as type, returns True iff M is a generator
        matrix in standard form.  This is the default behaviour.
        Using 'parity' or 'p' as type, returns True iff M is a parity-check
        matrix in standard form.
        """
        rows = len(M)
        cols = len(M[0])
        if rows > cols:
            return False
        for row in M:
            if len(row) != cols:
                return False
        
        if type == "generator" or type == "g":
            for i in range(rows):
                for j in range(rows):
                    if (i == j and M[i][j] != 1) or (i != j and M[i][j] != 0):
                        return False
            return True
        elif type == "parity" or type == "p":
            for i in range(rows):
                for j in range(rows):
                    offset = cols - rows
                    if (i == j and M[i][j+offset] != 1):
                        return False
                    if (i != j and M[i][j+offset] != 0):
                        return False
            return True
        else:
            raise ValueError("type argument must be either 'g', 'generator', 'p', or 'parity'")
        
    def rref(self, M, verbose=None):
        """
        Returns a reduced row echelon form matrix that is row equivalent to
        matrix M.
        """
        rows = len(M)
        cols = len(M[0])
        for row in M:
            if len(row) != cols:
                raise ValueError("Matrix not valid, check row lengths")
        M = copy.deepcopy(M)
        M = self.identity(M)
        
        def exchange_rows(m1, m2):
            if m1 != m2:
                (M[m1], M[m2]) = (M[m2], M[m1])
                MSG = "Exchange rows %s and %s."
                self._v_printM(M, MSG % (m1+1,m2+1), verbose)
        
        def add_row(m1, a, m2):
            M[m2] = self.add(M[m2], self.scale_vec(a, M[m1]))
            MSG = "Added %s times row %s to row %s."
            self._v_printM(M, MSG % (a,m1+1,m2+1), verbose)
            
        def pivot_down(m, n):
            pivot = M[m][n]
            if pivot == 0:
                raise Exception("There has been a terrible error in RREF")
            first = True
            for i in range(m + 1, rows):
                if first:
                    MSG = "PLAN: Pivot down from position (%s, %s)" 
                    self._v_printM(M, MSG % (m+1, n+1), verbose)
                    first = False
                below = M[i][n]
                if below != 0:
                    multiplier = self.negative(
                        self.mult_scalar(self.mult_inverse(pivot, False), below)
                        )
                    add_row(m, multiplier, i)
        
        def pivot_up(m, n):
            pivot = M[m][n]
            if pivot == 0:
                raise Exception("There has been a terrible error in RREF")
            first = True
            for i in range(m - 1, -1, -1):
                if first:
                    MSG = "PLAN: Pivot up from position (%s, %s)"
                    self._v_printM(M, MSG%(m+1, n+1), verbose)
                    first = False
                above = M[i][n]
                if above != 0:
                    multiplier = self.negative(
                        self.mult_scalar(self.mult_inverse(pivot, False), above)
                        )
                    add_row(m, multiplier, i)
        
        def reduce_row(m, n):
            pivot = M[m][n]
            scale = self.mult_inverse(pivot, False)
            M[m] = self.scale_vec(scale, M[m])
            self._v_printM(M, "Scale row %s by %s." % (m+1, scale), verbose)
        
        self._v_printM(M, "Original matrix.", verbose)
        pivots = []
        m = 0; n = 0
        while (m < rows and n < cols):
            # Find a pivot and ensure pivot is at M[m][n]
            success = False
            for i in range(m, rows):
                if M[i][n] != 0:
                    exchange_rows(m, i)
                    success = True
                    break
            if not success:
                # No pivot in column n
                n += 1
                continue
            
            pivot_down(m, n)      # Pivot down
            reduce_row(m, n)      # Scale so that pivot == 0
            pivots.append((m, n)) # Remember pivot
            
            m += 1; n += 1

        for pivot in pivots:
            pivot_up(*pivot)
            
        return M
    
    def rank(self, M):
        """Returns the rank (dimension of rowspace) of matrix M"""
        Mref = self.rref(M)
        return len([v for v in Mref if any(v)])
    
    def encode(self, G, w):
        """
        Encodes word w using generator matrix G and returns the result.
        
        Keyword arguments:
        w -- word of length k (a list of elements in field)
        G -- a k x n matrix
        """
        if len(w) != len(G):
            raise ValueError("Input word is wrong length.")
        G = copy.deepcopy(G)
        for i in range(len(G)):
            G[i] = self.scale_vec(w[i], G[i])
        return reduce(self.add, G)
    
    def _tranpose(self, M):
        """Returns transpose of matrix M"""
        return list(map(list, zip(*M)))
        
    def _v_print(self, str, verbose, end="\n"):
        """Prints if verbose is on"""
        verbose = self.verbose if verbose is None else verbose
        if verbose:
            print(str, end=end)
    
    def _v_printM(self, M, str, verbose):
        """Prints a formatted matrix if verbose is on"""
        verbose = self.verbose if verbose is None else verbose
        if verbose:
            first = True
            print("")
            for row in M:
                print("|", end = " ")
                for el in row:
                    print(el, end = " ")
                print("|", end = "")
                if first:
                    print("   "+str, end = "")
                    first = False
                print("")
        

# This function is only used for sizes of finite fields, which tend to be not
# too large, so the O(sqrt(n)) performance is ok
def _is_prime(n):
    """
    Return True iff n is prime
    Source: http://stackoverflow.com/a/1801446
    """
    if n == 2:
        return True
    if n == 3:
        return True
    if n % 2 == 0:
        return False
    if n % 3 == 0:
        return False

    i = 5
    w = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += w
        w = 6 - w
    return True
