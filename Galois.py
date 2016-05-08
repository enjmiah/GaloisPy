from functools import reduce

"""
GaloisPy
A Python library for computations involving finite fields

@author enjmiah / jerryyin.me
"""

"""
A finite field.  You can create a field with optional argument verbose as True, 
    which will have some methods print to the console the steps it took to get 
    to the answer.  This is useful if you wish to know *how* to get to an answer
    just as much as what the answer is.  You can change the verbose member any
    time after initialization as well.
    
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
        """Returns an equivalent scalar in the field, if possible"""
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

            if x == 0 or y == 0:
                return y
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

    def mult_inverse(self, a):
        """Returns the multiplicative inverse of scalar a"""
        if a == 0:
            raise ZeroDivisionError()
        elif a == 1:
            return 1
        elif self._modular:
            return self._prime_field_mult_inverse(a % self.size, False)
        elif self.size == 4:
            if a == "a":
                return "b"
            elif a == "b":
                return "a"
            else:
                raise TypeError("Value %s was not in GF(4)." % a)
        else:
            raise NotImplementedError();

    def _prime_field_mult_inverse(self, a, verbose=None): # TODO: verbose
        """
        Calculate the inverse of a in the field GF(n) where n is prime, using the
        Euclidean Algorithm
        """
        t = 0
        r = self.size
        newt = 1
        newr = a
        while newr != 0:
            quotient = r // newr
            (t, newt) = (newt, t - quotient * newt)
            (r, newr) = (newr, r - quotient * newr)

        if r > 1:
            return "a is not invertible"
        if t < 0:
            t += self.size

        return t

    def _mult_vec(self, u, v):
        """Do component-wise multiplication on vectors u and v."""
        return list(map(self.mult_scalar, u, v))

    def add_vec(self, u, v):
        """Add two vectors u and v and returns the result."""
        return list(map(self.add_scalar, u, v))

    def scale_vec(self, a, v):
        """Multiplies vector v by scalar a."""
        def f(b):
            return self.mult_scalar(a, b)
        return list(map(f, v))

    def is_lin_indep(self, u, v):
        """Determine whether vectors u and v are linearly independent."""
        for i in self.elements:
            if u == self.scale_vec(i, v):
                return False
        return True

    # TODO: generalize to work with sets of more than two vectors
    def is_lin_indep_set(self, S):
        """Determine whether a set is linearly independent."""
        for u in S:
            for v in S:
                if not self.is_lin_indep(u, v):
                    return False
        return True

    def dot_vec(self, u, v):
        """Return the inner (dot) product of vectors u and v"""
        return reduce(self.add_scalar, self._mult_vec(u, v), 0);

    def is_generator_matrix(self, A):
        """
        Returns true iff A is a generator matrix (has linearly independent rows)
        """
        return self.is_lin_indep(A)
        
    def is_pc_matrix(self, A, B, verbose=None):
        """
        Returns true iff B is a parity-check matrix of A, where A and B are
        arrays of row vectors (also arrays).
        """
        for u in A:
            for v in B:
                if self.dot_vec(u, v) != 0:
                    self._v_print(u, " and ", v, " are not orthogonal.")
                    return False
        return True
    
    def create_pc_matrix(self, G, verbose=None):
        """Returns the parity-check matrix for generator matrix G"""
        # TODO
        pass
        
    def create_gen_matrix(self, H, verbose=None):
        """Returns the generator matrix for parity-check matrix H"""
        # TODO: (remember to delete zero rows)
        pass
    
    def is_standard_form(self, M, type="generator"):
        """
        Using 'generator' or 'g' as type, returns True iff M is a generator
            matrix in standard form.  This is the default behaviour.
            Using 'parity' or 'p' as type, returns True iff M is a parity-check
            matrix in standard form.
        """
        # TODO
        pass
        
    def rref(self, M, verbose=None):
        """
        Converts a matrix M into row echelon form.
            Note: Unlike most of the other functions in this module, rref() is
            in-place
        """
        rows = len(M)
        cols = len(M[0])
        for row in M:
            if len(row) != cols:
                raise ValueError("Matrix not valid, check row lengths")
        for i in range(rows):
            for j in range(cols):
                M[i][j] = self.identity(M[i][j])
        
        def exchange_rows(m1, m2):
            if m1 != m2:
                (M[m1], M[m2]) = (M[m2], M[m1])
        
        def add_row(m1, a, m2):
            M[m2] = self.add_vec(M[m2], self.scale_vec(a, M[m1]))
            
        def pivot_down(m, n):
            pivot = M[m][n]
            if pivot == 0:
                raise Exception("There has been a terrible error in RREF")
            for i in range(m + 1, rows):
                below = M[i][n]
                if below != 0:
                    multiplier = self.negative(
                        self.mult_scalar(self.mult_inverse(pivot), below)
                        )
                    add_row(m, multiplier, i)
        
        def pivot_up(m, n):
            pivot = M[m][n]
            if pivot == 0:
                raise Exception("There has been a terrible error in RREF")
            for i in range(m - 1, -1, -1):
                above = M[i][n]
                if above != 0:
                    multiplier = self.negative(
                        self.mult_scalar(self.mult_inverse(pivot), above)
                        )
                    add_row(m, multiplier, i)
        
        def reduce_row(m, n):
            pivot = M[m][n]
            M[m] = self.scale_vec(self.mult_inverse(pivot), M[m])
        
        # Turn M into an upper triangular matrix
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
            
            self._v_printM(M, verbose)
            
            m += 1; n += 1

        for pivot in pivots:
            pivot_up(*pivot)
            
        self._v_printM(M, verbose)
    
    def rank(self, M):
        """Returns the rank of matrix M"""
        # TODO
        pass
    
    def encode(self, G, c):
        """Encodes codeword c using generator matrix G"""
        # TODO
        pass
        
    def _v_print(self, str, verbose):
        """Prints if verbose is on"""
        verbose = self.verbose if verbose is None else verbose
        if verbose:
            print(str)
    
    def _v_printM(self, M, verbose):
        """Prints a formatted matrix if verbose is on"""
        verbose = self.verbose if verbose is None else verbose
        if verbose:
            print("")
            for row in M:
                print("|", end = " ")
                for el in row:
                    print(el, end = " ")
                print("|")
        

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
