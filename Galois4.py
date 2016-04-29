from functools import reduce

"""
   Galois4
   A Python library for computations involving finite fields

   @author enjmiah / jerryyin.me
"""


"""A finite field."""
class GF:
   def __init__(self, size):
      self.size = size
      self.modular = _is_prime(size)

# SCALARS
GF4 = GF(4)
GF4.a = 2
GF4.b = 3
a = 2
b = 3

# TODO: generalize to other fields
def mult_scalar(x, y):
   """Multiply two scalars x and y and return the result."""
   if (x == 0 or y == 0):
      return 0
   elif x == 1:
      return y
   elif y == 1:
      return x

   if (x == GF4.a and y == GF4.a):
      return GF4.b
   elif (x == GF4.a and y == GF4.b) or (x == GF4.b and y == GF4.a):
      return 1
   elif (x == GF4.b and y == GF4.b):
      return GF4.a

   raise ValueError("Passed values %d, %d not in GF(4)." % (x, y))

def add_scalar(x, y):
   """Add two scalars x and y and return the result."""
   if x == y:
      return 0
   if y < x:
      t = y
      y = x
      x = t

   if x == 0:
      return y
   elif (x == a and y == b):
      return 1
   elif (x == 1 and y == b):
      return GF4.a
   elif (x == 1 and y == a):
      return GF4.b

def exp_scalar(a, n):
   """Returns a**n over the appropriate finite field."""
   if n == 0:
      return 1
   elif n == 1:
      return a
   else:
      return mult_scalar(a, exp_scalar(a, n-1))

def add_inverse(a):
   """Returns the additive inverse of scalar a"""
   return a

# TODO
def mult_inverse(a):
   """Returns the multiplicative inverse of scalar a"""
   pass

def __mult_vec(u, v):
   """Do component-wise multiplication on vectors u and v."""
   return map(mult_scalar, u, v)

def add_vec(u, v):
   """Add two vectors u and v and returns the result."""
   return map(add_scalar, u, v)

def scale_vec(a, v):
   """Multiplies vector v by scalar a."""
   def f(b):
      return mult_scalar(a, b)
   return map(f, v)

def is_lin_indep(u, v):
   """Determine whether vectors u and v are linearly independent."""
   for i in [0, 1, a, b]:
      if u == scale_vec(i, v):
         return False
   return True

# TODO: generalize to work with sets of more than two vectors
def is_lin_indep_set(S):
   """Determine whether a set is linearly independent."""
   for u in S:
      for v in S:
         if (not is_lin_indep(u, v)):
            return False
   return True

def dot_vec(u, v):
   """Return the inner (dot) product of vectors u and v"""
   return reduce(add_scalar, __mult_vec(u, v), 0);

def is_pc_matrix(A, B):
   """
   Returns true if B is a parity-check matrix of A, where A and B are arrays of
      row vectors (also arrays).
   """
   for u in A:
      for v in B:
         if (dot_vec(u, v) != 0):
            return False
   return True

# This function is only used for sizes of finite fields, which tend to be not too large,
# so the O(sqrt(n)) performance is pretty ok
def _is_prime(n):
   """
   Return True iff n is prime
   Source: http://stackoverflow.com/a/1801446
   """
   if n <= 1:
      return False
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
