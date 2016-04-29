from math import *

def inverse(a, n):
  """
  Calculate the inverse of a in the field GF(n) where n is prime, using the
  Euclidean Algorithm
  """
  t = 0
  r = n
  newt = 1
  newr = a

  while newr != 0:
    quotient = floor(r / newr)
    (t, newt) = (newt, t - quotient * newt)
    (r, newr) = (newr, r - quotient * newr)

  if r > 1:
    return "a is not invertible"
  if t < 0:
    t += n

  return t
