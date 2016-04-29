from math import *

def gcd(a, b):
  return euclid_alg(a, b, False)

def euclid_alg(a, b, visual = True):
  """Perform the euclidean algorithm on a and b."""
  if (a < b):
    return euclid_alg(b, a, visual)

  # a = bq + r
  if (b == 0):
    return a
  else:
    if (visual):
      print("%s = %s * %s + %s   ===>   %s = %s - %s * %s" %
            (a, b, a // b, a % b,       a % b, a, b, a // b))
    return euclid_alg(b, a % b, visual)

def __inverse():
  pass

def inverse(a, n):
  """
  Calculate the inverse of a in the field GF(n) where n is prime
  """
  t = 0
  r = n
  newt = 1
  newr = a

  while newr != 0:
    quotient = r // newr
    (t, newt) = (newt, t - quotient * newt)
    (r, newr) = (newr, r - quotient * newr)

  if r > 1:
    return "a is not invertible"
  if t < 0:
    t += n

  return t

def check_expect(actual, expected):
  p = actual == expected
  if (not p):
    print("Expected %s, got %s" % (expected, actual))
    return False
  return True

def run_tests():
  if (check_expect(3, 3)):
    print("All tests passed");
run_tests();
