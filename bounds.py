from math import *

"""Calculate various bounds for the size of (n, M, d) codes."""

class bcolors:
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'

def choose(n, k):
  return (factorial(n) / (factorial(k) * factorial(n-k)))

def hasDuplicates(list):
  return len(list) != len(set(list))

def hamming_bound(n, d, q = 2):
  try:
    summation = 0
    for i in range(0, floor((d - 1)/2) + 1):
      summation += choose(n, i) * ((q - 1)**i)
    return q**n / summation
  except ZeroDivisionError:
    print(bcolors.WARNING + "Division by zero" + bcolors.ENDC)
    return -1

def singleton_bound(n, d, q = 2):
  return q**(n - d + 1)

def lower_bound(n, d, q = 2):
  try:
    summation = 0
    for m in range(0, d - 1 + 1):
      summation += choose(n, m) * ((q - 1)**m)
    return ceil(q**n / summation)
  except ZeroDivisionError:
    return -1

def bound_info(n, d, q = 2):
  print(u"Singleton bound: A_q ≤ %d" % singleton_bound(n, d, q))
  hamming = hamming_bound(n, d, q)
  if (hamming != -1):
    print(u"Hamming bound:   A_q ≤ %d" % hamming)
  lower = lower_bound(n, d, q)
  if (lower != -1):
    print(u"Lower bound:     A_q ≥ %s" % lower)

