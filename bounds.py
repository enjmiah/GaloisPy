from math import *
import sys

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

def tests(n = -1, d = -1, q = -1):
  try:
    Z25 = [
      "00000", "01101", "10110", "11011",
      "00001", "01100", "10111", "11010",
      "00010", "01111", "10100", "11001",
      "00100", "01001", "10010", "11111",
      "01000", "00101", "11110", "10011",
      "10000", "11101", "00110", "01011",
      "10001", "11100", "00111", "01010",
      "00011", "01110", "10101", "11000"
      ]
    if (hasDuplicates(Z25)):
      print("Test failed A.")
      raise Exception
    # print(set(Z25))
    for i in range(0, 2**5):
      str = "{0:b}".format(i)
      str = "0000" + str
      str = str[len(str) - 5:len(str)]
      Z25.insert(len(Z25), str)
      if (not hasDuplicates(Z25)):
        print("Test failed B, found " + str)
        raise Exception
      Z25.remove(str)
    print("All tests passed.")
  except Exception:
    print("Tests failed.")
# tests()

def testAssn5():
  C = [
    "000000", "100101", "010111", "001110", "110010", "101011", "011001", "111100",
    "100000", "000101", "110111", "101110", "010010", "001011", "111001", "011100",
    "010000", "110101", "000111", "011110", "100010", "111011", "001001", "101100",
    "001000", "101101", "011111", "000110", "111010", "100011", "010001", "110100",
    "000100", "100001", "010011", "001010", "110110", "101111", "011101", "111000",
    "000010", "100111", "010101", "001100", "110000", "101001", "011011", "111110",
    "000001", "100100", "010110", "001111", "110011", "101010", "011000", "111101",
    "000011", "100110", "010100", "001101", "110001", "101000", "011010", "111111"
  ]
  print(hasDuplicates(C))
  Z2 = [0, 1]
  for u1 in Z2:
    for u2 in Z2:
      for u3 in Z2:
        for u4 in Z2:
          for u5 in Z2:
            for u6 in Z2:
              str = ("%d%d%d%d%d%d" % (u1, u2, u3, u4, u5, u6))
              if str not in C:
                print(str)
testAssn5()
