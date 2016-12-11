# This will simply calculate the real A-discriminant amoeba 
# from the input A matrix as a list. 
# When A\subset Z^n with |A|=n+3.
#
# If you want the real A-discriminant for the support
#
#     [ 6 0 0 0 3 1 ]
# A = [ 0 3 1 6 0 0 ]
#     [ 1 1 1 0 0 0 ]
#
# in the box [-10,10]^2, then you use the command:
# A = [[6,0,0,0,3,1],[0,3,1,6,0,0],[1,1,1,0,0,0]]
# show(amoeba(A, 1/100), xmin=-10, ymin=-10, xmax=10, ymax=10)
#
# The parameters for the function amoeba are as follows:
# A     - The support of the family. Each element of A is a row of the support
#         matrix.
# angle - The smallest angle to differentiate. The smaller the more accurate the
#         resulting graph.
#

var('n')

# This will not run as quickly if the entries of B are 
# irrational. That is, it will leave the irrational part, 
# irrational.
def f_from_B(B):
  linears = map(lambda b: b[0]*cos(n)+b[1]*sin(n), B)
  logs = map(lambda l: l.abs().log(), linears)
  mult = lambda b: map(lambda i: b[i]*logs[i], 
      range(len(logs)))
  f0 = sum(mult(map(lambda b: b[0], B)), 0)
  f1 = sum(mult(map(lambda b: b[1], B)), 0)
  return lambda m: [f0.subs({n:1.0*m}), 
                    f1.subs({n:1.0*m})]

# This gets the zeros of the linear forms of B and converts 
# them to angles (in [0,pi) ) Then it sorts them.
def get_asymptotes(B):
  sgn = lambda b: -1 if b < 0 else 1
  at2 = lambda b: arctan2(sgn(b[0])*b[0], -sgn(b[0])*b[1])
  BB = map(lambda b: at2(b), B)
  BB.sort()
  return BB

# This returns a collection of lines representing the A-discriminant amoeba
# corresponding to B.
# B     - The B matrix.
# count - The number of pieces to break the domain into.
#
# This could be made more accurate if the cusps were found. This is not a
# difficult task, but I have not included it for brevity. I have not found any
# examples where it is noticeably different.
#
# All it does is break the domain into arcs between the asymptotes, then it
# breaks the lines into count many pieces.
def amoeba_from_B(B, count):
  ass = get_asymptotes(B)
  f = f_from_B(B)
  if min(ass) > 1e-6:
    ass = [0] + ass
  if math.pi - max(ass) > 1e-6:
    ass.append(float(math.pi))
  lines = []
  for i in range(len(ass)-1):
    s = ass[i] + 1/count
    e = ass[i+1] - 1/count
    k = ceil((e-s)/math.pi/2*count)+2
    step = (e-s)/(k-1)
    l = map(lambda i: f(s+i*step), range(k))
    lines.append(line(l))
  return sum(lines)

# This gets the B matrix, which is the kernel of \hat{A}.
def get_B_list(A):
  Ah = A + [[1]*len(A[0])]
  Am = matrix(Ah).transpose()
  Bm = Am.integer_kernel().basis_matrix().transpose()
  B = map(lambda b: list(b), list(Bm))
  return B

# This returns the point that poly would be mapped to in the 
# discriminant amoeba.
#
# if log is false then poly is expected to already be the logarithm
# of the coordinates.
def get_point(A, poly, log = False):
  B = get_B_list(A)
  if log:
    pp = lambda q: math.log(abs(q))
  else:
    pp = lambda q: q
  return map(lambda i: sum(map(lambda j: B[j][i]*pp(poly[j]), range(len(poly)))), [0,1])
  
# This gets a polynomial (coefficient vector) that maps to the point, point.
#
# if exp is False then it returns the logarithms of the coefficients.
def get_polynomial(A, point, exp = False):
  B = get_B_list(A)
  BB = matrix(B)
  BBB = BB.transpose()*BB
  fm = matrix(point)*BBB.inverse()*BB.transpose()
  pows = list(list(fm)[0])
  if exp:
    return map(lambda p: math.exp(p), pows)
  return pows

def amoeba(A, angle):
  B = get_B_list(A)
  return amoeba_from_B(B, ceil(2*math.pi/angle+1))
