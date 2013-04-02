# This will simply calculate the real A-discriminant amoeba 
# from the input A matrix as a list.

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

def get_B_list(A):
  Ah = A + [[1]*len(A[0])]
  Am = matrix(Ah).transpose()
  Bm = Am.integer_kernel().basis_matrix().transpose()
  B = map(lambda b: list(b), list(Bm))
  return B

def amoeba(A, angle):
  B = get_B_list(A)
  return amoeba_from_B(B, ceil(2*math.pi/angle+1))
