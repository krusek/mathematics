# This will simply calculate the p-adic A-discriminant amoeba 
# from the input A matrix as a list. 
# When A\subset Z^n with |A|=n+3.
#
# If you want the 3-adic A-discriminant for the support
#
#     [ 6 0 0 0 3 1 ]
# A = [ 0 3 1 6 0 0 ]
#     [ 1 1 1 0 0 0 ]
#
# in the box [-4,25]x[-25,4], then you use the command:
# A = [[6,0,0,0,0,3,1],[0,3,1,6,0,0],[1,1,1,0,0,0]]
# show(amoeba(A, 3), xmin=-4, ymin=-25, xmax=25, ymax=4, aspect_ratio=1)
#
# The parameters for the function amoeba are as follows:
# A     - The support of the family. Each element of A is a row of the support
#         matrix.
# prime - The prime to use for Qp.
#

def amoeba_from_B(B, p):
  Z = filter(lambda b: b[1][0] != 0, enumerate(B))
  Q = Qp(p)
  Bm = matrix(B)
  rval = []
  for z in Z:
    bb = B[z[0]]
    BB = map(lambda b: [b[0]/bb[0], b[1]-b[0]*bb[1]/bb[0]], B)
    QL = map(lambda b: [Q(b[0]).valuation(), Q(b[1]).valuation()], BB)
    FQ = lambda r: map(lambda b: min(b[0]+r,b[1]), QL)
    FI = lambda x: (matrix(FQ(x))*Bm).list()

    zv = map(lambda b: b[1]-b[0], QL)
    zv = filter(lambda zz: zz != Infinity and zz != -Infinity, zv)
    zv = [-100, 100] + list(set(zv))
    zv.sort()
    
    l1 = line([FI(zv[0]), FI(zv[1])], color='green')
    ll = line([FI(zv[-2]), FI(zv[-1])], color='green')
    rval.append(l1)
    rval.append(ll)
    if len(zv) > 2:
      l = line(map(lambda i: FI(zv[i]), range(1, len(zv)-1)))
      rval.append(l)
  return rval

def get_B_list(A, transpose=False):
  Ah = A + [[1]*len(A[0])]
  Am = matrix(Ah).transpose()
  Bm = Am.integer_kernel().basis_matrix()
  if transpose:
    Bm = Bm.transpose()
  B = map(lambda b: list(b), list(Bm))
  return B

# This attempts to get a linear combination of elements of B that give 
# the basis vector e1 or e0, depending on whether index = 0, 1
def get_ei(B, index):
  f = [0]*len(B)
  ii = (index + 1) % 2
  if index == 0:
    det = lambda i, j: B[i][0]*B[j][1]-B[i][1]*B[j][0]
  else:
    det = lambda i, j: B[i][1]*B[j][0]-B[i][0]*B[j][1]
  for i in range(len(B)):
    if B[i][index] == 0: continue
    for j in range(i+1, len(B)):
      if B[j][index] == 0: continue
      dj = det(i,j)
      if dj == 0:
        continue
      if dj == -1:
        i, j = j, i
        dj = 1
      fj = [0]*len(B)
      fj[i] = B[j][ii]
      fj[j] = -B[i][ii]
      if dj == 1:
        return ff
      for k in range(j+1, len(B)):
        if B[k][index] == 0: continue
        dk = det(i, k)
        if dk == 0:
          continue
        if dk == -1:
          i, k = k, i
          dk = 1
        fk = [0]*len(B)
        fk[i] = B[k][ii]
        fk[k] = -B[i][ii]
        if dj == 1:
          return fk
        g = xgcd(dj, dk)
        if g[0] == 1:
          return map(lambda i: fj[i]*g[1]+fk[i]*g[2], range(len(fj)))
  return None  

# This attempts to find a polynomial (coefficient bector) that will map to 
# the point p in the discriminant amoeba.
#
# If p is None, then only the valuations are given.
def get_polynomial(A, point, p = None):
  B = get_B_list(A, True)
  e0 = get_ei(B, 0)
  e1 = get_ei(B, 1)
  if e0 == None or e1 == None:
    print "FAILED"
    return None
  p0 = map(lambda ee: ee*point[0], e0)
  p1 = map(lambda ee: ee*point[1], e1)
  ee = map(lambda i: p0[i]+p1[i], range(len(e0)))
  if p == None:
    return ee
  else:
    return map(lambda e: p^e, ee)

# This returns the point that the input polynomial (coefficient vector) 
# maps to.
#
# if p is None then poly is expected to be a list of powers. That is, the 
# actual polynomial is expected to be [p^poly[i] for i in range(len(poly))]
def get_point(A, poly, p = None):
  B = get_B_list(A, True)
  if p != None:
    Q = Qp(p)
    pp = lambda q: Q(q).valuation()
  else:
    pp = lambda q: q
  return map(lambda i: sum(map(lambda j: B[j][i]*pp(poly[j]), range(len(poly)))), [0,1])

def amoeba(A, p):
  B = get_B_list(A, True)
  return amoeba_from_B(B, p)
