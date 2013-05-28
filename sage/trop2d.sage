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

def amoeba(A, p):
  B = get_B_list(A, True)
  return amoeba_from_B(B, p)
