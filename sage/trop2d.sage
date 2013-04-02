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

    l = map(lambda i: FI(zv[i]), range(len(zv)))
    l = line(l)
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