def getlist(n, i, l, ex=[0]):
  for j in range(i, n):
    if ex.count(j) > 0:
      continue
    if l == 1:
      yield [j]
    else:
      for k in getlist(n, j+1, l-1):
        yield [j] + k

def get_ieq(l, i, j):
  # l[i] >= l[j]
  ieq = [0] * len(l)
  ieq[i] = -1
  ieq[j] = 1
  ieq[0] = -l[i]+l[j]
  return ieq

def get_tropical_line_complement(ll):
  l = [ll[-1]] + ll[:-1]
  I = filter(lambda j: l[j] != Infinity, range(len(l)))
  p = []
  for ii in range(len(I)):
    i = I[ii]
    for jj in range(ii + 1, len(I)):
      j = I[jj]
      eqs1 = [get_ieq(l, i, j)]
      eqs2 = [get_ieq(l, j, i)]
      for kk in range(len(I)):
        k = I[kk]
        if kk == ii or kk == jj:
          continue
        eqs1.append(get_ieq(l, i, k))
        eqs2.append(get_ieq(l, j, k))
      pp = Polyhedron(ieqs=eqs1)
      if p.count(pp) == 0:
        p.append(pp)
      pp = Polyhedron(ieqs=eqs2)
      if p.count(pp) == 0:
        p.append(pp)
  return p 

def get_verts(bound, n):
  if n == 0:
    yield []
  else:
    for v in get_verts(bound, n-1):
      yield [bound] + v
      yield [-bound] + v

def amoeba_from_B(B, prime, bound=10):
  verts = map(lambda v: v, get_verts(bound, len(B[0])-1))
  amoeba = []
  Q = Qp(prime)
  val = lambda l: map(lambda q: Q(q).valuation(), l)
  Bm = matrix(B)
  for II in getlist(len(B), 0, len(B[0])-1, [0]):
    print II
    p = [Polyhedron(verts)]
    I = II + [0]
    D = map(lambda i: B[i], I)
    Dm = matrix(D)
    print Dm
    if Dm.determinant() == 0:#Dm.is_invertible() == False:
      continue
    BD = Bm*Dm.inverse()
    BL = map(lambda l: list(l), list(BD))
    QL = map(lambda l: val(l), BL)
    mn = lambda l, r: min(map(lambda i: l[i]+r[i], range(len(l))))
    FQ = lambda r: map(lambda l: mn(l,r), QL)
    FI = lambda x: (matrix(FQ(x + [0]))*Bm).list()

    CC = map(lambda l: get_tropical_line_complement(l), QL)
    print len(CC)
    CC = filter(lambda c: c != [], CC)
    print len(p)
    for C in CC:
      p = map(lambda c: map(lambda pp: c.intersection(pp), p), C)
      p = sum(p, [])
      p = filter(lambda pp: pp.dim() == C[0].dim(), p)
    print len(p)
    for pp in p:
      v = pp.vertices()
      v = map(lambda vv: FI(vv), v)
      amoeba.append(Polyhedron(v))
  return amoeba

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
  return amoeba_from_B(B, p, 100)
