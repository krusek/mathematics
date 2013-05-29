This contains various files associated to my doctoral dissertation.

# sage

This is a directory containing sage code that was included in the appendices of
my doctoral dissertation.

## amoeba.sage

This file creates the reduced A-discriminant contour for n-variate,
(n+3)-nomials. Code to generate the contour for the cubic family,
`a+bx+cx^2+dx^3` is:

    load "amoeba.sage"
    A = [[0,1,2,3]]
    a = amoeba(A, 100)
    show(a, xmin=-10, xmax=10, ymin=-10, ymax=10)

For a multivariate polynomial, like `ax^6+by^3+cy+dty^6+etx^3+ftx` in `Q_p[x,y,t]` use

    A = [[6,0,0,0,3,1],[0,3,1,6,0,0],[0,0,0,1,1,1]]


## trop2d.sage

This file creates the reduced p-adic A-discriminant amoeba for
n-variate, (n+3)-nomials. Code to generate the amoeba for the cubic family,
`a+bx+cx^2+dx^3` when `p=3` is:

    load "trop2d.sage"
    A = [[0,1,2,3]]
    a = amoeba(A, 3)
    show(a, xmin=-4, xmax=25, ymin=-25, ymax=4)

For a multivariate polynomial, like `ax^6+by^3+cy+dty^6+etx^3+ftx` in `Q_p[x,y,t]` use

    A = [[6,0,0,0,3,1],[0,3,1,6,0,0],[0,0,0,1,1,1]]

## tropdd.sage

This file creates the reduced p-adic A-discriminant amoeba for
n-variate, (n+k)-nomials for general k. Code to generate the amoeba
for the quartic family (`k=4`) when `p=3` in the box `[-10,10]^3` is:

    load "tropdd.sage"
    A = [[0,1,2,3,4]]
    box = Polyhedron(list(get_verts(10,3)))
    a = amoeba(A, 3)
    show(sum(lambda b: box.intersection(b).show(), a))


