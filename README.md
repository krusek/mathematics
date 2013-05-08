This contains various files associated to my doctoral dissertation.

# sage

This is a directory containing sage code that was included in the appendices of
my doctoral dissertation.

## amoeba.sage

This file creates the reduced A-discriminant contour for n-variate,
(n+k)-nomials. Code to generate the contour for the cubic family,
`a+bx+cx^2+dx^3` is:

    A = [0,1,2,3]
    a = amoeba(A, 100)
    show(a, xmin=-10, xmax=10, ymin=-10, ymax=10)

More to come

