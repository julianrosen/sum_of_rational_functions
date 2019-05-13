from sage.matrix.special import vandermonde
from IPython.display import display, Math, Latex

from sage.matrix.special import vandermonde
from IPython.display import display, Math, Latex

def elementary(i,n,p):
    """
    Expression for the i-th elementary symmetric function
    in 1/1,1/2,...,1/(p-1) modulo p^n
    """
    if i == 0:
        return 1
    if n <= 1:
        return 0
    if n == 2 and i%2 == 1:
        return 0
    number_of_terms = i+n-(1 if (i+n)%2==0 else 2)
    binomial_list = [binomial(j*p,p)/j for j in range(1,number_of_terms+1)]
    M = vandermonde(range(number_of_terms),QQ)
    Mi = M.inverse()
    return expand((Mi*vector(binomial_list))[i]/p^i)

def power_sum(i,n,p):
    """
    Expression for the i-th power sum symmetric function
    in 1/1,1/2,...,1/(p-1) modulo p^n
    Comes from Newton's formula
    """
    if i == 1:
        return elementary(1,n,p)
    T = sum((-1)^(j-1)*elementary(j,n,p)*power_sum(i-j,n,p) for j in range(1,i))
    return T + (-1)^(i-1)*i*elementary(i,n,p)

def sum_0_to_p(f,n,p):
    # Returns a polynomial in p, p^{-1}, and binomial
    # coefficients {kP choose P}, whose value is 
    # congruent to \sum_{i=0}^{P-1} f(i) modulo P^n
    K = parent(f)
    x = K.gens()[0]
    pf = f.partial_fraction_decomposition()
    var('m')
    F = sum(pf[0].subs({x:m}),m,0,p-1) # Polynomial part
    for s in pf[1]: # Each s corresponds to a term u / (x+a)^i
        L = s.factor()
        u = L.unit() # Constant factor associated with term
        L = list(L)[0]
        a = L[0]-x
        if a not in NN:
            print a
            print "Poles need to be negative integers"
            return None
        i = -L[1] # Exponent in denominator
        F += u*power_sum(i,n,p)
        F -= u*sum(1/j^i for j in range(1,a)) # Correction for "missed" initial terms
        F += u/p^i # Term x+a = p
        for j in range(1,a): # Terms with x+a > p
            F += u*sum(binomial(-i,q)*p^q/j^(q+i) for q in range(n))
    return F

def nice_latex(F):
    s = str(latex(F))
    s=s.replace('\\, p \\choose','p \\choose')
    return s

def test(f,n,P):
    """
    Test the computation
    """
    var('p')
    K = parent(f)
    x = K.gens()[0]
    direct_summation = sum(f.subs({x:j}) for j in range(P))
    F = sum_0_to_p(f,n,p)
    binomial_evaluation = F.subs({p:P})
    r = direct_summation - binomial_evaluation
    return Integer(r.numerator())%P^n == 0

def pretty_congruence(f,m):
    """
    Displays the congruence in a pretty fashion
    """
    var('p,n')
    K = parent(f)
    x = K.gens()[0]
    F = sum_0_to_p(f,m,p)
    S = r'\sum_{n=1}^{p-1}' + latex(f.subs({x:n})) + r'\equiv' + nice_latex(F) + '\\mod p^{%i}'%m
    display(Math(S))
    return None
    
    
# Example:
# R.<x> = PolynomialRing(QQ)
# K = R.fraction_field()
# f = x^2/((x+3)*(x+2))
# pretty_congruence(f,3)