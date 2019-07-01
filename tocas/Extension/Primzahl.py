import random
import math

from Tocas import Ganzzahlring, Polynomring, GanzzahlRestklassenring

import Extension.polynomring_extension
from Extension.PolynomRestklassenring import PolynomRestklassenringElement, PolynomRestklassenring

def miller_rabin(n, k=40):
    if not isinstance(n, int):
        raise TypeError('Argumente nicht vom Typ int')

    if n <= 1:
        return False

    if n == 2:
        return True
    if n % 2 == 0:
        return False

    s = 0
    t = n - 1
    while t % 2 == 0:
        s += 1
        t //= 2

    ZnZ = GanzzahlRestklassenring(n)

    for _ in range(k):
        a = ZnZ.element(random.randrange(1, n - 1))
        a = a ** t

        if a == ZnZ.eins:
            continue
        for _ in range(s):
            if a == -ZnZ.eins:
                break
            elif a == ZnZ.eins:
                a = a ** 2
        else:
            return False

    return True

def _ord(r, n):
    if Ganzzahlring.ext_ggT(r, n)[0] > 1:
        return math.inf


    ZrZ = GanzzahlRestklassenring(r)

    k = 1
    while ZrZ.element(n) ** k != ZrZ.eins:
        k += 1

    return k


def _phi(r):
    phi = 0

    for i in range(1, r):
        if Ganzzahlring.ext_ggT(r, i)[0] == 1:
            phi += 1

    return phi


def aks(n):
    from Extension.PolynomRestklassenring import PolynomRestklassenring

    if not isinstance(n, int):
        raise TypeError('Argumente nicht vom Typ int')

    if n <= 1:
        return False

    for a in range(2, math.ceil(math.sqrt(n))):
        b = 2
        while a ** b <= n:
            if n == a ** b:
                return False
            b += 1

    r = 1
    while _ord(r, n) <= math.log2(n) ** 2:
        r += 1

    for a in range(2, r + 1):
        ggt = Ganzzahlring.ext_ggT(a, n)[0]

        if ggt > 1 and ggt < n:
            return False

    if n <= r:
        return True

    K_X_n = Polynomring(GanzzahlRestklassenring(n))
    K_X_n_f = PolynomRestklassenring(K_X_n.variable ** r - K_X_n.eins)

    for a in range(1, math.floor(math.sqrt(_phi(r)) * math.log2(n) + 1)):
        if K_X_n_f.element([a, 1]) ** n != K_X_n_f.erzeuger ** n + K_X_n_f.element([a]):
            return False

    return True


def primes(n):
    divisors = [d for d in range(2, n//2+1) if n % d == 0]
    values =  [d for d in divisors if \
            all(d%od != 0 for od in divisors if od != d)]
    if values:
        return values

    else:
        values.append(n)
        return values



def rabin_test_reducible(polynom):
    # monic polynomial f of degree n
    # p_1, ..., p_k all distinct prime divisors of n

    from .polynomring_extension import _polynomring_ext_ggt
    from Tocas import PolynomringElement
    from .PolynomRestklassenring import PolynomRestklassenringElement

    if not isinstance(polynom.ring.basisring, GanzzahlRestklassenring):
        raise TypeError("Nur GanzahlRestklassenring als Basisring erlaubt.") 

    try:
        q = polynom.ring.basisring.modulus
    except AttributeError:
        raise Attribute("Der Rabit-Test braucht ein Galoi-Feld als Basisring.")

    n = polynom.grad
    x = PolynomringElement([0, 1], polynom.ring)

    f = polynom
    p_x = []

    k =  primes(n)
    for j, l in enumerate(k):
        v = int(n/l)
        if v not in p_x:
            p_x.append(v)

    for i in p_x:
        print(x, q, i, f)
        h = (((x**q)**i) - x ) % f
        print(h)
        a, u, v = _polynomring_ext_ggt(f.ring, f, h)
        ggt = a*u + h*v 
        print("Erstes G: ", ggt)
        if ggt != f.ring.eins:
            return True

    g = ((((x) ** q) ** n) - x) % f
    print("Zweites G: ", g)

    return g != g.ring.null


def get_p_pow_n_elements(p, n, ring):
    from itertools import product
    rest = []
    for rest_polynom in product(range(p), repeat=n):
    #rest.append(list(rest_polynom))
        rest.append(PolynomringElement(list(rest_polynom), ring))
    return rest

def check_is_prime(n):
    from math import sqrt; from itertools import count, islice

    return n > 1 and all(n%i for i in islice(count(2), int(sqrt(n)-1)))

def check_multiple_inverse(values, modulus):
    for value in values[1:]:
        a, u, v = Ganzzahlring.ext_ggT(value, modulus)
        if ((value*u) + (v*modulus)) != 1:
            return False

    return True


def is_prime_field(p, values, modulus):
    if check_is_prime(p) and  check_multiple_inverse(values, modulus):
        return True
    return False

def is_prime_field_2(p, values, modulus, polynom):
    pass

def endlicher_koerper(p: int, n: int):

    if n <= 1:

        Z = GanzzahlRestklassenring(p)
        values = [i for i in range(p)]

        if is_prime_field(p, values, p):
            for value in values:
                Z.elements.append(GanzzahlRestklassenringElement(value, p))
            return Z

        else:
            raise TypeError("Es ist ein Primkoerper!")

    else:
        Z = GanzRestklassenring(p)
        P = PolynomRestklassenring(Z)
        modulus = PolynomringElement([1 for i in range(n+1)])
        PolyRest = PolynomRestklassenring(modulus)
        ele_list = get_p_pow_n_elements(p, n, P)
        for ele in ele_list:
            poly = PolynomringElement(ele, P)
            PolyRest.elements.append(PolynomRestklassenringElement(poly, PolyRest))
        return PolyRest

    """
    is_prime_field = check_is_prime_field()
    num_of_elements = p ** n


    if is_prime_field:
        Z = GanzzahlRestklassenring(p)
        P = Polynomring(Z)
        modulus = PolynomringElement([1 for i in range(n+1)], P)
        PolyRest = PolynomRestklassenring(modulus)
        ele_list = get_p_pow_n_elements(p, n, P)
        for ele in ele_list:
            poly = PolynomringElement(ele, P)
            PolyRest.elements.append(PolynomRestklassenringElement(poly, PolyRest))
        return PolyRest

    else:
        pass
    """

def get_rest_from_polynomring(polynom):
    from itertools import product

    if not isinstance(polynom.ring.basisring, GanzzahlRestklassenring):
        raise TypeError("Ich brauche ein Ganzzahlrestklassenring")

    rest = []
    degree = polynom.grad
    ring = polynom.ring
    modulus_degree = polynom.ring.basisring.modulus
    n_pow_n = modulus_degree ** degree
    print(modulus_degree)
    for rest_polynom in product(range(modulus_degree), repeat=degree):
        #rest.append(list(rest_polynom))
        rest.append(PolynomringElement(list(rest_polynom), ring))
    return rest



from Tocas import GanzzahlRestklassenring, RingTupel, PolynomringElement, Ganzzahlring, GanzzahlRestklassenringElement
