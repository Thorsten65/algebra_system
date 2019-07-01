import random
from Tocas import Ring, RingElement, Polynomring, PolynomringElement, RingTupel, GanzzahlRestklassenring
from Extension import polynomring_extension
from Extension import Primzahl

class PolynomRestklassenring(Ring):
    def __init__(self, f: PolynomringElement):
        if not isinstance(f, PolynomringElement):
            raise RuntimeError("Das angegebene Objekt ist kein Polynom.")
        
        if f.grad > 0 and f.koeffizienten.koeffizienten[f.grad] != 1:
            for d in range(f.grad + 1):
                f.koeffizienten.koeffizienten[d] /= f.koeffizienten.koeffizienten[f.grad]

        self.modulus = f
        self.null = PolynomRestklassenringElement(f.ring.null, self)
        self.eins = PolynomRestklassenringElement(f.ring.eins, self)
        self.erzeuger = PolynomRestklassenringElement(f.ring.variable, self)
        self.elements = []

        if self.modulus.ring.basisring.ist_endlicher_koerper and  \
           not Primzahl.rabin_test_reducible(f):
            self.ist_endlicher_koerper = True
        self._frier()

        
    def __str__(self):
        return "{}/({})".format(self.modulus.ring, self.modulus.drucke_element())

    def __eq__(self, other):
        if not super().__eq__(other):
            return False
        return (self.modulus == other.modulus)

    def element(self, a):
        return PolynomRestklassenringElement(a, self)

    def random(self):
        if isinstance(self.modulus.basisring, GanzzahlRestklassenring):
            r = RingTupel((self.modulus.grad + 1) *
                          [self.modulus.basisring.null])

            for d in range(self.modulus.grad + 1):
                r.koeffizienten[d] = random.randint(
                    0, self.modulus.basisring.modulus - 1)

            return PolynomRestklassenringElement(r, self)
        elif isinstance(self.modulus.basisring, PolynomRestklassenring):
            r = RingTupel((self.modulus.grad + 1) *
                          [self.modulus.basisring.null])

            for d in range(self.modulus.grad + 1):
                r.koeffizienten[d] = self.modulus.basisring.random()

            return PolynomRestklassenringElement(r, self)
        else:
            raise TypeError(
                "random ist nur für PolynomRestklassenringe mit Restklassenringen als Basisring implementiert")


class PolynomRestklassenringElement(RingElement):
    def __init__(self, p, r: PolynomRestklassenring):
        self.ring = r

        if isinstance(p, PolynomRestklassenringElement):
            if p.ring.modulus.basisring != r.modulus.basisring:
                raise RuntimeError(
                    "Polynom und Ring nicht vom selben Basisring.")

            if p.ring.modulus % self.ring.modulus == self.ring.modulus.ring.null:
                self.wert = p.wert % self.ring.modulus

            else:
                raise RuntimeError("Die Moduli passen nicht zusammen.")

        else:
            if not isinstance(p, PolynomringElement):
                p = PolynomringElement(p, r.modulus.ring)

        self.wert = p % self.ring.modulus
        self._frier()

    def __str__(self):
        return '{} in {}'.format(self.drucke_element(), self.ring)

    def drucke_element(self):
        return '{}'.format(self.wert.drucke_element())

    def __eq__(self, other):
        if not super().__eq__(other):
            return False
        return self.wert == other.wert

    def __neg__(self):
        return PolynomRestklassenringElement(-self.wert, self.ring)

    def __add__(self, other):
        if isinstance(other, RingElement) and other.ring == self.ring.modulus.basisring:
            return self + other * self.ring.eins

        super().__add__(other)

        return PolynomRestklassenringElement(self.wert+other.wert, self.ring)

    def __rmul__(self, other):
        super().__rmul__(other)

        if type(other) == int:
            return PolynomRestklassenringElement(RingElement.intmult(other, self.wert), self.ring)

        if not isinstance(other, RingElement):
            raise TypeError(
                "Der erste Faktor ist keine Zahl und kein Ringelement.")

        if other.ring == self.ring.modulus.basisring:
            return PolynomRestklassenringElement(other*self.wert.koeffizienten, self.ring)

        if self == self.ring.null or other == self.ring.null:
            return self.ring.null

        r = PolynomRestklassenringElement(
            (self.wert.grad + 1)*[self.wert.basisring.null], self.ring)

        for d in range(other.wert.grad, 0, -1):
            r += other.wert.koeffizient(d) * self
            r = PolynomRestklassenringElement(RingTupel(
                [r.wert.basisring.null] + r.wert.koeffizienten.koeffizienten), self.ring)

        r += other.wert.koeffizient(0) * self

        return r

    def invers(self):
        _, u, _ = Polynomring.ExtGGT(self.wert, self.ring.modulus)
        return PolynomRestklassenringElement(u, self.ring)

    def endlicher_koerper():
        pass
