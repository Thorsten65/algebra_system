from AbstrakteRinge import *
from Polynomringe import *
from random import randint
from GanzzahlRestklassenringe import *
from BruchzahlRing import *
# from itertools import product


class PolynomRestklassenring(Ring):

    def __init__(self, ring, polynom):
        self.basisring = ring.basisring
        self.ring = ring
        self.modulus = polynom
        self.eins = PolynomRestklassenringElement(polynom.ring.eins, self)
        self.null = PolynomRestklassenringElement(polynom.ring.null, self)
        self._frier()

    def get_rest(self, polynom):
        return self.get_rest_from_polynomring(polynom)

    # def get_rest_from_polynomring(self, polynom):
    #    rest = []
    #    degree = polynom.grad
    #    for rest_polynom in product(range(2), repeat=degree):
    #        rest.append(RingTupel(list(rest_polynom)))
    #    return rest

    def element(self, element):
        return PolynomRestklassenringElement(element, self)

    def __str__(self):
        return "{}/({})".format(self.modulus.ring, self.modulus.drucke_element())

    def random(self):
        if isinstance(self.modulus.basisring, GanzzahlRestklassenring):
            r = RingTupel((self.modulus.grad+1)*[self.modulus.basisring.null])
            for d in range(self.modulus.grad+1):
                r.koeffizienten[d] = randint(0, self.modulus.basisring.modulus-1)
            return PolynomRestklassenringElement(r, self)

        elif isinstance(self.modulus.basisring, PolynomRestklassenring):
            r = RingTupel((self.modulus.grad + 1)*[self.modulus.basisring.null])
            for d in range(self.modulus.grad + 1):
                r.koeffizienten[d] = self.modulus.basisring.random()
            return PolynomRestklassenringElement(r, self)

        else:
            raise TypeError("random ist nur für PolynomRestklassenringe mit "
                            "Restklassenringen als Basisring implementiert")


class PolynomRestklassenringElement(RingElement):
    def __init__(self, ringelement, polynomrestklassenring):
        if not isinstance(polynomrestklassenring, PolynomRestklassenring):
            raise TypeError('Ring muss ein PolynomRestklassenring sein.')

        self.ring = polynomrestklassenring
        if ringelement == polynomrestklassenring.ring.eins:
            self.element = polynomrestklassenring.ring.eins

        elif ringelement == polynomrestklassenring.ring.null:
            self.element = polynomrestklassenring.ring.null

        else:
            self.element = self.modulos(ringelement, self.ring.modulus)

    def modulos(self, ringelement, modulus):
        if isinstance(ringelement, PolynomRestklassenringElement):
            if ringelement.ring.modulus.basisring != self.ring.modulus.basisring:
                raise RuntimeError("Ring des Elementes stimmen nicht ueberein.")

            q, r = Polynomring.polynom_division(self.ring.modulus,
                                                ringelement.ring.modulus)
            if not r == self.ring.modulus.ring.null:
                raise RuntimeError("Die Moduli stimmen nicht ueberein.")

            ringelement = ringelement.element

        elif isinstance(ringelement, RingElement):

            if isinstance(ringelement, PolynomringElement):
                ringelement = ringelement

            elif not (ringelement.ring == self.ring.basisring):
                raise RuntimeError("Ring des Elementes stimmen nicht ueberein")

            else:
                ringelement = PolynomringElement([ringelement], self.ring.ring)

        elif isinstance(ringelement, tuple):
            ring_tupel = RingTupel(list(ringelement))
            ringelement = PolynomringElement(ring_tupel, self.ring.modulus.ring)

        elif isinstance(ringelement, RingTupel):
            ringelement= PolynomringElement(ringelement, self.ring.modulus.ring)

        if isinstance(ringelement, int):
            ringelement = PolynomringElement(ringelement, self.ring.modulus.ring)


        element = Polynomring.polynom_division(ringelement, modulus)[1]
        return element

    def __radd__(self, other):
        if isinstance(other, int):
            other = other*self.ring.eins

        if isinstance(other, RingElement) and other.ring == self.ring.modulus.basisring:
            return self+other*self.ring.eins

        return PolynomRestklassenringElement(self.element+other.element, self.ring)

    def __rmul__(self, other):

        if isinstance(other, int):
            ring_element = RingElement.intmult(other, self.element)
            return PolynomRestklassenringElement(ring_element, self.ring)

        if other.ring == self.ring.modulus.basisring:
            return PolynomRestklassenringElement(self.element*other.element, self.ring)

        if self == self.ring.null or other == self.ring.null:
            return self.ring.null

        return PolynomRestklassenringElement(self.element*other.element, self.ring)

    def __eq__(self, other):
        if not super().__eq__(other):
            return False
        return self.element == other.element

    def __neg__(self):
        return PolynomRestklassenringElement(-self.element, self.ring)

    def __str__(self):
        return "{} in {}".format(self.drucke_element(), self.ring)

    def drucke_element(self):
        return "{}".format(self.element.drucke_element())

    def invers(self):
        q, r = Polynomring.ext_ggt(self.element, self.ring.modulus)
        return PolynomRestklassenringElement(q, self.ring)

