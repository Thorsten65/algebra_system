from AbstrakteRinge import *


class Bruchzahlring(Ring):
    """Ring der rationalen Zahlen"""

    def __init__(self):
        self.null = BruchzahlringElement(0, 1, self)
        self.eins = BruchzahlringElement(1, 1, self)

        self._frier()

    def __str__(self):
        return "Q"

    def __eq__(self, other):
        return super().__eq__(other)

    def element(self, p, q=1):
        if isinstance(p, BruchzahlringElement):
            return p
        if not (isinstance(p, int) and isinstance(q, int)):
            raise TypeError("Elemente sind keine ganzen Zahlen.")

        return BruchzahlringElement(p, q, self)


class BruchzahlringElement(RingElement):
    """Rationale Zahl"""

    def __init__(self, a, b, ring=None):
        if not (isinstance(a, int) and isinstance(b, int)):
            raise TypeError(
                "Die angegebenen Parameter sind keine ganzen Zahlen.")

        if b < 0:
            a = -a
            b = -b

        ggt = Ganzzahlring.ext_ggT(a, b)[0]
        self.a = a // ggt
        self.b = b // ggt
        if not ring:
            ring = Bruchzahlring()
        self.ring = ring

        self._frier()

    def __str__(self):
        return '{} in {}'.format(self.drucke_element(), self.ring)

    def drucke_element(self):
        return '{}/{}'.format(self.a, self.b)

    def drucke_element_mit_klammern(self):
        return '({}/{})'.format(self.a, self.b)

    def __eq__(self, other):
        if not super().__eq__(other):
            return False
        return self.a == other.a and self.b == other.b

    def __neg__(self):
        return BruchzahlringElement(-self.a, self.b, self.ring)

    def __add__(self, other):
        if not isinstance(other, int):
            super().__add__(other)

        if type(other) == int:
            return BruchzahlringElement(self.a + other * self.b, self.b, self.ring)
        elif isinstance(other, BruchzahlringElement):
            return BruchzahlringElement(self.a * other.b + other.a * self.b, self.b * other.b, self.ring)

    def __rmul__(self, other):
        super().__rmul__(other)

        if type(other) == int:
            return BruchzahlringElement(self.a * other, self.b, self.ring)
        elif isinstance(other, BruchzahlringElement):
            return BruchzahlringElement(self.a * other.a, self.b * other.b, self.ring)

    def invers(self):
        if self.a > 0:
            return BruchzahlringElement(self.b, self.a, self.ring)
        else:
            return BruchzahlringElement(-self.b, -self.a, self.ring)


