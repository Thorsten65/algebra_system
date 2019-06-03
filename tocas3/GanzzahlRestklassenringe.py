from AbstrakteRinge import *

# Zufall:

from random import randint


class GanzzahlRestklassenring(Ring):

    """Instanziierbare Klasse"""


    def __init__(self,modulus : int):
    
        if not isinstance(modulus,int):
            raise RuntimeError("Das angegebene Objekt ist keine ganze Zahl.")
        if modulus <= 0:
            raise RuntimeError("Der Modulus ist nicht positiv.")
        self.modulus = modulus
        self.null = GanzzahlRestklassenringElement(0,self)
        self.eins = GanzzahlRestklassenringElement(1,self)

        self._frier()

    def __str__(self):
        
        return "Z/{0}Z".format(self.modulus)
    
        
    def __eq__(self,other):
        
        if not super().__eq__(other):
            return False
        return (self.modulus == other.modulus)


    def element(self,a):
        return GanzzahlRestklassenringElement(a,self)


    def zufaellig(self):
        
        return GanzzahlRestklassenringElement(randint(0,self.modulus-1),self.modulus)

    


class GanzzahlRestklassenringElement(RingElement):
 
    """Instanziierbare Klasse"""
     
    def __init__(self,a,n):
        
        if isinstance(n,int):
            if n <= 0:
                raise RuntimeError("Der Modulus ist nicht positiv.")
            self.ring =GanzzahlRestklassenring(n)

        elif isinstance(n,GanzzahlRestklassenring):
            self.ring = n

        else:
            raise RuntimeError("Das zweite angegebene Objekt ist keine Zahl und kein Ganzzahl-Restklassenring.")

        if type(a) == int:

            self.wert = a % self.ring.modulus

        elif type(a) == GanzzahlRestklassenringElement:

            if a.ring.modulus % self.ring.modulus == 0:
                self.wert = a.wert % self.ring.modulus

            else:
                raise RuntimeError("Die Moduli passen nicht zusammen.")

        else:
            raise RuntimeError("Das erste angegebene Objekt ist keine Zahl und kein Ganzzahl-Restklassenringelement.")

        self._frier()


    def drucke_element(self):

        return "[{0}]".format(self.wert)



    def __eq__(self,other):

        if not super().__eq__(other):
            return False
        return self.wert == other.wert

    

    # Jetzt kommt das Überladen / Definieren der arithmetischen Operatoren

    
    def __neg__(self):

        return GanzzahlRestklassenringElement(-self.wert, self.ring)


    def __radd__(self,other):
        
        super().__radd__(other)

        if type(other) == int:
            return GanzzahlRestklassenringElement(self.wert+other,self.ring)

        return GanzzahlRestklassenringElement(self.wert+other.wert,self.ring)

        

    def __rmul__(self,other):
 
        super().__rmul__(other)

        # Der eine Faktor ist a:
        # Dieser Faktor ist ein Ringelement:
        
        a = self
                
        # Der andere Faktor ist b
        # Die Eingabe ist ja: other*self
        # Die Bedingung lautet: other muss entweder eine Zahl sein
        # oder ein Ganzzahl-Restklassenringelement, dessen modulus
        # ein Vielfaches von dem von a ist:
        # Dieser Faktor wird (zunächst) als Zahl dargestellt. 

        if type(other) == int:
            b = other
        elif (isinstance(other,GanzzahlRestklassenringElement) and 
                            other.ring.modulus % self.ring.modulus == 0):
            b = other.wert % self.ring.modulus
        else:
            raise RuntimeError("Die Elemente können nicht multipliziert werden.")
 
        # Umwandlung in 2-adische Schreibweise:
        
        zweiadisch = zwei_adisch(b)

        # Jetzt kommt double-and-add


        c = self.ring.null

        for i in range (0,len(zweiadisch)-1):
            if zweiadisch[i] == '1':
                c = c + a
            c = c + c

        # Am Ende muss man noch einmal addieren, ohne zu multiplizieren.
     
        if (len(zweiadisch) > 0 
                and zweiadisch[len(zweiadisch)-1] == '1'):
            c = c + a

        return c



    def invers(self):

        a = self.wert
        b = self.ring.modulus

        u, s = GanzzahlRestklassenringElement(1,self.ring), GanzzahlRestklassenringElement(0,self.ring)

        # Sehr kompakte GGT-Modulo-Berechnung
        while b!=0:
            q=a//b
            a, b = b, a-q*b
            u, s = s, u-q*s

        if a != 1:
            raise InvertierungsFehler(self)

        return u

