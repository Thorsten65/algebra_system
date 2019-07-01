from AbstrakteRinge import *

# Es werden die Klassen Polynomring und PolynomringElement definiert.



class Polynomring(Ring):

    """Instanziierbare Klasse"""
    
    def __init__(self,R : Ring,variablenname = "x"):

        if not isinstance(R,Ring):
            raise TypeError()   
             
        self.basisring = R
        
        self.null = PolynomringElement(self.basisring.null,self)
        self.eins = PolynomringElement(self.basisring.eins,self)

        
        if not type(variablenname) == str:
            raise TypeError()
        
        self.variablenname = variablenname

        self.variable = PolynomringElement(RingTupel([self.basisring.null, \
                            self.basisring.eins]), self)

        self._frier()
        
        
    def __str__(self):

            return "{0}[{1}]".format(self.basisring.__str__(),self.variablenname)

    
    def __eq__(self,other):

        if not super().__eq__(other):
            return False
        return (self.basisring == other.basisring)



    def element(self,*info):
        return PolynomringElement(*(info+(self,)))

 
    def monom (self,exp):

        if not type(exp) == int:
            raise TypeError
        if exp < 0:
            raise RuntimeError("Der Exponent ist negativ.")

        # Achtung:
        # Man muss die Eins und die Null aus dem Basisring nehmen.
        tup = exp*[self.basisring.null] + [self.basisring.eins]
        return PolynomringElement(RingTupel(tup),self)


    
class PolynomringElement(RingElement):

    """Instanziierbare Klassen"""
    
    def __init__(self,koeffizienten,polyring : Ring):

        if not isinstance(polyring,Polynomring):
            raise TypeError("Das zweite Objekt ist kein Polynomring.")

        # Die folgende Definition folgt der allgemeinen Spezifikation für Ring
        self.ring = polyring
        self.basisring = polyring.basisring

        if isinstance(koeffizienten,RingElement):
        # Es gibt hier zwei Fälle:
        # Das Element ist aus dem Basisring.
        # Dann wird es eine Konstante. Die Koeffizienten des neuen
        # Elementes sind gegeben durch das 1-Tupel koeffizienten
        # Oder:
        # Das Element ist aus einem zu self.ring vergleichbaren Ring.
        # Dies ist dann natürlich auch ein Polynomring.
        # Dann muss es von einem Polynomring in einen isomorphen Polynomring
        # überführt werden.
        # Das heißt genau, dass die koeffizienten überführt werden müssen.
        # Diese sind dann koeffizienten.koeffizienten.
        
            if koeffizienten.ring == self.basisring:
                self.koeffizienten = RingTupel([koeffizienten],self.basisring)
            elif koeffizienten.ring == self.ring:
                self.koeffizienten = RingTupel(koeffizienten.koeffizienten,self.basisring)
            else:
                raise RuntimeError ("Das Ringelement kann nicht in den Polynomring abgebildet werden.")
                
        elif isinstance(koeffizienten,int):
            self.koeffizienten = RingTupel([koeffizienten*self.basisring.eins])

        elif isinstance(koeffizienten,RingTupel) or \
                    type(koeffizienten) == tuple or \
                    type(koeffizienten) == list:
            self.koeffizienten = RingTupel(koeffizienten,self.basisring)
            
        else:
            raise TypeError("Das erste Objekt muss vom Typ int, RingElement, RingTupel, tuple oder list sein.")


        # Die Tupel, die Polynome repräsentieren, sollen mit einem
        # Nicht-Nulleintrag enden.
        # Damit ist dann der Grad des Polynoms gleich der Länge
        # des Tupels.
            
        self.koeffizienten.auslaufende_nullen_loeschen()
        self.grad = self.koeffizienten.laenge-1

        self._frier()
            
            
    def drucke_element(self):

        ausgabe = ""       

        for i in range(0,self.grad+1):

            if not self.koeffizienten[i] == self.basisring.null:

                if not ausgabe == "":
                    ausgabe = ausgabe + " + " 


                if self.basisring == Z:
                    ausgabe = ausgabe + str(self.koeffizienten[i])
                else:
                    ausgabe = ausgabe + self.koeffizienten[i].drucke_element_mit_klammern()

                if i >= 1:
                    ausgabe = ausgabe + "*" + self.ring.variablenname

                if i >= 2:
                    ausgabe = ausgabe + "^" + str(i)
                    
        if ausgabe == "":
            if self.basisring == Z:
                ausgabe = "0"

            else:
                ausgabe = self.basisring.null.drucke_element()

        return ausgabe


    def drucke_element_mit_klammern(self):
    
        rueck = self.drucke_element()
    
        i = 0
        while self.koeffizienten[i] == self.basisring.null:
            i += 1
            
        if i < self.grad:
            rueck = "(" + rueck + ")"
    
        return rueck
    


    def __eq__(self,other):
    
        if super().__eq__(other) == False:
            return False
        
        return self.koeffizienten == other.koeffizienten
        
    
    def koeffizient(self,i):

        if not type(i) == int:
            raise TypeError("Eingabe ist keine Zahl.")

        if i >= 0 and i <= self.grad:
            return self.koeffizienten[i]
        else:
            return self.basisring.null
 
    
            
    def __neg__(self):

            return PolynomringElement(-self.koeffizienten,self.ring)


    def __radd__(self,other):

        super().__radd__(other)

        if type(other) == int:
            other = other*self.ring.eins

        g = max(self.grad,other.grad)
        
        koeff1 = RingTupel.zusammenfuegen(self.koeffizienten,
                            RingTupel((g - self.grad)*[self.basisring.null]))
        koeff2 = RingTupel.zusammenfuegen(other.koeffizienten,
                            RingTupel((g - other.grad)*[other.basisring.null]))
         

        return PolynomringElement(koeff1 + koeff2,self.ring)

        
    def __rmul__(self,other):

        if not RingElement.test(other):
            raise TypeError("Der erste Faktor ist keine Zahl und kein Ringelement.")

        if not (type(other) == int or RingElement.ring(other) == self.basisring or RingElement.ring(other) == self.ring):
            raise TypeError("Die Elemente liegen in inkompatiblen Ringen.")
        
        if type(other) == int or RingElement.ring(other) == self.basisring:
            return(PolynomringElement(other*self.koeffizienten,self.ring))
        
        if self == self.ring.null or other == self.ring.null:
            return self.ring.null

        # Jetzt geht's los:
            
        tup = RingTupel((self.grad + other.grad+1)*[self.basisring.null])

        for d in range(0,self.grad+other.grad+1):
            for i in range(0,d+1):
                tup.koeffizienten[d] = tup.koeffizienten[d] + \
                           self.koeffizient(i) * other.koeffizient(d-i)
       
        return PolynomringElement(tup,self.ring)


    def invers(self):
            
        if not self.grad == 0:
            raise RuntimeError(
                "Der Grad des zu invertierenden Elementes muss null sein.")

        # Das zu invertierende Ringelement:
        element = self.koeffizienten[0]
            
        if self.basisring == Z:
            # Über Z kann man nur +1 und -1 invertieren, das Ergebnis ist \
            # dann das Element selbst.
            if not (element == 1 or element == -1):
                raise RuntimeError("Das Element ist nicht invertierbar.")
            return self
        else:                
            return PolynomringElement(element.invers(),self.ring)

        
    
class PolynomringElement(RingElement):

    """Instanziierbare Klassen"""
    
    def __init__(self,koeffizienten,polyring : Ring):

        if not isinstance(polyring,Polynomring):
            raise TypeError("Das zweite Objekt ist kein Polynomring.")

        # Die folgende Definition folgt der allgemeinen Spezifikation für Ring
        self.ring = polyring
        self.basisring = polyring.basisring

        if isinstance(koeffizienten,RingElement):
        # Es gibt hier zwei Fälle:
        # Das Element ist aus dem Basisring.
        # Dann wird es eine Konstante. Die Koeffizienten des neuen
        # Elementes sind gegeben durch das 1-Tupel koeffizienten
        # Oder:
        # Das Element ist aus einem zu self.ring vergleichbaren Ring.
        # Dies ist dann natürlich auch ein Polynomring.
        # Dann muss es von einem Polynomring in einen isomorphen Polynomring
        # überführt werden.
        # Das heißt genau, dass die koeffizienten überführt werden müssen.
        # Diese sind dann koeffizienten.koeffizienten.
        
            if koeffizienten.ring == self.basisring:
                self.koeffizienten = RingTupel([koeffizienten],self.basisring)
            elif koeffizienten.ring == self.ring:
                self.koeffizienten = RingTupel(koeffizienten.koeffizienten,self.basisring)
            else:
                raise RuntimeError ("Das Ringelement kann nicht in den Polynomring abgebildet werden.")
                
        elif isinstance(koeffizienten,int):
            self.koeffizienten = RingTupel([koeffizienten*self.basisring.eins])

        elif isinstance(koeffizienten,RingTupel) or \
                    type(koeffizienten) == tuple or \
                    type(koeffizienten) == list:
            self.koeffizienten = RingTupel(koeffizienten,self.basisring)
            
        else:
            raise TypeError("Das erste Objekt muss vom Typ int, RingElement, RingTupel, tuple oder list sein.")


        # Die Tupel, die Polynome repräsentieren, sollen mit einem
        # Nicht-Nulleintrag enden.
        # Damit ist dann der Grad des Polynoms gleich der Länge
        # des Tupels.
            
        self.koeffizienten.auslaufende_nullen_loeschen()
        self.grad = self.koeffizienten.laenge-1

        self._frier()
            
            
    def drucke_element(self):

        ausgabe = ""       

        for i in range(0,self.grad+1):

            if not self.koeffizienten[i] == self.basisring.null:

                if not ausgabe == "":
                    ausgabe = ausgabe + " + " 


                if self.basisring == Z:
                    ausgabe = ausgabe + str(self.koeffizienten[i])
                else:
                    ausgabe = ausgabe + self.koeffizienten[i].drucke_element_mit_klammern()

                if i >= 1:
                    ausgabe = ausgabe + "*" + self.ring.variablenname

                if i >= 2:
                    ausgabe = ausgabe + "^" + str(i)
                    
        if ausgabe == "":
            if self.basisring == Z:
                ausgabe = "0"

            else:
                ausgabe = self.basisring.null.drucke_element()

        return ausgabe


    def drucke_element_mit_klammern(self):
    
        rueck = self.drucke_element()
    
        i = 0
        while self.koeffizienten[i] == self.basisring.null:
            i += 1
            
        if i < self.grad:
            rueck = "(" + rueck + ")"
    
        return rueck
    


    def __eq__(self,other):
    
        if super().__eq__(other) == False:
            return False
        
        return self.koeffizienten == other.koeffizienten
        
    
    def koeffizient(self,i):

        if not type(i) == int:
            raise TypeError("Eingabe ist keine Zahl.")

        if i >= 0 and i <= self.grad:
            return self.koeffizienten[i]
        else:
            return self.basisring.null
 
    
            
    def __neg__(self):

            return PolynomringElement(-self.koeffizienten,self.ring)


    def __radd__(self,other):

        super().__radd__(other)

        if type(other) == int:
            other = other*self.ring.eins

        g = max(self.grad,other.grad)
        
        koeff1 = RingTupel.zusammenfuegen(self.koeffizienten,
                            RingTupel((g - self.grad)*[self.basisring.null]))
        koeff2 = RingTupel.zusammenfuegen(other.koeffizienten,
                            RingTupel((g - other.grad)*[other.basisring.null]))
         

        return PolynomringElement(koeff1 + koeff2,self.ring)

        
    def __rmul__(self,other):

        if not RingElement.test(other):
            raise TypeError("Der erste Faktor ist keine Zahl und kein Ringelement.")

        if not (type(other) == int or RingElement.ring(other) == self.basisring or RingElement.ring(other) == self.ring):
            raise TypeError("Die Elemente liegen in inkompatiblen Ringen.")
        
        if type(other) == int or RingElement.ring(other) == self.basisring:
            return(PolynomringElement(other*self.koeffizienten,self.ring))
        
        if self == self.ring.null or other == self.ring.null:
            return self.ring.null

        # Jetzt geht's los:
            
        tup = RingTupel((self.grad + other.grad+1)*[self.basisring.null])
        for d in range(0,self.grad+other.grad+1):
            for i in range(0,d+1):
                tup.koeffizienten[d] = tup.koeffizienten[d] + \
                           self.koeffizient(i) * other.koeffizient(d-i)
       
        return PolynomringElement(tup,self.ring)


    def invers(self):
            
        if not self.grad == 0:
            raise RuntimeError(
                "Der Grad des zu invertierenden Elementes muss null sein.")

        # Das zu invertierende Ringelement:
        element = self.koeffizienten[0]
            
        if self.basisring == Z:
            # Über Z kann man nur +1 und -1 invertieren, das Ergebnis ist \
            # dann das Element selbst.
            if not (element == 1 or element == -1):
                raise RuntimeError("Das Element ist nicht invertierbar.")
            return self
        else:                
            return PolynomringElement(element.invers(),self.ring)	    
