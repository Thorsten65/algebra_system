from AbstrakterAnfang import *
from copy import deepcopy


# Es werden abstrakte Klassen "Ring" und 
# "RingElement" definiert.
# Die Eingaben werden oftmals nur angegeben, um die Idee
# klarzumachen.
# Dies wird dann in Unterklassen überschrieben.

# Man beachte:
# Eine Klasse mit einer als "abstactmethod" (@abstractmethod) deklarierten
# Methode kann nicht instanziiert werden.
# Aber so eine Methode kann trotzdem mit super(). .. von einer
# Unterklasse aus aufgerufen werden.

# Hiernach wird dann noch die instanziierbare Klassen Ganzzahlring
# und RingTupel definiert.




class Ring(EinfrierbaresObjekt):

    """Abstrakte Klasse"""

    
    @abstractmethod
    def __init__(self,info):
        
        Ring.null, Ring.eins = None, None
      

    def __eq__(self,other):

        """Das ist ein sehr rudimentärer Gleichheitstest.
        Der reicht dann nur bei Ganzzahlring aus."""



        if not type(self) == type(other):
            return False

        return True


    def __call__(self,*info):
    
        return self.element(*info)

 
    @abstractmethod
    def element(self,*info):

        """Abstrakte Methode
        
        Für einen Ring R soll mit R.element(info)
        eine durch info definierte Instanz der 
        entsprechenden von RingElement abgeleiteten Klasse
        zurückgegeben werden.
        Siehe z.B. in Restklassenring."""

        pass

    
    def tupel(self,*koeffizienten):

        return RingTupel(*(koeffizienten+(self,)))
    


    
class RingElement(EinfrierbaresObjekt):

    """Abstrakte Klasse"""

    @abstractmethod
    def __init__(self,elementinfo,ringinfo,zusatz=None):
        
        self.ring = None


       
    @abstractmethod
    def drucke_element(self):
		
        """Abstrakte Methode
        Die Beschreibung soll immer so aussehen:

		  Ringelement  in Ring
	    
	    Für die Ausgabe von "Ringelement" ist drucke_element zuständig.
	    Dies muss noch implementiert werden.
	    Die Zusammensetzung ist dann immer gleich.
	    Deshalb ist __str__ nicht abstrakt.
	    (Und __repr__ ist schon in MeinABCObjekt implementiert.)
	    Daneben gibt es noch die Methode drucke_element_mit_klammern.
	    Dies wird von der Ausgabemethode für Polynomringelement
	    aufgerufen.
	    Es werden nur dann Klammern um das Element gesetzt,
	    wenn dies notwendig ist.
	    In diesem Sinne wird dies erstmal durch drucke_element implementiert."""

        pass


    def drucke_element_mit_klammern(self):

        """Diese Methode sollte eigentlich "drucke_element_mit_klammern_wenn_notwendig"
        heißen.
        
        Sie soll benutzt werden, wenn es für eine Multiplikation self*b
        notwendig ist, self einzuklammern. Dies ist der Fall für
        Polynome.
        
        Hier wird einfach nur auf drucke_element verwiesen."""

        return self.drucke_element()


    def __str__(self):
        
        return self.drucke_element() + "  in " + self.ring.__str__()

    

    
    @abstractmethod   
    def __eq__(self,other):

        if not type(self) == type(other):
            return False

        if not self.ring == other.ring:
            return False

        return True
        
        
    def umgebung(self):
        return self.ring

    # +x ist x:
    def __pos__(self):
        
        return self

    
    @abstractmethod   
    def __neg__(self):
        pass

    def __add__(self,other):
		
        """Addition der Form self + other

       Die Reihenfolge wird einfach umgedreht.
       Das Ziel ist: Es sollen sowohl Elemente des gleichen Rings
       als auch Ringelemente mit ganzen Zahlen addiert werden
       können."""
       
        return self.__radd__(other)


    def __radd__(self,other):

        """Addition der Form other + self
        
        Dies soll weiter implementiert werden."""
    
        if type(other) == int:
            return True
    
        if not isinstance(other,RingElement):
            raise TypeError("Der erste Summand ist kein Ringelement.")
               
        if not (type(self) == type(other)):
            raise TypeError("Die Elemente sind nicht aus vergleichbaren Ringen und nicht ganze Zahlen.")

        if not (self.ring == other.ring):
            raise RuntimeError("Die Ringe stimmen nicht überein.")

        return True


    def __sub__(self,other):

        """Subtraktion der Form self - other

        Die Reihenfolge wird umgedreht und das dann das
        Negative ausgegeben."""
		
        return -(self.__rsub__(other))
        

    def __rsub__(self,other):

        """Subtraktion der Form other - self
        
        mittels (-self)."""
        
        return other + (-self)
    

    
    def __mul__(self,other):
		
        """Multiplikation der Form self*other

        Hier werden nur die Fälle other eine ganze Zahl
        und other ein RingElement behandelt.
        
        Hierfür wird rmul aufgerufen.

        Für a*e, e ganze Zahl wird die Reihenfolge umgedreht.
        
        Für a*b, b Ringelement, wird die Reihenfolge beibehalten.
        
        Dies soll nicht weiter implementiert werden. Das lasst
        dann offen, was für a*x, x kein Ringelement, passiert."""

        
        if not (type(self) == type(other) or type(other) == int):
            return NotImplemented
            
        if type(other) == int:
            return self.__rmul__(other)

        if type(other) == type(self):
            return other.__rmul__(self)

    
    @abstractmethod
    def __rmul__(self,other):

        """Multiplikation der Form other*self"""
        
        if not (type(other) == int or isinstance(other,RingElement)):
            raise TypeError("Das erste Objekt ist keine Zahl und kein Ringelement.")

        return True


    
    @abstractmethod
    def invers(self):

        """Abstrakte Methode
        
        Könnte in einer Implementierung auch "geht nicht"
        zurückgeben."""

        pass


    # Wenn man invertieren kann, kann man auch
    # dividieren:
    
    def __truediv__(self,other):

        """Division der Form self / other"""
        
        if not (type(self) == type(other) or type(other) == int):
            raise TypeError("Der Divisor ist keine Zahl und kein Ringelement.")
        if type(self) == type(other):
            return self * (other ** -1)

        if type(other) == int:
            return self * (self.ring.element(other) ** -1)


    def __rtruediv__(self,other):

        """Division der Form other / self"""

        if not (type(other) == int):
            raise TypeError("Der Divident ist keine Zahl und kein Ringelement.")

        return other * (self ** -1)


 
    def __pow__(self,exponent):

        """Dies ist für a^b, Eingabe: a ** b"""

    
        if not isinstance(exponent,int):
            raise TypeError()
        
        if exponent < 0:
            self = self.invers()
            exponent = -exponent
            
        # Jetzt kommt square-and-multiply
        # Der Exponent wird 2-adisch dargestellt.
        # Die Funktion dafür ist darunter.
        
        a = self
 
        zweiadisch = zwei_adisch(exponent)

        c = self.ring.eins

        for i in range (0,len(zweiadisch)-1):
            if zweiadisch[i] == '1':
                c = a*c
            c = c*c

        # Am Ende muss man noch einmal addieren, ohne zu multiplizieren.
     
        if len(zweiadisch) > 0 and zweiadisch[len(zweiadisch)-1] == '1':
            c = a*c

        return c


        
    @staticmethod        
    def intmult(n,ele):

        if not isinstance(ele,RingElement):
            return TypeError("Zweites Objekt ist kein Ringelement.")

        if not type(n) == int:
            return TypeError("Das erste Objekt ist keine ganze Zahl.")

        zweiadisch = zwei_adisch(n)

        c = ele.ring.null

        for i in range (0,len(zweiadisch)-1):
            if zweiadisch[i] == '1':
                c = c + ele
            c = c + c

        # Am Ende muss man noch einmal addieren, ohne zu multiplizieren.
     
        if (len(zweiadisch) > 0 
                and zweiadisch[len(zweiadisch)-1] == '1'):
            c = c + ele

        return c



    @staticmethod     

    def test(testelement):
		
        """Mit dieser Methode kann getestet werden, ob ein Objekt vom 
        Typ RingElement oder int ist,
        d.h. ob es im mathematischen Sinne ein Ringelement ist."""


        return isinstance(testelement,RingElement) or type(testelement) == int

    
    @staticmethod     

    def ring(ringelement):

        """Unter Eingabe eines Ringelements im mathematischen Sinn
        (d.h. ein Objekt vom Typ RingElement oder int) wird der Ring
        zurückgegeben.
        
        Hiermit kann eine Fallunterscheidung nach RingElement 
        und int vermieden werden."""

        if not (isinstance(ringelement,RingElement) or type(ringelement) == int):
            raise TypeError("Die Eingabe ist kein Ringelement.")

        if type(ringelement) == int:
            return Z

        if isinstance(ringelement,RingElement):
            return ringelement.ring
            


    
class InvertierungsFehler(ArithmeticError):

    """Fehlerklasse
    
    Für Fehler bei Invertierung."""

    def __init__(self,element):

        self.message = element
        super().__init__("{0} ist nicht invertierbar.".format(element.__str__())
)
    


def zwei_adisch(a):

    """Eine natürliche Zahl wird in zwei-adische Darstellung umgewandelt.
    Die Reihenfolge der bits ist "so wie man schreibt",
    d.h. die letzte Stelle gibt an, ob gerade oder ungerade"""


    if not isinstance(a,int):
        raise TypeError()
 
    if a < 0:
        raise RuntimeError("Die Zahl ist negativ.")

    s = ""
    while a != 0:
        bit = a % 2
        a = a//2
        s = str(bit) +s
    return s




class Ganzzahlring(Ring):

    """Instanziierbare Klasse Ganzzahlring
    Da der Typ int schon existiert, gibt es (leider)
    keine Klasse GanzzahlringElement."""


    def __init__(self): 
        
        self.null, self.eins = 0,1
    
        self._frier()
        
         
    def __str__(self):
        return "Z"


    def element(self,ele):
        if not type(ele) == int:
            raise TypeError("Element ist keine ganze Zahl.")
        return ele

    
    
    @staticmethod

    def ext_ggT(a, b):

        """Berechnung von ggT(a,b) und u,v mit ua + vb = ggT(a,b):
        'geklaut' von http://www.inf.fh-flensburg.de/lang/krypto/algo/euklid.htm
        Die Methode fuer ElementRestklassenring.invers() ist hiervon adaptiert."""


        if not (isinstance(a,int) and isinstance(b,int)):
            raise TypeError()
            
        u, v, s, t = 1, 0, 0, 1
        while b!=0:
            q=a//b
            a, b = b, a-q*b
            u, s = s, u-q*s
            v, t = t, v-q*t
        return a, u, v


    

# Und dann initialisieren wir gleich mal Z:
Z = Ganzzahlring()



class RingTupel(EinfrierbaresObjekt):

    """Instanziierbare Klasse"""
    
    def __init__(self,eingabe1,eingabe2 = None, ring = None):

        # Es gibt mehrere sinnvolle Eingaben:
        
        # eingabe1 : Liste (list) von Ringelementen, möglicherweise eingabe2 oder ring: Ring
        # eingabe1 : Tupel (tuple) von Ringelementen, möglicherweise eingabe2 oder ring: Ring
        # Hier kann "Ringelement" entweder eine Instanz von RingElement
        # oder von int sein.

        # eingabe1 : RingTupel, möglicherweise eingabe2 : Ring

        # alles drei soll dann ein Tupel von Elementen in dem Ring ergeben

        # standardmäßig ist der Ring durch das erste Ringelement gegeben
        
        # eingabe1 : Ringelement, eingabe2 : nicht-negative ganze Zahl,
        #    möglicherweise ring : Ring
        # Dies soll ein Tupel der Länge eingabe2 mit einem Eintrag eingabe1
        # als Element von Ring sein.
        # Standardäßig ist der Ring durch das Ringelement gegeben.

        if isinstance(eingabe2,Ring) and ring == None:

            ring = eingabe2
            eingabe2 = None


        # Jetzt erstmal die Ringelemente kopieren:
#        eingabe1 = deepcopy(eingabe1)
#       unnötig

        self.ring = ring
            
        # Jetzt ist nur noch sinnvoll:
        # eingabe1: Liste oder Tupel von Ringelementen oder RingTupel,
        #                          eingabe2 : None,
        # oder:
        # eingabe1: Ringelement, eingabe2 : int

        # Ferner ist nun self.ring entweder None oder eben der angegebene Ring.

        # Jetzt geht's in zwei Schritten:
        # 1. Eine Liste oder ein Tupel koeffizienten wird initiiert.
        # 2. Aus dieser Liste wird mit dem richtigen Ring das Tupel hergestellt.

        # 1. Schritt:
        if eingabe2 == None:
        # Nur eine Eingabe (bis auf möglicherweise den Ring).
        # Dann sollte diese eine Liste oder ein Tupel oder vom Typ RingTupel sein.
            
            if type(eingabe1) == list:
            # Die eine Eingabe ist eine Liste.
            # Dann wird sie einfach übertragen.
                koeffizienten = eingabe1
            
            elif type(eingabe1) == tuple:
            # Die eine Eingabe ist ein Tupel.
            # Dann wird es in eine Liste verwandelt und übertragen.
            
                koeffizienten = []
                for i in range(0,len(eingabe1)):
                    koeffizienten += [eingabe1[i]]


            elif type(eingabe1) == RingTupel:
                koeffizienten = eingabe1.koeffizienten
                    
            else:
                raise TypeError("Es wurde keine Länge angegeben. Deshalb wurde das erste Objekt als eine Liste oder ein Tupel von Ringelementen (RingElement oder int) oder vom Typ RingTupel erwartet. Dies war nicht der Fall.")


        elif type(eingabe2) == int:
        # Jetzt wird eine Liste der Länge eingabe2 vom Element eingabe1 erzeugt

            if eingabe2 < 0:
                RuntimeError("Die angegebene Länge ist negativ.")

            koeffizienten = []
            for i in range(0,eingabe2):
                koeffizienten += [eingabe1]

                
        else:
            raise RuntimeError("Zu diesen Daten kann kein Objekt vom Typ RingTupel erzeugt werden.")            
            
        # 2. Schritt
        # Jetzt gibt es koeffizienten
            
        self.laenge = len(koeffizienten)
        self.koeffizienten = [] # Hier kommen dann die Koeffizienten rein.
        
        # Wenn die Länge gleich Null ist, ist Schluss.
        if self.laenge == 0:
            if self.ring == None:
                self.ring = Z

        else:
        # Jetzt ist koeffizienten ein nicht-leeres Tupel.
            
        # Die Einträge sollten entweder vom Typ int
        # oder von einem Typ abgeleitet von RingElement sein.
            for i in range(0,self.laenge):
                if not RingElement.test(koeffizienten[i]):
                    raise TypeError("Nicht alle Koeffizienten sind Ringelemente.")


            # Wenn nun self.ring noch nicht vorhanden ist, wird es über
            # den ersten Koeffizienten definiert.

            if self.ring == None:

                if isinstance(koeffizienten[0],RingElement):
                    self.ring = koeffizienten[0].ring
                else: # D.h. koeffizienten[0] ist eine ganze Zahl
                    self.ring = Z
            
            # Jetzt können die Koeffizienten in diesen Ring abgebildet werden.
            # INFO (koeffizient is a PolynomRingElement)
            for i in range(0,self.laenge):
                self.koeffizienten.append(self.ring.element(koeffizienten[i]))


        self._frier()

        
    def __str__(self):
        
        rueck = str(self.laenge) + "-er Tupel von Elementen aus dem Ring {0}: \n\n(".format(self.ring) 
        for i in range(0,self.laenge):

            if self.ring == Z:
                rueck = rueck + str(self.koeffizienten[i])

            else:
                rueck = rueck + self.koeffizienten[i].drucke_element()

            if not i == self.laenge-1:
                rueck = rueck + ","
        
        rueck = rueck + ")"

        return rueck
                
    
    def __eq__(self,other):
        
        if not type(self) == type(other):
            return False
        
        return self.koeffizienten == other.koeffizienten


    def __getitem__(self,index):

        if not type(index) == int:
            raise RuntimeError("Die Eingabe für den Index ist keine ganze Zahl.")
        
        if index < 0 or index >= self.laenge:
            raise RuntimeError("Der Index ist außerhalb des Bereichs.")
    
        return self.koeffizienten[index]

    
    
    def __setitem__(self,index,wert):

        if not type(index) == int:
            raise RuntimeError("Die Eingabe für den Index ist keine ganze Zahl.")
        
        if index < 0 or index >= self.laenge:
            raise RuntimeError("Der Index ist außerhalb des Bereichs.")

        if not RingElement.test(wert):
            raise RuntimeError("Das Objekt auf der rechten Seite ist kein Ringelement.")

        self._tau()
        self.koeffizienten[index] = self.ring.element(wert)
        self._frier()

        
    def __pos__(self):

            return self
        
    def __neg__(self):
        
        rueck = RingTupel(self.ring.null,self.laenge)
        # Das ist das Nulltupel
            
        for i in range(0,self.laenge):
            rueck.koeffizienten[i] = -self.koeffizienten[i]
        
        return rueck

            
    def __add__(self,other):
        
        if not isinstance(other,RingTupel):
            raise TypeError("Der zweite Summand ist kein RingTupel.")
        
        if not self.laenge == other.laenge:
            raise RuntimeError("Die Längen sind nicht gleich.")

        try:
            rueck = RingTupel(other,self.ring)
        except:
            return NotImplemented
            
        for i in range(0,self.laenge):
            rueck[i] += self[i]
        
        return rueck

    
    def __sub__(self,other):
 
        if not isinstance(other,RingTupel):
            return NotImplemented
        
        return self + (-other)
 
    
    def __rmul__(self,other):

        if not RingElement.test(other):
            return NotImplemented
        
        #Wie immer wird hier other*self betrachtet.
                
        # Wir wollen von links mit ganzen Zahlen
        # und mit Ringelementen multiplizieren können.
        # Wenn das Tupel ganze Zahlen enthält, wollen wir ein
        # Tupel von Ringelementen (entsprechend other) erhalten.
        
        if self.ring == Z:

            ring = RingElement.ring(other)

            rueck = RingTupel(ring.null,self.laenge)

            for i in range(0,self.laenge):
                rueck[i] = self[i]*other

        else:

            # other wird nach self.ring abgebildet:
            try:
                multiplikator = self.ring.element(other)
            except:
                return NotImplemented

            rueck = RingTupel(self.ring.null,self.laenge)
            
            for i in range(0,self.laenge):
                rueck[i] = multiplikator*self[i]
                                        
        return rueck


    def auslaufende_nullen_loeschen(self):

        self._tau()

        while self.laenge > 0 and self.koeffizienten[-1] == self.ring.null:
            del self.koeffizienten[-1]
            self.laenge -= 1

        self._frier()

        
    @staticmethod
    def zusammenfuegen(ringtupel1,ringtupel2):

        if not (isinstance(ringtupel1,RingTupel) and isinstance(ringtupel2,RingTupel)):
            raise TypeError() 

        if ringtupel1 == RingTupel([]):
            return ringtupel2
        
        if ringtupel2 == RingTupel([]):
            return(ringtupel1)
        
        if not ringtupel1.ring == ringtupel2.ring:
            raise TypeError("Die Ringe der Tupel sind nicht vergleichbar.")

        return RingTupel(ringtupel1.koeffizienten+ringtupel2.koeffizienten)

