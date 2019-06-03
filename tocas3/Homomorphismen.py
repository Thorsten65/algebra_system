from GanzzahlRestklassenringe import *
from Polynomringe import *


class Homomorphismus(EinfrierbaresObjekt):

    def __init__ (self,quelle,ziel,element = None,basishom = None):

        # erstmal Attribute setzen:
        self.quelle = quelle
        self.ziel = ziel
        self.element = element
        self.basishom = basishom
    
        if not (isinstance(self.quelle,Ring) and isinstance(self.ziel,Ring)):
            raise TypeError("Die ersten beiden Objekte müssen Ringe sein.")


        # zwei kanonische Homomorphismen:
        # 1. Isomorphismus bei als gleich geltenden Ringen:
        # 2. Die Einbettung R -> R[x]:

        if (element == None and basishom == None) and \
                         (self.quelle == self.ziel or
                         (type(self.ziel) == Polynomring and self.ziel.basisring == self.quelle)):

            def homfunktion(e):
                if self.quelle == Z and type(e) == int:
                    pass
                elif isinstance(self.quelle,Ring) and isinstance(e,RingElement) and e.ring == self.quelle:
                    pass
                else:
                    raise RuntimeError( \
                            "Das Element passt nicht zur Quelle des Homomorphismus.")
                return ziel.element(e)

        else:
        
            # Ansonsten erstmal checken, dass das überhaupt anwendbar ist:
        
            if not (type(self.quelle) == Ganzzahlring or type(self.quelle) == GanzzahlRestklassenring or \
                                     type(self.quelle) == Polynomring):
                raise RuntimeError("Homomorphismus ist für diesen Typ von Ringen nicht implementiert.")

            # Fallunterscheidung nach type(quelle):
            # - Ganzzahlring
            # - Restklassenring
            # - Polynomring

            if type(self.quelle) == Ganzzahlring:
        
                if self.element != None or self.basishom != None:
                    raise RuntimeError("Der erste Ring ist Z, deshalb wird nur noch ein weiterer Ring benötigt.")

                def homfunktion(e):
                    if not type(e) == int:
                        raise TypeError("Die Eingabe muss eine ganze Zahl sein.")
                    return self.ziel.element(e)

        
            elif type(self.quelle) == GanzzahlRestklassenring:

                if self.element != None or basishom != None:
                    raise RuntimeError("Der erste Ring ist ein Ganzzahl-Restklassenring, deshalb wird neben dem zweiten Ring nichts weiter benötigt.")

                if not ((type(self.ziel) == Ganzzahl-Restklassenring and \
                                    self.quelle.modulus % self.ziel.modulus == 0) or \
                     (type(self.ziel) == Polynomring and self.ziel.basisring == self.quelle)):
                    raise RuntimeError("Konstruktion des Homomorphismus nicht möglich.")
                
                def homfunktion(e):
                    if not isinstance(e,RingElement) or e.ring != self.quelle:
                        raise RuntimeError("Das Element passt nicht zur Quelle des Homomorphismus.")
                    return self.ziel.element(e)

            elif type(self.quelle) == Polynomring:
                if element == None:
                    raise RuntimeError("Da die Quelle ein Polynomring ist, muss noch ein Ringelement angegeben werden.")
                if not ((self.ziel == Z and type(self.element) == int) or \
                       (isinstance(self.element,RingElement) and self.element.ring == self.ziel)):
                    raise RuntimeError("Das dritte Objekt muss ein Element aus dem zweiten Ring sein.")

            # Manchmal gibt es kanonische Homomorphismen vom Grundring in das Ziel (basishom)
            # Wenn basisring und ziel als gleich gelten oder
            # wenn quelle und ziel als gleich gelten 
            # (d.h.: ziel gleich quelle.basisring[x]) oder 
            # wenn der Basisring Z oder ein Ganzzahl-Restklassenring ist,
            # können wir das schon definierte Homomorphismus benutzen.
            # Genau diese Homomorphismen sind Beginn implementiert worden.
            
                if self.basishom == None:
                    if self.quelle.basisring == self.ziel or self.quelle == self.ziel or self.quelle.basisring == Z or isinstance(self.quelle.basisring,GanzzahlRestklassenring):
                        self.basishom = Homomorphismus(self.quelle.basisring,self.ziel)
                
                    else:
                        raise RuntimeError("Es wird noch ein Ringhomomorphismus vom Grundring zum Ziel benötigt.")

                if not isinstance(self.basishom,Homomorphismus):
                    raise TypeError("Das vierte Objekt muss ein Ringhomomorphismus (Typ Homomorphismus oder KringelHomomorphismus) sein.")

                if not (self.basishom.quelle == self.quelle.basisring and self.basishom.ziel == self.ziel):
                    raise RuntimeError("Der Grundring-Homomorphismus passt nicht zum Basisring der Quelle oder nicht zum Ziel.")
                
                def homfunktion(e):
                    if not isinstance(e,RingElement) or e.ring != quelle:
                        raise RuntimeError("Das Element passt nicht zur Quelle des Homomorphismus.")
                    rueck = self.ziel.null
                    exp_element = self.ziel.eins
                    for i in range(0,e.grad+1):
                        rueck += self.basishom.anwenden(e.koeffizient(i))*exp_element
                        exp_element = exp_element*self.element
                    return rueck

        # Dies wird auf jeden Fall durchgeführt:
        self.anwenden = homfunktion

        self._frier()


    def __call__(self,other):
        return self.anwenden(other)
		
        
    def __eq__(self,other):

        if not type(other) == type(self):
            return False

        return (self.quelle == other.quelle and self.ziel == other.ziel and \
            self.element == other.element and self.basishom == other.basishom)

    
    def __str__(self):

        rueck = "Ringhomomorphismus von {0} zu {1}".format(self.quelle,self.ziel)
        if (self.element != None and isinstance(self.quelle,Polynomring)):
            if type(self.element) == int:
                ele = self.element
            else:
                ele = self.element.drucke_element()
            rueck = rueck + "\nzu {0} |-> {1}".format(self.quelle.variablenname,ele)
            
            if not self.basishom == None:
                rueck = rueck + "\nund zum Homomorphismus auf dem Grundring:\n{0}".format(self.basishom)

        else:
            rueck = "Kanonischer " + rueck

        return rueck

    # Verknüfung von hom's
    
    def __mul__(self,other):

        return KringelHomomorphismus(self,other)


#------------------------------------------------
# Instanziierbare Klasse KringenHomomorphismus

# Hiermit können Homomorphismen vernüpft werden.
# Man beachte, dass dies zu einer sehr schwachen
# Vergleichbarkeit führt:
# Z.B.: Für phi : R -> R führt
# (phi * phi) * phi == phi * (phi * phi)
# immer zu "False".

class KringelHomomorphismus(Homomorphismus):

    def __init__(self,hom1,hom2):

        if not (isinstance(hom1,Homomorphismus) and isinstance(hom2,Homomorphismus)):
            raise TypeError("Eines der Objekte ist kein Ringhomomorphismus.")
        if not hom1.quelle == hom2.ziel:
            raise RuntimeError("Die Homomorphismen können nicht verknüpft werden.")

        self.hom1 = deepcopy(hom1)
        self.hom2 = deepcopy(hom2)

        self.quelle = self.hom2.quelle
        self.ziel = self.hom1.ziel
        
    def anwenden(self,e):
        
        return self.hom1.anwenden(self.hom2.anwenden(e))

    
    def __str__(self):

        return "Verknüpfte Ringhomomorphismen von {0} nach {1}".format(self.quelle,self.ziel)

    
    # eine sehr rudimentäre Vergleichbarkeit:
    
    def __eq__(self,other):

        if not type(self) == type(other):
            return False
        return self.hom1 == other.hom1 and self.hom2 == other.hom2

    # (Die Verknüpfung "*" wird von Homomorphismus übernommen.)

    
