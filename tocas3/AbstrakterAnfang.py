from abc import ABC, abstractmethod


def typ_beschreibung(objekt):
	
    """Diese Funktion gibt zu einer Instanz einer Klasse den Namen der Klasse als String zurück.
    Anders ausgedrückt: Es gibt den Namen (!) des Typs zurück
    Denn: type(objekt) ist identisch zu objekt.__class__ """

    return objekt.__class__.__name__




class MeinABCObjekt(ABC):

    """Abstrakte Klasse
    
    Der Grund für diese rudimentäre abstrakte Klasse ist dieser:
    
    __str__ wird mit print(.) ausgeben.
    __repr__ wird ausgegeben, wenn man einfach
    den Objektnamen eingibt.
    Nun wird von __repr__ auf __str__ verwiesen."""


    @abstractmethod
    def __init__(self):
        pass

          
    @abstractmethod
    def __str__(self):
        pass
    
    # Hier wird also immer dasselbe ausgegeben.
       
    def __repr__(self):
   
        return self.__str__()
   



class EinfrierbaresObjekt(MeinABCObjekt):

    """Abstrakte Klasse
	
	Die folgende Idee zum "Einfrieren" von Klassen habe ich hier "geklaut"
	https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init
	
	Zum Einfrieren benutzt man dann das Kommando self._frier().
	Auftauen kann man dann mit self._tau()"""
	


    @abstractmethod
    def __init__(self):
        pass
       
    _gefroren = False

    def __setattr__(self, key, value):
        if self._gefroren and not key == "_gefroren":
            raise RuntimeError("{0}-Objekt ist eingefroren.".format(TypBeschreibung(self)))
        object.__setattr__(self, key, value)

    def _frier(self):
        self._gefroren = True

    def _tau(self):
        self._gefroren = False

