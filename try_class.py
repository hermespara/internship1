#!/usr/bin/python3.6
class Prelevement:
    def __init__(self):
        self.tissu = None
        self.analysis = None


class Individu:
    def __init__(self):
        self.age = None
        self.death = None
        self.prelevements = []
        self.id = None
GTEX1117F = Individu()

GTEX1117F.age = 13
print(GTEX1117F)
