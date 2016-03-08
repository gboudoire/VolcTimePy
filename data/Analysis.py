# -*- coding: utf-8 -*-
'''
Created on 26 oct. 2015

:copyright:
    Guillaume Boudoire (guillaume.boudoire@gmail.com)
        Laboratoire GeoSciences Réunion
        Université de La Reunion
        Observatoire Volcanologique du Piton de La Fournaise
        Institut de Physique du Globe de Paris
    Patrice Boissier (boissier@ipgp.fr)
        Observatoire Volcanologique du Piton de La Fournaise
        Institut de Physique du Globe de Paris
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
'''


class Analysis(object):
    '''
    classdocs
    '''

    def __init__(self, mineral):
        '''
        Constructor
        '''
        self.measures = []
        self.mineral = mineral
