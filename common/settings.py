import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import gridspec
from matplotlib import gridspec as gs
import numpy as np
import pandas as pd

def set_arial_as_font():
    font = 'Arial'
    family = 'sans-serif'

    matplotlib.rcParams['font.{}'.format(family)] = font
    matplotlib.rcParams['font.family'] = family

class Colors:
    def normalize_rgb(dictionary):
        for color in dictionary.keys():
            r, g, b = dictionary[color]
            dictionary[color] = (r / 255., g / 255., b / 255.)
        return dictionary

    almanac = {
        'red': (167, 12, 54)
    }

    google = {
        'red': (255, 19, 1),
        'light red 3': (244, 204, 204), 'light red 2': (234, 153, 153), 'light red 1': (224, 102, 102),
        'dark red 1': (204, 13, 1), 'dark red 2': (153, 7, 0), 'dark red 3': (102, 3, 0),

        'yellow': (255, 253, 2),
        'light yellow 3': (255, 242, 204), 'light yellow 2': (255, 229, 153), 'light yellow 1': (255, 217, 102),
        'dark yellow 1': (241, 194, 50), 'dark yellow 2': (191, 144, 0), 'dark yellow 3': (127, 96, 0),

        'green': (4, 252, 1),
        'light green 3': (217, 234, 211), 'light green 2': (182, 215, 168), 'light green 1': (147, 196, 125),
        'dark green 1': (106, 168, 79), 'dark green 2': (56, 118, 30), 'dark green 3': (39, 78, 19),

        'grey': (204, 204, 204), 'white': (255, 255, 255),
        'light grey 3': (243, 243, 243), 'light grey 2': (239, 239, 239), 'light grey 1': (217, 217, 217),
        'dark grey 1': (183, 183, 183), 'dark grey 2': (153, 153, 153), 'dark grey 3': (102, 102, 102),
    }
    
    tableau10 = {
        'blue': (78, 121, 167), 'orange': (242, 142, 43), 'red': (225, 87, 89),
        'cyan': (118, 183, 178), 'green': (89, 161, 79), 'yellow': (237, 201, 72),
        'purple': (176, 122, 161), 'pink': (225, 157, 167), 'brown': (156, 117, 95),
        'grey': (186, 176, 172), 'white': (240, 240, 240)}
    
    grey = {
        0: (121, 112, 110), 1: (153, 143, 140), 2: (186, 176, 172), 3: (230, 230, 230) 
    }

    greengradient = {
        0: (230, 254, 230), 1: (184, 228, 184), 2: (139, 203, 139), 3: (94, 177, 95),
        4: (51, 152, 52), 5: (15, 127, 18)
    }

    almanac = normalize_rgb(almanac)
    google = normalize_rgb(google)
    tableau10 = normalize_rgb(tableau10)
    grey = normalize_rgb(grey)
    greengradient = normalize_rgb(greengradient)
    
    comut_figure = [tableau10['white'], tableau10['blue'], tableau10['orange']]
