# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:38:25 2024

@author: ChatGPT
"""

import random

# Hoeveel getallen heb je?
totaal = 163
# Hoeveel random getallen wil je?
aantal = 10

def generate_random_numbers():
    # Genereer 10 willekeurige getallen tussen 1 en 163
    random_numbers = [random.randint(1, totaal) for i in range(aantal)]
    return random_numbers

# Roep de functie aan om de lijst met willekeurige getallen te verkrijgen
random_numbers_list = generate_random_numbers()

# Druk de resultaten af
print("Lijst met willekeurige getallen:", random_numbers_list)