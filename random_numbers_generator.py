# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:38:25 2024

@author: ChatGPT
"""


# Hoeveel getallen heb je?
totaal = 163
# Hoeveel random getallen wil je?
aantal = 10


import random

def generate_unique_random_numbers():
    # Genereer 10 unieke willekeurige getallen tussen 1 en 163
    unique_random_numbers = random.sample(range(1, totaal + 1), aantal)
    return unique_random_numbers

# Roep de functie aan om de lijst met unieke willekeurige getallen te verkrijgen
unique_random_numbers_list = generate_unique_random_numbers()

# Druk de resultaten af
print(f"Lijst met {aantal} unieke willekeurige getallen:", unique_random_numbers_list)