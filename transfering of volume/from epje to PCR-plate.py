# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:05:07 2023

@author: rdebeer
"""
#VARIABLES TO SET!!!
##=============================================================================
number_of_samples = 24                #enter the number of samples 
volume_of_sample = 30                 #enter the volume of the sample
starting_tip = 'A1'                   #enter the starting tip of either the p20 or p200 tips

##=============================================================================
# IMPORT STATEMENTS
#==============================================================================
#Import opentrons protocol API v2
from opentrons import protocol_api

#Import other modules
import math                         #math to do some calculations (rounding up)

##=============================================================================
# METADATA
##=============================================================================
metadata = {
    'protocolName': 'from epje to PCR-plate.py',
    'author': 'RB <rob.de.beer@nioz.nl>',
    'description': ('transfering samples from 1.5ml to a PCR-plate'),
    'apiLevel': '2.12'}

def run(protocol: protocol_api.ProtocolContext):
    """
    pooling replicate PCR reactions into one of the reaction tubes
    
    """

##=============================================================================
# LOADING LABWARE AND PIPETTES
##=============================================================================
    #Loading pipettetips   
    if volume_of_sample > 20:
        tips_1 = protocol.load_labware(
                'opentrons_96_filtertiprack_200ul',  
                11,                                  
                '200tips_1')
        tips_2 = protocol.load_labware(
                'opentrons_96_filtertiprack_200ul',  
                10,                                  
                '200tips_2')
    else:    
        tips_1 = protocol.load_labware(
                'tipone_96_tiprack_20ul',  
                11,                                  
                'tipone_20tips_1')
        tips_2 = protocol.load_labware(
                'tipone_96_tiprack_20ul',  
                10,                                  
                'tipone_20tips_2')
   
    #Loading labware - destination plate
    destination = protocol.load_labware(
                'biorad_96_wellplate_200ul_pcr',
                6,
                'sample_plate_96')
    
    #Loading labware - source
    #Calculating how many sample racks are needed
    sample_racks = math.ceil(number_of_samples/24)
    if sample_racks >= 1:
        sample_tubes_1 = protocol.load_labware(
            'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
            5,                                                       
            'sample_tubes_1')                                        
        if sample_racks >= 2:
            sample_tubes_2 = protocol.load_labware(
                'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
                8,                                                       
                'sample_tubes_2')                                        
            if sample_racks >= 3:
                sample_tubes_3 = protocol.load_labware(
                    'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
                    4,                                                       
                    'sample_tubes_3')                                        
                if sample_racks >= 4:
                    sample_tubes_4 = protocol.load_labware(
                        'opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap',
                        7,                                                      
                        'sample_tubes_4')
    
    #Loading pipettes
    if volume_of_sample >= 20:
        pipette = protocol.load_instrument(
            'p300_single_gen2',             
            'right',                        
            tip_racks=[tips_1, tips_2])
    else:
        pipette = protocol.load_instrument(
            'p20_single_gen2',                  
            'left',                             
            tip_racks=[tips_1, tips_2])      
      
##=============================================================================
# SETTING LOCATIONS============================================================
##=============================================================================
    #Setting starting tip
    pipette.starting_tip = tips_1.well(starting_tip)     
    
    
    PCR1_rows = destination.rows_by_name()
    PCR1_wells = []
    for row in PCR1_rows:
        for well in row:
            PCR1_wells.append(well)
            
    protocol.command(PCR1_wells)         