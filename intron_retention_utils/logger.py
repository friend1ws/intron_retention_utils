#! /usr/bin/env python

import logging

def get_logger():

    logger = logging.getLogger("intron_retention")
    logger.propagate = False
 
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
     
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %H:%M:%S')

    ch.setFormatter(formatter)

    logger.addHandler(ch)

    return logger

