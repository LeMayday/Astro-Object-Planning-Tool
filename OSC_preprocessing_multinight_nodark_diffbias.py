# script to calibrate and stack multiple nights of astrophotography data 
# Author: Mayday
# Date: 8/20/2023
# starter code from https://siril.org/tutorials/pysiril/

import sys
import os
import tkinter as tk
from tkinter import filedialog
from pysiril.siril   import *
from pysiril.wrapper import *

###############################################################################
# initial file structure:
# [work directory]
# -> S0i (replace i with a number starting at 1 for each night)
#   -> flats
#   -> lights
#   -> biases
# =============================================================================
# final file structure:
# [work directory]
# -> S0i (repace i with a number starting at 1 for each night)
#   -> flats
#   -> lights
#   -> biases
# -> process
#   -> S0i
#       -> flats
#           master_flat.fit
#       -> lights
#           light.seq
#       -> biases
#           master_bias.fit
#   pp_light_group.seq
# result.tif
###############################################################################

# define command blocks for creating masters and processing lights
def master_bias(input_dir, output_dir):
    siril.cd(input_dir)
    siril.convert('bias' , out=output_dir, fitseq=True)
    siril.cd(output_dir)
    # default Winsorized is used if omitted
    siril.stack('bias', type='rej', sigma_low=3, sigma_high=3, norm='no')
    
def master_flat(input_dir, output_dir):
    siril.cd(input_dir)
    siril.convert('flat', out=output_dir, fitseq=True)
    siril.cd(output_dir)
    siril.calibrate('flat', bias='bias_stacked')
    siril.stack('pp_flat', type='rej', sigma_low=3, sigma_high=3, norm='mul')
    
def calibrate_lights(input_dir, output_dir):
    siril.cd(input_dir)
    siril.convert('light', out=output_dir, fitseq=True)
    siril.cd(output_dir)
    siril.calibrate('light', bias='bias_stacked', flat='pp_flat_stacked',
                    cfa=True, equalize_cfa=True, debayer=True)

# register and stack full sequence of calibrated light frames
def stack_lights(input_dir, output_dir):
    siril.cd(input_dir)
    siril.register('pp_light_group')
    siril.stack('r_pp_light_group', type='rej', sigma_low=3, sigma_high=3, norm='addscale', 
                output_norm=True, out=output_dir + '/result')
    siril.close()
    
###############################################################################
# create Siril instance
app = Siril()

# folder chooser
# how to keep filedialog on top: https://stackoverflow.com/a/74643784
# https://blog.filestack.com/thoughts-and-knowledge/python-file-picker-best-practices/
window = tk.Tk()
window.withdraw()
window.attributes('-topmost', True)
workdir = filedialog.askdirectory(mustexist=True, parent=window, title="Select Work Directory")
while not workdir:
    workdir = filedialog.askdirectory(mustexist=True, parent=window, title="Select Work Directory")

# https://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python
with open(workdir + '/log.txt', 'w') as sys.stdout:
    try:
        # start wrapper, Siril, and set preferences
        siril = Wrapper(app)
        app.Open()
        siril.set16bits()
        siril.setext('fit')

        # set folder structure
        process_dir = workdir + '/process'
        # bias_process_dir is the output directory for converted biases, not for input
        # bias_process_dir = process_dir + '/biases'
        session_dirs = []

        # https://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
        for folder in next(os.walk(os.path.join(workdir, '.')))[1]:
            if folder[0] == 'S' and folder[1:2].isdigit():
                session_dirs.append(folder)

        # prepare master bias
        # master_bias(workdir + '/biases', bias_process_dir)

        # prepare master flats and calibrate light frames for each session
        light_seqs = []
        for session in session_dirs:
            master_bias(workdir + '/' + session + '/biases', process_dir + '/' + session)
            master_flat(workdir + '/' + session + '/flats', process_dir + '/' + session)
            calibrate_lights(workdir + '/' + session + '/lights', process_dir + '/' + session)
            light_seqs.append(process_dir + '/' + session + '/pp_light')
        
        # merge documentation says that the final item is the output sequence name
        light_seqs.append('pp_light_group')
        siril.cd(process_dir)
        # https://docs.python.org/3/tutorial/controlflow.html#tut-unpacking-arguments
        siril.merge(*light_seqs)
        stack_lights(process_dir, workdir)

    except Exception as e:
        print("\n**** ERROR *** " + str(e) + "\n" )    

    # close Siril and delete Siril instance
    app.Close()
    del app
    print("Done.")