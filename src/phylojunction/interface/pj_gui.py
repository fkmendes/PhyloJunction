"""pj_gui.py: PhyloJunction's GUI"""

# NOTE: running the GUI from a conda environment messes up the fonts

from genericpath import isfile
import os
import typing as ty
import re
import PySimpleGUI as sg # type: ignore
import pyperclip # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import pandas as pd
import pickle # type: ignore
# from screeninfo import get_monitors
# import pyautogui
from tabulate import tabulate # type: ignore
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg # type: ignore
from matplotlib.figure import Figure # type: ignore
from matplotlib.widgets import Slider, Button, RadioButtons

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.inference.revbayes.rb_inference as rbinf
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.readwrite.pj_read as pjread
import phylojunction.plotting.pj_plot as pjplot
import phylojunction.interface.cmd.cmd_parse as cmd
import phylojunction.utility.exception_classes as ec
import phylojunction.data.tree as pjdt

__author__ = "Fabio K. Mendes"
__email__ = "f.mendphylojunctiones@wustl.edu"

########################
# Deprecated functions #
########################
# def _add_to_clipboard(text):
#     import tempfile
#     with tempfile.NamedTemporaryFile("w") as fp:
#         fp.write(text)
#         fp.flush()
#         command = "xclip < {}".format(fp.name)
#         os.system(command)

# def _draw_figure(canvas, figure):
#     figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
#     figure_canvas_agg.draw()
#     figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
#     return figure_canvas_agg

# def _delete_fig_agg(fig_agg):
#     fig_agg.get_tk_widget().forget()
#     plt.close('all')

# def _get_scaling():
    # # called before window created
    #     root = sg.tk.Tk()
    #     scaling = root.winfo_fpixels('1i')/72
    #     root.destroy()

    #     return scaling

###############
# Entry point #
###############
def call_gui():
    
    sg.theme("LightGray2")
    
    ###################
    # Inner functions #
    ###################
    def parse_cmd_lines(cmd_line_list, cmd_line_hist, a_pgm, wdw):
        valid_cmd_line = None
        
        for line in cmd_line_list:
            line = line.strip() # removing whitespaces from left and right
            
            try:
                valid_cmd_line = cmd.cmdline2pgm(a_pgm, line)
        
            except Exception as e:
                if not event == "Simulate":
                    if e.__context__:
                        wdw["-ERROR_MSG-"].Update(e.__context__) # __context__ catches the innermost exception message
                    else:
                        wdw["-ERROR_MSG-"].Update(e)

            if valid_cmd_line:
                cmd_line_hist.append(valid_cmd_line)
                wdw["-HIST-"].update("\n".join(cmd_line_hist))
                wdw["-PGM-NODES-"].update([node_name for node_name in pgm_obj.node_name_val_dict])
                wdw["-COMPARISON-PGM-NODES-"].update([nd.node_name for nd in pgm_obj.get_sorted_node_pgm_list() if nd.is_sampled])


    def initialize_axes(fig: Figure):
        ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        return ax


    def clean_disable_everything(cmd_line_hist, msg):
        """
        Destroy PGM object, reset all panels (except history panel).
        History panel is updated with 'reset PGM' line.
        """
        pgm_obj = pgm.ProbabilisticGraphicalModel() # new PGM

        window["-COPY-VALUE-"].update(disabled=True)
        window["-COPY-ALL-"].update(disabled=True)
        window["-ITH-REPL-"].update(1, values=[1], disabled=True)
        window["-ITH-SAMPLE-"].update(1, values=[1], disabled=True)
        window["-ALL-SAMPLES-"].update(True, disabled=True)
        window["-ONE-SAMPLE-"].update(False, disabled=True)
        window["-PGM-NODES-"].update([])
        window["-PGM-NODE-DISPLAY-"].update("")
        window["-PGM-NODE-STAT-"].update("")
        window["-COMPARISON-PGM-NODES-"].update([])
        window["-COMPARISON-PGM-NODE-STATS-"].update([])
        
        cmd_line_hist.append("\n" + msg + "\n")
        window["-HIST-"].update("\n".join(cmd_line_hist))

        node_display_fig.clf() # clear everything
        ax = initialize_axes(node_display_fig)
        node_display_fig.canvas.draw() # updates canvas

        comparison_fig.clf() # clear everything
        comparison_ax = initialize_axes(comparison_fig)
        comparison_fig.canvas.draw() # updates canvas

        return pgm_obj, ax, comparison_ax


    def draw_node_pgm(axes, node_pgm, sample_idx=None, repl_idx=0, repl_size=1):
        return node_pgm.plot_node(axes, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)

    
    def selected_node_read(pgm_obj, node_name):
        node_pgm = pgm_obj.get_node_pgm_by_name(node_name)
        # display_node_pgm_value_str = pg.get_display_str_by_name(node_name)
        sample_size = len(node_pgm) # this is n_sim inside sampling distribution classes
        repl_size = node_pgm.repl_size
        
        return node_pgm, sample_size, repl_size

    
    def selected_node_display(wdw, pgm_obj, node_pgm, do_all_samples, sample_idx=None, repl_idx=0, repl_size=1):
        display_node_pgm_value_str = str()
        display_node_pgm_stat_str = str()

        # first we do values
        # we care about a specific sample and maybe a specific replicate
        if not sample_idx == None and not do_all_samples:
            start = sample_idx * repl_size
            end = start + repl_size
            display_node_pgm_value_str = node_pgm.get_start2end_str(start, end) # values
            display_node_pgm_stat_str = node_pgm.get_node_stats_str(start, end, repl_idx) # summary stats
        
        # we get all samples
        else:
            # just calling __str__
            display_node_pgm_value_str = pgm_obj.get_display_str_by_name(node_pgm.node_name)
            # getting all values
            display_node_pgm_stat_str = node_pgm.get_node_stats_str(0, len(node_pgm.value), repl_idx) # summary stats
        
        wdw["-PGM-NODE-DISPLAY-"].update(display_node_pgm_value_str)
        wdw["-PGM-NODE-STAT-"].update(display_node_pgm_stat_str)

        
    def selected_node_plot(fig_obj, node_pgm, do_all_samples, sample_idx=None, repl_idx=0, repl_size=1):
        """
        Plot pgm node on 'node_display_fig_axes' (Axes object) scoped to 'call_gui()',
        then update canvas with new plot
        """
        try:
            # if a tree
            if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                draw_node_pgm(node_display_fig_axes, node_pgm, sample_idx=sample_idx, repl_idx=repl_idx)
            # when not a tree
            else:
                if do_all_samples: 
                    draw_node_pgm(node_display_fig_axes, node_pgm, repl_size=repl_size)
                else:
                    draw_node_pgm(node_display_fig_axes, node_pgm, sample_idx=sample_idx, repl_size=repl_size)
        # when it's deterministic
        except:
            draw_node_pgm(node_display_fig_axes, node_pgm)

        fig_obj.canvas.draw()

    
    def do_selected_node(pgm_obj, wdw, fig_obj, node_name, do_all_samples=True):
        """
        Given selected node name, display its string representation and
        plot it on canvas if possible
        """
        node_pgm, sample_size, repl_size = selected_node_read(pgm_obj, node_name)

        # updates spin window with number of elements in this node_pgm
        # window["-ITH-VAL-"].update(values=[x for x in range(1, sample_size + 1)]) # can only select the number of values this node contains
        wdw["-ITH-SAMPLE-"].update(values=[x for x in range(1, sample_size + 1)]) # can only select the number of values this node contains

        if type(node_pgm.value) == list:
            if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                wdw["-ITH-REPL-"].update(values=[x for x in range(1, repl_size + 1)]) # can only select the number of values this node contains
            else:
                wdw["-ITH-REPL-"].update(disabled=True)
        
        # makes sure we stay at the right spin element
        # sample_idx = int(window["-ITH-VAL-"].get()) - 1 # (offset)
        sample_idx = int(wdw["-ITH-SAMPLE-"].get()) - 1 # (offset)
        repl_idx = int(wdw["-ITH-REPL-"].get()) - 1 # (offset)

        # updating node values on window happens inside
        selected_node_display(wdw, pgm_obj, node_pgm, do_all_samples, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)
    
        # plotting to canvas happens inside
        selected_node_plot(fig_obj, node_pgm, do_all_samples, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)

        return node_pgm

    ###################### 
    # Development screen #
    ######################
    # call sg.Window.get_screen_size()
    # my_width, my_height = 3840, 3240 # LG monitor home on Ubuntu
    # my_width, my_height = 1440, 900 # macbook pro built-in
    # my_width, my_height = 1920, 1080 # LG monitor office
    # my_width, my_height = 1366, 768 # dell xps 13
    # scaling_old = _get_scaling() # 1.0 when I call it on 1440, 900
    # print("scaling_old = " + str(scaling_old))
    
    ##############
    # New screen #
    ##############
    # alternative
    # pyautogui.size()
    # ws, hs = list(), list()
    # for m in get_monitors():
    #     ws.append(m.width)
    #     hs.append(m.height)
    # width = max(ws)
    # height = max(hs)
    # width, height = sg.Window.get_screen_size() # alternative (will return built-in dimensions)
    # print("width = " + str(width) + " height = " + str(height))

    # scaling = scaling_old * min(width / my_width, height / my_height) 
    # print("scaling = " + str(scaling))

    # scaling should not be necessary -- I will keep the code
    # above just in case...
    # sg.set_options(scaling=scaling)






    ####################################
    # ~+~+~+~+~ Layout START ~+~+~+~+~ #
    ####################################
    
    menu_def = [
        ["&File", ["&Open script", "&Load model", "---", "&Save data to", "Save &model to", "Save &inference specs. to", "Save &history as", "Save &figure as", "---", "E&xit"]]
    ]

    # gray96
    node_value_layout = [
        [
            sg.Multiline(key="-PGM-NODE-DISPLAY-", font=("Courier", 14), disabled=True, background_color="gray96", size=(102,10))
        ]
    ]

    node_stat_layout = [
        [
            sg.Multiline(key="-PGM-NODE-STAT-", font=("Courier", 14), disabled=True, background_color="gray96", size=(102,10))
        ]
    ]


    layout_main = [

        [ sg.Menu(menu_def, key="-TOP-MENU-", font=("helvetica", 14)) ],

        [ sg.Push(),

          sg.Frame( 
            layout = [
                [
                    sg.TabGroup(
                        [
                            [ sg.Tab("Value(s)", node_value_layout, font=("Helvetica", 14)) ],
                            [ sg.Tab("Summary stats.", node_stat_layout, font=("Helvetica", 14)) ]
                        ],
                        font=("Helvetica", 14),
                    ),

                    sg.VPush(),

                    sg.Column(
                        layout = [
                            [ sg.Radio("All samples", "-WHICH-SAMPLES-", key="-ALL-SAMPLES-", default=True, enable_events=True, disabled=True, font=("helvetica 14")) ],
                            [ sg.Radio("One sample", "-WHICH-SAMPLES-", key="-ONE-SAMPLE-", default=False, enable_events=True, disabled=True, font=("helvetica 14")) ],
                            [ sg.Text("Sample #", font=("helvetica 14")) ],
                            [ sg.Spin([1], initial_value=1, key="-ITH-SAMPLE-", enable_events=True, readonly=True, disabled=True, font=("helvetica 14"), size=(8,1)) ],
                            [ sg.Text("Replicate #", font=("helvetica 14"), size=(8,1)) ],
                            [ sg.Spin([1], initial_value=1, key="-ITH-REPL-", enable_events=True, readonly=True, disabled=True, font=("helvetica 14"), size=(8,1)) ],
                            [ sg.VPush() ],
                            [ sg.Button("Copy values", key="-COPY-VALUE-", size=(10,1), disabled=True, disabled_button_color=("gray70", "gray85"), font=("helvetica 14")) ],
                            [ sg.Checkbox("Copy all", key="-COPY-ALL-", default=False, font=("helvetica 14")) ]
                        ],
                        element_justification="center",
                        expand_y="true" # necessary for VPush() to have an effect
                    )           
                ]
            ],
            title="SELECTED NODE",
            title_location="n",
            relief=sg.RELIEF_SOLID,
            border_width=1,
            font="helvetica 12 bold"
          ),

          sg.Push()
        ],

        [ sg.Push(),
        
          sg.Canvas(key="-CANVAS-", background_color="white", size=(1200,450)),
          
          sg.Push(),

          sg.Frame(
            layout = [
                [
                    sg.Column(
                        # gray96 for Listbox
                        layout = [
                            [ sg.Listbox([], key="-PGM-NODES-", font=("helvetica", 16), enable_events=True, background_color="mint cream", size=(15,16)) ],
                            [ sg.Button("Reset", key="-DESTROY-PGM-", enable_events=True, disabled=True, disabled_button_color=("gray70", "gray85"), font=("helvetica 14")) ]
                        ],
                        element_justification="center",
                        background_color="mint cream",
                    )
                ]
            ],
            title="MODEL",
            title_color="steel blue",
            title_location="n",
            relief=sg.RELIEF_FLAT,
            background_color="mint cream",
            font="helvetica 16 bold"
          ),

          sg.Push()
        ],

        [ sg.Frame(
            layout = [
                [ sg.Input(key="-CMD-", font=("Helvetica", 16), text_color="white", background_color="gray20", focus=True, do_not_clear=False, size=(111,1))  ]
            ],
            title="Command prompt",
            relief=sg.RELIEF_FLAT,
            # pad=(5,20),
            font="helvetica 14"
            )
        ],

        [ sg.Button("Enter command", visible=False, bind_return_key=True) ],

        [ sg.Frame(
            layout = [
                [ sg.Multiline(key="-ERROR_MSG-", font="helvetica 16", text_color="red", background_color="gray96", disabled=True, size=(110,3)) ]
            ],
            title="Log",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ]
    ]
    

    layout_history = [
        [ sg.Frame(
            layout = [
                [ sg.Multiline(key="-HIST-", font=("helvetica", 16), disabled=True, background_color="gray96", size=(110,40)) ]
            ],
            title="Command history",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ]
    ]

    ################
    # Layout infer #
    ################
    layout_infer = [
        [
            sg.Frame(
            layout = [
                [ sg.Multiline(key="-INFERENCE-SPEC-", font=("helvetica", 16), disabled=True, background_color="gray96", size=(110,38)) ]
            ],
            title="Inference specification",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ],
        [
            sg.Text("Inference platform", font=("Helvetica", 14)),
            sg.Combo(["RevBayes", "BEAST 2"], key="-PROGRAM-", font=("Helvetica", 14)),
            sg.Text("Chain length", font=("Helvetica", 14)), sg.Input(size=(10, 1), default_text="1000000", font=("Helvetica", 14), key="-CHAIN-"),
            sg.Button("See", key="-GEN-INF-", font=("Helvetica", 14)),
            sg.Spin([x for x in range(1, 1)], 1, key="-ITH-INF-", enable_events=True, readonly=True, disabled=True, font=("helvetica 14"), size=(8,1))
            # sg.InputText(key="-INFERENCE-SAVE-TO-DUMMY-", enable_events=True, visible=False), # dummy element to trigger folder browse
            # sg.FolderBrowse(button_text="Choose folder", initial_folder="./", key="-INFERENCE-DIR-", size=(13,1), font=("Helvetica", 14)),
            # sg.Button("Save specification(s)", key="-INFERENCE-SAVE-", size=(12,1), font=("Helvetica", 14))
            # sg.InputText(key="-SCRIPT-SAVE-DUMMY-", enable_events=True, visible=False),
            # sg.FileSaveAs(button_text="Save history as", initial_folder="./", key="-SCRIPT-SAVE-AS-", font=("Helvetica", 14), button_color=("gold")),
        ]
    ]

    ##################
    # Layout compare #
    ##################
    layout_compare = [
        [
            
            sg.Column(
                # gray96 for Listbox
                layout = [
                    [ sg.VPush() ],
                    [ 
                        sg.Frame(
                            layout = [
                                [
                                    sg.Listbox([], key="-COMPARISON-PGM-NODES-", font=("helvetica", 16), enable_events=True, size=(15,16), background_color="linen")
                                ],
                            ],
                            title="SAMPLED NODES",
                            title_location="n",
                            relief=sg.RELIEF_FLAT,
                            border_width=1,
                            font="helvetica 12 bold",
                            background_color="antiquewhite"
                        )
                    ],
                    [ sg.Checkbox("Replicate avgs.", key="-REPL-AVG-", default=False, font=("helvetica 14")) ],
                    [ sg.VPush() ],
                    [ 
                        sg.Frame(
                            layout = [
                                [
                                    sg.Listbox([], key="-COMPARISON-PGM-NODE-STATS-", font=("helvetica", 16), enable_events=True, size=(15,16))
                                ],
                            ],
                            title="SUMMARY STAT.",
                            title_location="n",
                            relief=sg.RELIEF_FLAT,
                            border_width=1,
                            font="helvetica 12 bold"
                        )
                    ],
                    [ sg.VPush() ]
                ],
                element_justification="left",
                expand_y="true" # necessary for VPush() to have an effect
            ),

            sg.Column(
                layout = [
                    [ sg.Button("Compare to .csv (...)", key="-LOAD-COMPARISON-CSV-", font=("Helvetica", 14)) ],
                    [ sg.Multiline(key="-COMPARE-TO-", font=("Courier", 12), disabled=True, background_color="gray96", size=(108,18)) ],
                    [ sg.Canvas(key="-COMPARISON-CANVAS-", background_color="white", size=(1200,800)) ],
                    [
                        sg.Button("Draw", key="-DRAW-VIOLIN1-", font=("Helvetica", 14)),
                        sg.Button("Save plot to", key="-SAVE-VIOLIN-TO-", font=("Helvetica", 14))
                    ]
                ],
                element_justification="center"
            )
        ]
    ]


    ###################
    # Layout validate #
    ###################
    layout_validate = [
        []
    ]


    #####################
    # Layout save specs #
    #####################
    layout_save = [
        [
            sg.Text("Save with prefix", font=("Helvetica", 14)), sg.Input(size=(10, 1), font=("Helvetica", 14), key="-PREFIX-")
        ]
    ]

    ################
    # Top tab menu #
    ################
    tabgrp = [
        [
            sg.TabGroup(
                [
                    [ sg.Tab("PGM", layout_main, font=("Helvetica", 14)) ],
                    [ sg.Tab("History", layout_history, font=("Helvetica", 14)) ],
                    [ sg.Tab("Infer", layout_infer, font=("Helvetica", 14)) ],
                    [ sg.Tab("Compare", layout_compare, font=("Helvetica", 14)) ],
                    [ sg.Tab("Validate", layout_validate, font=("Helvetica", 14)) ],
                    [ sg.Tab("Save specs.", layout_save, font=("Helvetica", 14)) ]
                ],
                font=("Helvetica", 14)
            )
        ]
    ]
    
    ##################################
    # ~+~+~+~+~ Layout END ~+~+~+~+~ #
    ##################################





    ###############
    # Main window #
    ###############
    
    # margins=(0,0)
    window = sg.Window("PhyloJunction", tabgrp, finalize=True, resizable=True, keep_on_top=False, element_justification="c", location=(0,0)) 
    window['-CMD-'].Widget.config(insertbackground="white")
    window["-REPL-AVG-"].set_tooltip("test")
    # window.Size = (1050, 760) # resolution-dependent

    # Main screen figure
    fig_agg = None
    node_display_fig = Figure(figsize=(11,4.5))
    node_display_fig_axes = initialize_axes(node_display_fig)

    # Link matplotlib to PySimpleGUI Graph
    fig_agg = FigureCanvasTkAgg(node_display_fig, window["-CANVAS-"].TKCanvas)
    plot_widget = fig_agg.get_tk_widget()
    plot_widget.grid(row=0, column=0)

    # Compare figure
    fig_agg2 = None
    comparison_fig = Figure(figsize=(11,4.5))
    comparison_fig_axes = initialize_axes(comparison_fig)

    # Link matplotlib to PySimpleGUI Graph
    fig_agg2 = FigureCanvasTkAgg(comparison_fig, window["-COMPARISON-CANVAS-"].TKCanvas)
    plot_widget2 = fig_agg2.get_tk_widget()
    plot_widget2.grid(row=0, column=0)

    ##############
    # Event loop #
    ##############
    pgm_obj: pgm.ProbabilisticGraphicalModel = pgm.ProbabilisticGraphicalModel()
    cmd_history: ty.List[str] = []
    cmd_line: str = ""
    script_out_fp: str = ""
    data_out_dir: str = ""
    inf_out_dir: str = ""

    # for comparison tab
    comparison_csv_fp: str = ""
    scalar_output_stash: ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]] = []
    tree_output_stash: ty.List[ty.Dict[str, pd.DataFrame]] = []
    pj_comparison_df: pd.DataFrame = pd.DataFrame()
    other_comparison_df: pd.DataFrame = pd.DataFrame()
    
    while True:
        # general parameters for window
        event, values = window.read()

        #######################
        # Try to read command #
        #######################
        try:
            cmd_line = values["-CMD-"]
        # user closed window without typing command line
        except:
            pass


        ########
        # Quit #
        ########
        if cmd_line in ("q()", "quit", "exit", "close", "bye"):
            event = "Quit"
        
        if event in (None, "Quit", "Quit2", "Exit"):
            break


        #######################################################
        # Copying value of selected node in "Created node(s)" #
        #######################################################
        elif event == "-COPY-VALUE-":
            value_str = str()
            
            if values["-COPY-ALL-"]:
                value_str = pgm_obj.get_display_str_by_name(node_pgm.node_name)
            else:
                value_str = values['-PGM-NODE-DISPLAY-']
            
            pyperclip.copy(str(value_str)) # requires xclip (Linux) / pbcopy (OS X) are installed

        elif event == "Simulate":
            pass


        ####################
        # Destroying model #
        ####################
        elif event == "-DESTROY-PGM-":
            # pgm_obj and node_display_fig_axes are overwritten with new clean objects
            # cmd_history is updated inside
            pgm_obj, node_display_fig_axes, comparison_fig_axes = clean_disable_everything(cmd_history, "## Reset taking place at this point. Previous model, if it existed, is now obliterated") # clean figure, reset axes


        ####################################
        # Select node in "Created node(s)" #
        ####################################
        elif event == "-PGM-NODES-" or event == "-ALL-SAMPLES-" or event == "-ONE-SAMPLE-":

            # if nodes have been created and selected
            if values["-PGM-NODES-"]:
                selected_node_pgm_name = values["-PGM-NODES-"][0]
                do_all_samples = window["-ALL-SAMPLES-"].get() # True or False

                # if selected node is tree, we do not want to show all trees on display by default
                try:
                    if isinstance(pgm_obj.get_node_pgm_by_name(selected_node_pgm_name).value[0], pjdt.AnnotatedTree):
                        do_all_samples = False
                except: pass # the value of the node_pgm might be an MacroevolStateDependentRateParameter, which is not subscriptable, so we pass
                
                node_pgm = do_selected_node(pgm_obj, window, node_display_fig, selected_node_pgm_name, do_all_samples=do_all_samples)
                
                # we enable value copying as soon as a node is clicked
                window["-COPY-VALUE-"].update(disabled=False)
                window["-COPY-ALL-"].update(disabled=False)

                # if there is a chance for replicates to exist, we enable the one-sample radio button
                if node_pgm.is_sampled:
                    window["-ONE-SAMPLE-"].update(disabled=False)

                    # cycling through trees can only be done with "one-sample" radio button
                    if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                        window["-ALL-SAMPLES-"].update(disabled=True)
                        window["-ALL-SAMPLES-"].update(False)
                        window["-ONE-SAMPLE-"].update(True)
                        window["-ITH-SAMPLE-"].update(disabled=False)
                        window["-ITH-REPL-"].update(disabled=False)
                    # for all other stochastic nodes, cycling can be done through "all-samples"
                    else:
                        window["-ALL-SAMPLES-"].update(disabled=False)

                # if node was created by hand through assignment, radio buttons and cycling are disabled
                else:
                    window["-ALL-SAMPLES-"].update(disabled=True)
                    window["-ONE-SAMPLE-"].update(disabled=True)

            # if we're looking at all samples, all samples and replicates shown in histogram
            # and if tree, only one sample at a time is allowed
            if event == "-ALL-SAMPLES-":
                window["-ITH-REPL-"].update(disabled=True)
                window["-ITH-SAMPLE-"].update(disabled=True)
            
            # if one sample at a time, spin element works (only for sampled random variables)
            if event == "-ONE-SAMPLE-":
                window["-ITH-SAMPLE-"].update(disabled=False)

                # if we are looking at trees, we can cycle through replicates
                if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                    window["-ITH-REPL-"].update(disabled=False)
                # otherwise, all replicates will be visualized as histogram (no cycling allowed)
                else:
                    window["-ITH-REPL-"].update(disabled=True)


        #########################################
        # Going through scalars with replicates #
        #########################################
        elif event == "-ITH-SAMPLE-":
            # if nodes have been created and selected
            if values["-PGM-NODES-"]:
                selected_node_pgm_name = values["-PGM-NODES-"][0]
                do_all_samples = window["-ALL-SAMPLES-"].get() # True or False
                node_pgm = do_selected_node(pgm_obj, window, node_display_fig, selected_node_pgm_name, do_all_samples=do_all_samples)


        #######################
        # Going through trees #
        #######################
        elif event == "-ITH-REPL-":
            # # TODO: fix visualization of tree depending on radio button
            # sample_idx = values["-ITH-SAMPLE-"] - 1 # (offset)
            # repl_idx = values["-ITH-REPL-"] - 1 # (offset)
            
            # # only updates display if tree node is selected
            # if isinstance(node_pgm.value[0], AnnotatedTree):
            #     draw_node_pgm(node_display_fig_axes, node_pgm, sample_idx=sample_idx, repl_idx=repl_idx)
            #     node_display_fig.canvas.draw()

            # if nodes have been created and selected
            if values["-PGM-NODES-"]:
                selected_node_pgm_name = values["-PGM-NODES-"][0]
                do_all_samples = window["-ALL-SAMPLES-"].get() # True or False
                node_pgm = do_selected_node(pgm_obj, window, node_display_fig, selected_node_pgm_name, do_all_samples=do_all_samples)


        ##################
        # Reading script #
        ##################
        elif event == "Open script":
            _script_in_fp = sg.popup_get_file("Script to load", no_window=True, keep_on_top=True)
            
            if _script_in_fp:
                _, node_display_fig_axes, comparison_fig_axes = clean_disable_everything(cmd_history, "## Loading script at this point. Previous model, if it existed, is now obliterated") # clean figure, reset axes
                cmd_lines = pjread.read_text_file(_script_in_fp)
                parse_cmd_lines(cmd_lines, cmd_history, pgm_obj, window)

                # it model has at least one node, destroying is now a possibility
                if pgm_obj.n_nodes >= 1:
                    window["-DESTROY-PGM-"].update(disabled=False)
            
            else: pass # event canceled by user


        ###############
        # Loading PGM #
        ###############
        elif event == "Load model":
            _model_in_fp = sg.popup_get_file("Model to load", no_window=True, keep_on_top=True)
            
            if _model_in_fp:
                _, node_display_fig_axes = clean_disable_everything(cmd_history, "## Loading serialized model at this point. Previous model, if it existed, is now obliterated") # clean figure, reset axes
                pgm_obj = pjread.read_serialized_pgm(_model_in_fp)
                window["-PGM-NODES-"].update([node_name for node_name in pgm_obj.node_name_val_dict])
                window["-COMPARISON-PGM-NODES-"].update([nd.node_name for nd in pgm.get_sorted_node_pgm_list() if nd.is_sampled])
            
            else: pass # user canceled event


        ###############
        # Saving data #
        ###############
        elif event == "Save data to":
            try:
                data_out_dir = sg.popup_get_folder("Save to directory", no_window=True, keep_on_top=True)
            
            except: pass # event canceled by user
            
            prefix = values["-PREFIX-"]
                        
            # default value is "./"
            # if data_out_dir:
            pjwrite.dump_pgm_data(data_out_dir, pgm_obj, prefix=prefix)


        ###########################
        # Saving serialized model #
        ###########################
        elif event == "Save model to":
            try:
                model_out_dir = sg.popup_get_folder("Save to directory", no_window=True, keep_on_top=True)
            
            except: pass # event canceled by user
            
            prefix = values["-PREFIX-"]
            
            # pickling and saving PGM
            pjwrite.dump_serialized_pgm(model_out_dir, pgm_obj, prefix=prefix)


        ##################
        # Saving history #
        ##################
        elif event == "Save history as":
            try:
                script_out_fp = sg.popup_get_file("Save to file", save_as=True, no_window=True, keep_on_top=True) # directory for saving data files
                # script_out_fp = values['-SCRIPT-SAVE-DUMMY-'] # file path
            except: pass # canceled

            if script_out_fp:
                with open(script_out_fp, "w") as script_out:
                    pjwrite.write_text_output(script_out, cmd_history)


        #################
        # Saving figure #
        #################
        elif event == "Save figure as":
            try:
                fig_out_fp = sg.popup_get_file("Save figure as", save_as=True, no_window=True, keep_on_top=True) # directory for saving figure
            except: pass # canceled

            if fig_out_fp:
                pjwrite.write_fig_to_file(fig_out_fp, node_display_fig)




        # =-=-= INFERENCE =-=-= #

        #############################
        # Generating inference spec #
        #############################
        elif event == "-GEN-INF-":
            if pgm_obj.n_nodes >= 1:
                try:
                    mcmc_chain_length = int(values["-CHAIN-"])
                
                    if mcmc_chain_length == 0:
                        raise ec.InvalidMCMCChainLength("MCMC chain cannot have size zero.")
                except ec.InvalidMCMCChainLength as e:
                    sg.popup_error("Invalid MCMC chain length")

                all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = rbinf.pgm_obj_to_rev_inference_spec(pgm_obj, inf_out_dir, mcmc_chain_length=mcmc_chain_length)

                all_sims_spec_strs_list = pjwrite.get_write_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, write2file=False)

                # enabling and updating spin element
                window["-ITH-INF-"].update(values=[x for x in range(1, len(all_sims_model_spec_list) + 1)])
                window["-ITH-INF-"].update(1)
                window["-ITH-INF-"].update(disabled=False)

                window["-INFERENCE-SPEC-"].update(all_sims_spec_strs_list[0])
            
            else: pass


        #####################################
        # Rotating through simulation specs #
        #####################################
        elif event == "-ITH-INF-":
            sample_idx = int(window["-ITH-INF-"].get()) - 1 # (offset)
            window["-INFERENCE-SPEC-"].update(all_sims_spec_strs_list[sample_idx])


        #########################
        # Saving inference spec #
        #########################
        # elif event == "-INFERENCE-SAVE-":
        #     prefix = values["-PREFIX-"]
        #     try:
        #         mcmc_chain_length = int(values["-CHAIN-"])
                
        #         if mcmc_chain_length == 0:
        #             raise ec.InvalidMCMCChainLength("MCMC chain cannot have size zero.")
            
        #     except ec.InvalidMCMCChainLength as e:
        #         sg.popup_error("Invalid MCMC chain length")
            
        #     except Exception as e:
        #         sg.popup_error("MCMC chain length must be a number")

        #     all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = pjinf.pgm_obj_to_rev_inference_spec(pgm_obj, inf_out_dir, mcmc_chain_length=mcmc_chain_length, prefix=prefix)

        #     _ = pjio.get_write_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, prefix=prefix, write2file=True)


        #############################
        # Save-to inference scripts #
        #############################
        elif event == "Save inference specs. to":
            inf_out_dir = sg.popup_get_folder("Save to directory", no_window=True, keep_on_top=True) # directory for saving data files

            prefix = values["-PREFIX-"]
            try:
                mcmc_chain_length = int(values["-CHAIN-"])
                
                if mcmc_chain_length == 0:
                    raise ec.InvalidMCMCChainLength("MCMC chain cannot have size zero.")
            
            except ec.InvalidMCMCChainLength as e:
                sg.popup_error("Invalid MCMC chain length")
            
            except Exception as e:
                sg.popup_error("MCMC chain length must be a number")

            all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = rbinf.pgm_obj_to_rev_inference_spec(pgm_obj, inf_out_dir, mcmc_chain_length=mcmc_chain_length, prefix=prefix)

            _ = pjwrite.get_write_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, prefix=prefix, write2file=True)
            



        # =-=-= COMPARISON =-=-= #
        
        ###########################
        # Reading comparison .csv #
        ###########################
        elif event == "-LOAD-COMPARISON-CSV-":
            comparison_csv_fp = sg.popup_get_file("CSV file to load", no_window=True, keep_on_top=True)
            
            if comparison_csv_fp:
                if not os.path.isfile(comparison_csv_fp):
                    print("Could not find " + comparison_csv_fp) # TODO: add exception here later

                try:
                    other_comparison_df = pjread.read_csv_into_dataframe(comparison_csv_fp)
                    other_comparison_df_str = tabulate(other_comparison_df, other_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip()
                    window["-COMPARE-TO-"].update(other_comparison_df_str)
                except:
                    print("Could not load .csv file into pandas DataFrame. Later write an Exception for this")
            
            else: pass # event canceled by user


        elif event == "-COMPARISON-PGM-NODES-":
            selected_node_name: str = ""
            avg_over_repls = values["-REPL-AVG-"]
            
            if values["-COMPARISON-PGM-NODES-"]:
                selected_node_name = values["-COMPARISON-PGM-NODES-"][0]
                selected_node = pgm_obj.get_node_pgm_by_name(selected_node_name)
                selected_node_repl_size = selected_node.repl_size

                if not scalar_output_stash or not tree_output_stash:
                    scalar_output_stash, tree_output_stash = pjwrite.prep_data_df(pgm_obj)
                    _, scalar_value_df_dict, scalar_repl_summary_df = scalar_output_stash
                    _, tree_summary_df_dict, tree_repl_summary_df_dict = tree_output_stash
                    
                    # adding "program" column to PJ's pandas DataFrame's stored inside dicts
                    for repl_size in scalar_value_df_dict:
                        scalar_value_df_dict[repl_size].loc[:, "program"] = "PJ"

                    for tr_nd_name in tree_summary_df_dict:
                        tree_summary_df_dict[tr_nd_name].loc[:, "program"] = "PJ"
                    
                    for tr_w_repl_nd_name in tree_repl_summary_df_dict:
                        tree_repl_summary_df_dict[tr_w_repl_nd_name].loc[:, "program"] = "PJ"

                if isinstance(selected_node.value[0], (int, float, np.float64)):
                    if not avg_over_repls:
                        window["-COMPARISON-PGM-NODE-STATS-"].update([])
                        pj_comparison_df = scalar_value_df_dict[selected_node_repl_size]

                    else:
                        window["-COMPARISON-PGM-NODE-STATS-"].update(values=["average", "std. dev."])
                        pj_comparison_df = scalar_repl_summary_df

                elif isinstance(selected_node.value[0], pjdt.AnnotatedTree):
                    pj_comparison_df = tree_summary_df_dict[selected_node_name]
                    window["-COMPARISON-PGM-NODE-STATS-"].update(values=[tree_stat for tree_stat in pj_comparison_df.keys() if tree_stat not in ("program", "sample", "replicate")])
            
            else:
                print("If there are no sampled nodes in PJ model, we cannot compare against PJ. Later write an Exception for this")
                pass

        
        elif event == "-DRAW-VIOLIN1-":
            avg_over_repls = values["-REPL-AVG-"]
            comparison_fig.clf()
            comparison_fig_axes = initialize_axes(comparison_fig)
            comparison_fig.canvas.draw() # updates canvas

            if not pj_comparison_df.empty and not other_comparison_df.empty and values["-COMPARISON-PGM-NODES-"]:
                value_to_compare: str = ""
                
                # scalar
                if isinstance(selected_node.value[0], (int, float, np.float64)):
                    if not avg_over_repls:
                        value_to_compare = values["-COMPARISON-PGM-NODES-"][0]
                    
                    # averaging over replicates
                    else:
                        if values["-COMPARISON-PGM-NODE-STATS-"]:
                            value_to_compare = values["-COMPARISON-PGM-NODE-STATS-"][0] 

                        # summary stats were not picked, doing nothing
                        else:
                            pass
                
                # custom class object (e.g., AnnotatedTree)
                else:
                    if values["-COMPARISON-PGM-NODE-STATS-"]:
                        value_to_compare = values["-COMPARISON-PGM-NODE-STATS-"][0]
                    
                    # summary stats were not picked, doing nothing
                    else:
                        pass

                # debugging
                # print(tabulate(pj_comparison_df, pj_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())
                # print(tabulate(other_comparison_df, other_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())
                
                joint_dataframe = pjplot.join_dataframes(pj_comparison_df, other_comparison_df, value_to_compare=value_to_compare, summaries_avg_over_repl=False)

                # something went wrong, we clear comparison figure
                if joint_dataframe.empty: pass
                
                else:
                    pjplot.pjdraw.plot_violins(comparison_fig, comparison_fig_axes, joint_dataframe, "program", value_to_compare, xlab="Program", ylab=value_to_compare)
                
                # debugging
                # print(tabulate(joint_dataframe, joint_dataframe.head(), tablefmt="pretty", showindex=False).lstrip())

            else:
                print("I need both PJ's and the other simulator dataframes to plot violins.")

        
        elif event == "-SAVE-VIOLIN-TO-":
            try:
                compare_fig_out_fp = sg.popup_get_file("Save plot as", save_as=True, no_window=True, keep_on_top=True) # directory for saving figure
            
            except: pass # canceled

            if compare_fig_out_fp:
                pjwrite.write_fig_to_file(compare_fig_out_fp, comparison_fig)
        
        #######################################
        # If no other event, we parse command #
        #######################################
        else: 
            cmd_lines = [cl for cl in re.split("\n", cmd_line) if cl] # in case the user enters multiple lines, also ignores empty lines

            parse_cmd_lines(cmd_lines, cmd_history, pgm_obj, window)

            # it model has at least one node, we can destroy it
            if pgm_obj.n_nodes >= 1:
                window["-DESTROY-PGM-"].update(disabled=False)


    window.close()







# call GUI
if __name__ == "__main__":
    call_gui()