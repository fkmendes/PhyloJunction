"""pj_gui.py: PhyloJunction's GUI"""

# NOTE: running the GUI from a conda environment messes up the fonts

import os
import typing as ty
import re
import PySimpleGUI as sg # type: ignore
import pyperclip # type: ignore
import matplotlib.pyplot as plt # type: ignore
# from screeninfo import get_monitors
# import pyautogui
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg # type: ignore
from matplotlib.figure import Figure # type: ignore

# pj imports
import pgm.pgm as pgm
import inference.revbayes.rb_inference as rbinf
import readwrite.pj_write as pjwrite
import readwrite.pj_read as pjread
import interface.cmd.cmd_parse as cmd
import utility.exception_classes as ec
import data.tree as pjdt

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

def _add_to_clipboard(text):
    import tempfile
    with tempfile.NamedTemporaryFile("w") as fp:
        fp.write(text)
        fp.flush()
        command = "xclip < {}".format(fp.name)
        os.system(command)

def _draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
    return figure_canvas_agg

def _delete_fig_agg(_fig_agg):
    _fig_agg.get_tk_widget().forget()
    plt.close('all')


###############
# Entry point #
###############
def call_gui():
    
    sg.theme("LightGray2")
    
    def _get_scaling():
    # called before window created
        root = sg.tk.Tk()
        scaling = root.winfo_fpixels('1i')/72
        root.destroy()

        return scaling

    
    def _selected_node_read(pgm_obj, node_pgm_name):
        node_pgm = pgm_obj.get_node_pgm_by_name(node_pgm_name)
        # display_node_pgm_str = pg.get_display_str_by_name(node_pgm_name)
        sample_size = len(node_pgm) # this is n_sim inside sampling distribution classes
        repl_size = node_pgm.repl_size
        
        return node_pgm, sample_size, repl_size

    
    def _selected_node_display(wdw, pgm_obj, node_pgm, do_all_samples, sample_idx=None, repl_size=1):
        display_node_pgm_str = str()

        # we care about sample_idx
        if not sample_idx == None and not do_all_samples:
            start = sample_idx * repl_size
            end = start + repl_size 
            display_node_pgm_str = node_pgm.get_start2end_str(start, end)
        
        # we get all samples, and call __str__
        else:
            display_node_pgm_str =  pgm_obj.get_display_str_by_name(node_pgm.node_pgm_name)

        wdw["-PGM_NODE_DISPLAY-"].update(display_node_pgm_str)

    
    def _selected_node_plot(fig_obj, node_pgm, do_all_samples, sample_idx=None, repl_idx=0, repl_size=1):
        try:
            # if a tree
            if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                _draw_node_pgm(_axes, node_pgm, sample_idx=sample_idx, repl_idx=repl_idx)
            # when not a tree
            else:
                if do_all_samples: 
                    _draw_node_pgm(_axes, node_pgm, repl_size=repl_size)
                else:
                    _draw_node_pgm(_axes, node_pgm, sample_idx=sample_idx, repl_size=repl_size)
        # when it's deterministic
        except: _draw_node_pgm(_axes, node_pgm)

        fig_obj.canvas.draw()

    
    def _do_selected_node(pgm_obj, wdw, fig_obj, node_pgm_name, do_all_samples=True):
        node_pgm, sample_size, repl_size = _selected_node_read(pgm_obj, node_pgm_name)

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
        _selected_node_display(wdw, pgm_obj, node_pgm, do_all_samples, sample_idx=sample_idx, repl_size=repl_size)
    
        # plotting to canvas happens inside
        _selected_node_plot(fig_obj, node_pgm, do_all_samples, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)

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

    ###################
    # Inner functions #
    ###################
    def _parse_cmd_lines(cmd_line_list, cmd_line_hist, _pgm, wdw):
        valid_cmd_line = None
        
        for line in cmd_line_list:
            line = line.strip() # removing whitespaces from left and right
            
            try:
                valid_cmd_line = cmd.cmdline2pgm(_pgm, line)
        
            except Exception as e:
                if not event == "Simulate":
                    if e.__context__:
                        wdw["-ERROR_MSG-"].Update(e.__context__) # __context__ catches the innermost exception message
                    else:
                        wdw["-ERROR_MSG-"].Update(e)

            if valid_cmd_line:
                cmd_line_hist.append(valid_cmd_line)
                wdw["-HIST-"].update("\n".join(cmd_line_hist))
                wdw["-PGM-NODES-"].update([node_name for node_name in _pgm.node_name_val_dict])

    ##############
    # GUI layout #
    ##############
    
    menu_def = [
        ["&File", ["&Load script", "&Save history as", "Save &data to", "Save &inference specs. to", "---", "E&xit"]]
    ]

    layout = [

        [ sg.Menu(menu_def, key="-TOP-MENU-", font=("helvetica", 14)) ],

        [ sg.Push(),

          sg.Frame(
            layout = [
                [
                  sg.Multiline(key="-PGM_NODE_DISPLAY-", font=("Courier", 14), disabled=True, background_color="gray96", size=(102,6))
                ]
            ],
            title="Node value",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            ),

            sg.Push(),

            sg.Frame(
                layout = [
                    [
                    sg.Listbox([], key="-PGM-NODES-", font=("helvetica", 16), enable_events=True, background_color="gray96", size=(12,6))
                    ]
                ],
                title="Created node(s)",
                relief=sg.RELIEF_FLAT,
                font="helvetica 14"
            ),

            sg.Push()
        ],

        [ sg.Push(), sg.Canvas(key="-CANVAS-", background_color="white", size=(1200,450)), sg.Push(),

          sg.Frame(
              layout = [
                    [
                        sg.Column(
                            layout = [
                                [ sg.Checkbox("Copy all", key="-COPY-ALL-", default=False, font=("helvetica 14")) ],
                                [ sg.Button("Copy values", key="-COPY-VALUE-", size=(10,1), disabled=True, disabled_button_color=("gray70", "gray85"), font=("helvetica 14")) ], 
                                [ sg.Radio("All samples", "-WHICH-SAMPLES-", key="-ALL-SAMPLES-", default=True, enable_events=True, disabled=True, font=("helvetica 14")) ],
                                [ sg.Radio("One sample", "-WHICH-SAMPLES-", key="-ONE-SAMPLE-", default=False, enable_events=True, disabled=True, font=("helvetica 14")) ],
                                [ sg.Text("Sample #", font=("helvetica 14")) ],
                                [ sg.Spin([x for x in range(1, 1)], 1, key="-ITH-SAMPLE-", enable_events=True, readonly=True, disabled=True, font=("helvetica 14"), size=(8,1)) ],
                                [ sg.Text("Tree #", font=("helvetica 14"), size=(8,1)) ],
                                [ sg.Spin([x for x in range(1, 1)], 1, key="-ITH-REPL-", enable_events=True, readonly=True, disabled=True, font=("helvetica 14"), size=(8,1)) ]
                            ],
                            element_justification="center"
                        )
                    ]
                ],
              title="Selected node",
              font="helvetica 14"
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
                
                [ sg.Multiline(key="-ERROR_MSG-", font="helvetica 16", text_color="red", background_color="gray96", disabled=True, size=(110,5)) ]

            ],
            title="Log",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ],

        # sg.Input(key="-SCRIPT-READ-", enable_events=True, visible=False), # dummy element to trigger file browse
        # sg.FileBrowse(button_text="Open script file", target="-SCRIPT-READ-", font=("Helvetica", 14), button_color="skyblue1"),
        [
            sg.Text("Save with prefix", font=("Helvetica", 14)), sg.Input(size=(10, 1), font=("Helvetica", 14), key="-PREFIX-")
        ]

        # sg.InputText(key="-VALUES-SAVE-DUMMY-", enable_events=True, visible=False) # dummy element to trigger folder browse
        # sg.FolderBrowse(button_text="Save values in", initial_folder="./", key="-VALUES-SAVE-IN-", font=("Helvetica", 14), button_color=("gold"))
        # sg.Quit(font=("Helvetica", 14), button_color=("white", "firebrick")
    ]
    
    layout2 = [
        [ sg.Frame(
            layout = [
                [ sg.Multiline(key="-HIST-", font=("helvetica", 16), disabled=True, background_color="gray96", size=(110,37)) ]
            ],
            title="Command history",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ]
    ]

    layout3 = [
        [
            sg.Frame(
            layout = [
                [ sg.Multiline(key="-INFERENCE-SPEC-", font=("helvetica", 16), disabled=True, background_color="gray96", size=(110,36)) ]
            ],
            title="Inference specification",
            relief=sg.RELIEF_FLAT,
            font="helvetica 14"
            )
        ],
        [ sg.VPush() ],
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

    tabgrp = [
        [
            sg.TabGroup(
                [
                    [ sg.Tab("PGM", layout, font=("Helvetica", 14)) ],
                    [ sg.Tab("History", layout2, font=("Helvetica", 14)) ],
                    [ sg.Tab("Inference", layout3, font=("Helvetica", 14)) ]
                ],
                font=("Helvetica", 14)
            )
        ]
    ]


    ###############
    # Main window #
    ###############
    window = sg.Window("PhyloJunction", tabgrp, finalize=True, resizable=True, keep_on_top=False, element_justification="c") 
    window['-CMD-'].Widget.config(insertbackground="white")
    # window.Size = (1000, 750)

    _fig_agg = None
    _fig = Figure(figsize=(11,4.5))
    _axes = _fig.add_axes([0.25, 0.2, 0.5, 0.6])
    _axes.patch.set_alpha(0.0)
    _axes.xaxis.set_ticks([])
    _axes.yaxis.set_ticks([])
    _axes.spines['left'].set_visible(False)
    _axes.spines['bottom'].set_visible(False)
    _axes.spines['right'].set_visible(False)
    _axes.spines['top'].set_visible(False)

    # Link matplotlib to PySimpleGUI Graph
    _fig_agg = FigureCanvasTkAgg(_fig, window["-CANVAS-"].TKCanvas)
    plot_widget = _fig_agg.get_tk_widget()
    plot_widget.grid(row=0, column=0)

    def _draw_node_pgm(axes, node_pgm, sample_idx=None, repl_idx=0, repl_size=1):
        return node_pgm.get_gcf(axes, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)

    ##############
    # Event loop #
    ##############
    _pgm = pgm.ProbabilisticGraphicalModel()
    _cmd_history: ty.List[str] = []
    _cmd_line: str = ""
    _script_out_fp: str = ""
    _data_out_dir: str = ""
    _inf_out_dir: str = ""
    while True:
        # general parameters for window
        event, values = window.read()

        #######################
        # Try to read command #
        #######################
        try:
            _cmd_line = values["-CMD-"]
        # user closed window without typing command line
        except:
            pass

        ########
        # Quit #
        ########
        if _cmd_line in ("q()", "quit", "exit", "close", "bye"):
            event = "Quit"
        
        if event in (None, "Quit", "Quit2", "Exit"):
            break

        #######################################################
        # Copying value of selected node in "Created node(s)" #
        #######################################################
        elif event == "-COPY-VALUE-":
            value_str = str()
            
            if values["-COPY-ALL-"]:
                value_str = _pgm.get_display_str_by_name(node_pgm.node_pgm_name)
            else:
                value_str = values['-PGM_NODE_DISPLAY-']
            
            pyperclip.copy(str(value_str)) # requires xclip (Linux) / pbcopy (OS X) are installed

        elif event == "Simulate":
            pass

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
                    if isinstance(_pgm.get_node_pgm_by_name(selected_node_pgm_name).value[0], pjdt.AnnotatedTree):
                        do_all_samples = False
                except: pass # the value of the node _pgm might be an AtomicSSERateParameter, which is not subscriptable, so we pass
                
                node_pgm = _do_selected_node(_pgm, window, _fig, selected_node_pgm_name, do_all_samples=do_all_samples)
                
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
                node_pgm = _do_selected_node(_pgm, window, _fig, selected_node_pgm_name, do_all_samples=do_all_samples)

        #######################
        # Going through trees #
        #######################
        elif event == "-ITH-REPL-":
            # # TODO: fix visualization of tree depending on radio button
            # sample_idx = values["-ITH-SAMPLE-"] - 1 # (offset)
            # repl_idx = values["-ITH-REPL-"] - 1 # (offset)
            
            # # only updates display if tree node is selected
            # if isinstance(node_pgm.value[0], AnnotatedTree):
            #     _draw_node_pgm(_axes, node_pgm, sample_idx=sample_idx, repl_idx=repl_idx)
            #     _fig.canvas.draw()

            # if nodes have been created and selected
            if values["-PGM-NODES-"]:
                selected_node_pgm_name = values["-PGM-NODES-"][0]
                do_all_samples = window["-ALL-SAMPLES-"].get() # True or False
                node_pgm = _do_selected_node(_pgm, window, _fig, selected_node_pgm_name, do_all_samples=do_all_samples)

        ##################
        # Reading script #
        ##################
        elif event == "Load script":
            _pgm = pgm.ProbabilisticGraphicalModel() # resets _pgm
            _script_in_fp = sg.popup_get_file("Script to load", no_window=True, keep_on_top=True)
            
            cmd_lines = pjread.read_text_file(_script_in_fp)

            _parse_cmd_lines(cmd_lines, _cmd_history, _pgm, window)

        # original implmn had buttons
        # script_in_fp = values["-SCRIPT-READ-"]
        # elif event == "-SCRIPT-READ-":

        ###############
        # Saving data #
        ###############
        elif event == "Save data to":
            try:
                _data_out_dir = sg.popup_get_folder("Save to directory", no_window=True, keep_on_top=True)
            except: pass # canceled
            
            prefix = values["-PREFIX-"]
            
            # default value is "./"
            # if _data_out_dir:
            pjwrite.dump_pgm_data(_data_out_dir, _pgm, prefix=prefix)

        # original implmn had buttons
        # _data_out_dir = values["-VALUES-SAVE-DUMMY-"] # directory for saving data files
        # elif event == "-VALUES-SAVE-DUMMY-":

        ##################
        # Saving history #
        ##################
        elif event == "Save history as":
            try:
                _script_out_fp = sg.popup_get_file("Save to file", save_as=True, no_window=True, keep_on_top=True) # directory for saving data files
                # _script_out_fp = values['-SCRIPT-SAVE-DUMMY-'] # file path
            except: pass # canceled

            if _script_out_fp:
                with open(_script_out_fp, "w") as script_out:
                    pjwrite.write_text_output(script_out, _cmd_history)

        # original implmn had a button
        # elif event == "-SCRIPT-SAVE-DUMMY-":

        #############################
        # Generating inference spec #
        #############################
        elif event == "-GEN-INF-":
            if _pgm.n_nodes >= 1:
                try:
                    mcmc_chain_length = int(values["-CHAIN-"])
                
                    if mcmc_chain_length == 0:
                        raise ec.InvalidMCMCChainLength("MCMC chain cannot have size zero.")
                except ec.InvalidMCMCChainLength as e:
                    sg.popup_error("Invalid MCMC chain length")

                all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = rbinf.pgm_obj_to_rev_inference_spec(_pgm, _inf_out_dir, mcmc_chain_length=mcmc_chain_length)

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

        #     all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = pjinf.pgm_obj_to_rev_inference_spec(_pgm, _inf_out_dir, mcmc_chain_length=mcmc_chain_length, prefix=prefix)

        #     _ = pjio.get_write_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, prefix=prefix, write2file=True)

        #############################
        # Save-to inference scripts #
        #############################
        elif event == "Save inference specs. to":
            _inf_out_dir = sg.popup_get_folder("Save to directory", no_window=True, keep_on_top=True) # directory for saving data files

            prefix = values["-PREFIX-"]
            try:
                mcmc_chain_length = int(values["-CHAIN-"])
                
                if mcmc_chain_length == 0:
                    raise ec.InvalidMCMCChainLength("MCMC chain cannot have size zero.")
            
            except ec.InvalidMCMCChainLength as e:
                sg.popup_error("Invalid MCMC chain length")
            
            except Exception as e:
                sg.popup_error("MCMC chain length must be a number")

            all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = rbinf.pgm_obj_to_rev_inference_spec(_pgm, _inf_out_dir, mcmc_chain_length=mcmc_chain_length, prefix=prefix)

            _ = pjwrite.get_write_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, prefix=prefix, write2file=True)
            
        # original implmn had buttons
        # elif event == "-INFERENCE-SAVE-TO-DUMMY-":
        # _inf_out_dir = values["-INFERENCE-SAVE-TO-DUMMY-"]

        #######################################
        # If no other event, we parse command #
        #######################################
        else: 
            cmd_lines = [cl for cl in re.split("\n", _cmd_line) if cl] # in case the user enters multiple lines, also ignores empty lines

            _parse_cmd_lines(cmd_lines, _cmd_history, _pgm, window)


    window.close()

# call GUI
if __name__ == "__main__":
    call_gui()