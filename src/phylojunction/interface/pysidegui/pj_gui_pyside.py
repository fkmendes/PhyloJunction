import sys
import os
import re
import typing as ty
import numpy as np
import pandas as pd
from tabulate import tabulate # type: ignore
from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton, QFileDialog
from PySide6.QtGui import QAction
from PySide6.QtCore import QPropertyAnimation, QEasingCurve, QTimer

# pj imports #
from phylojunction.interface.pysidegui.content_main_window \
    import ContentGUIMainWindow
from phylojunction.pgm.pgm import ProbabilisticGraphicalModel
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.plotting.pj_organize as pjorg
import phylojunction.plotting.pj_draw as pjdraw
import phylojunction.readwrite.pj_read as pjread
import phylojunction.readwrite.pj_write as pjwrite
import phylojunction.data.tree as pjdt


class GUIModeling():

    pgm_obj: ProbabilisticGraphicalModel
    cmd_log_list: ty.List[str]

    def __init__(self):
        self.pgm_obj = ProbabilisticGraphicalModel()
        self.cmd_log_list = []


    def parse_cmd_update_pgm(self, cmd_line_list, gui_main_window_obj):
        valid_cmd_line = None
        print("\nReading following command lines:")
        
        for line in cmd_line_list:
            # removing whitespaces from left and right
            line = line.strip()
            print("  " + line)

            # side-effect in cmdline2pgm
            try:
                valid_cmd_line = cmdp.cmdline2pgm(self.pgm_obj, line)
        
            except Exception as e:
                gui_main_window_obj.ui.bottom_label_left.setText("Warning produced")

                # if not event == "Simulate":
                if e.__context__:
                    # __context__ catches innermost exception 
                    # update warnings page later, with e.__context__
                    print(e.__context__)
                    gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(e.__context__))
                else:
                    print(e)  # update warnings page with e
                    gui_main_window_obj.ui.ui_pages.warnings_textbox.setText(str(e))

            if valid_cmd_line:
                self.cmd_log_list.append(valid_cmd_line)
            

    def cmd_log(self):
        return "\n".join(self.cmd_log_list)
    
    def clear(self):
        self.pgm_obj = ProbabilisticGraphicalModel()


class GUIMainWindow(QMainWindow):

    input_script_filepath: str
    gui_modeling: GUIModeling
    is_avg_repl_check: bool
    pj_comparison_df: pd.DataFrame
    other_comparison_df: pd.DataFrame

    def __init__(self):
        super().__init__()

        self.setWindowTitle("PhyloJunction")

        self.init_top_menus()

        # setup main window #
        self.ui = ContentGUIMainWindow()

        # passing myself as argument to auxiliary setup class
        self.ui.setup_ui(self)

        #######################
        # Members for storage #
        # and operations      #
        #######################
        
        self.input_script_filepath = ""

        # model-related members #
        self.gui_modeling = GUIModeling()

        self.is_avg_repl_check = False
        

        #########################
        # Below we specify what #
        # happens when buttons  #
        # are pressed           #
        #########################

        # toggle button #
        self.ui.menu_button.clicked.connect(self.toggle_button_animation)

        # menu button #
        self.ui.pgm_button.clicked.connect(self.show_pgm_page)

        # pgm page button #
        self.ui.cmd_log_button.clicked.connect(self.show_cmd_log_page)

        # compare page button #
        self.ui.compare_button.clicked.connect(self.show_compare_page)

        # coverage page button #
        self.ui.covg_button.clicked.connect(self.show_pgm_page)

        # warnings page button #
        self.ui.warning_button.clicked.connect(self.show_warnings_page)


        ###################################
        # gui_pages buttons               #
        # (created by QtCreator/Designer) #
        ###################################

        # PGM page #
        # cmd line enter #
        self.ui.ui_pages.cmd_prompt.returnPressed.connect(self.parse_cmd_update_gui)

        # node list #
        self.ui.ui_pages.node_list.itemClicked.connect(self.do_selected_node_pgm_page)

        # radio button update #
        # self.ui.ui_pages.one_sample_radio.toggled.connect(lambda: self.refresh_selected_node_display_plot)
        # self.ui.ui_pages.all_samples_radio.toggled.connect(lambda: self.refresh_selected_node_display_plot)
        self.ui.ui_pages.one_sample_radio.clicked.connect(self.refresh_selected_node_display_plot)
        self.ui.ui_pages.all_samples_radio.clicked.connect(self.refresh_selected_node_display_plot)

        # spin boxes update #
        self.ui.ui_pages.sample_idx_spin.valueChanged.connect(self.refresh_selected_node_display_plot)
        self.ui.ui_pages.repl_idx_spin.valueChanged.connect(self.refresh_selected_node_display_plot)

        # clear model button #
        self.ui.ui_pages.clear_model.clicked.connect(self.clean_disable_everything)

        # Compare page #
        self.ui.ui_pages.compare_csv_button.clicked.connect(self.read_compare_csv)
        self.ui.ui_pages.compare_node_list.itemClicked.connect(self.do_selected_node_compare_page)
        self.ui.ui_pages.avg_replicate_check_button.clicked.connect(self.avg_repl_check)
        self.ui.ui_pages.draw_violins_button.clicked.connect(self.draw_violin)

        # show #
        self.show()

    # set up all top menus
    def init_top_menus(self):
        # setup menu bar #
        menubar = self.menuBar()
        
        # first item (pj menu)
        pj_menu = menubar.addMenu("PhyloJunction")

        # pj menu items
        # certain keywords cannot be used in
        # menubars on Mac, like "quit", "about"
        # https://stackoverflow.com/questions/73293306/string-containing-about-cannot-be-added-to-the-menu-in-pyside6
        about_action = pj_menu.addAction("About", self.print_about)
        about_action.setMenuRole(QAction.MenuRole.NoRole)
        pj_menu.addAction(about_action)

        quit_action = pj_menu.addAction("Quit", self.quit_app_button_clicked)
        quit_action.setMenuRole(QAction.MenuRole.NoRole)
        pj_menu.addAction(quit_action)

        # second item (file menu)
        file_menu = menubar.addMenu("File")

        # file menu items
        file_menu.addAction("Read script", self.read_execute_script)

        file_menu.addAction("Load model", self.print_about)
        file_menu.addSeparator()
        file_menu.addAction("Save data as", self.print_about)
        file_menu.addAction("Save model as", self.print_about)


    # verify which menu buttons are selected #
    def reset_selection(self):
        for btn in self.ui.left_menu.findChildren(QPushButton):
            try:
                btn.set_active(False)
            
            except:
                pass


    def toggle_button_animation(self):
        # get menu width #
        menu_width = self.ui.left_menu.width()

        # check width #
        width = 50  # base width
        if menu_width == 50:
            width = 180

        # start animation #
        self.animation = QPropertyAnimation(self.ui.left_menu, b"minimumWidth")
        self.animation.setStartValue(menu_width)
        self.animation.setEndValue(width)
        self.animation.setDuration(500)  # 500 ms
        self.animation.setEasingCurve(QEasingCurve.InOutCirc)
        self.animation.start()


    #########################################
    # Functions for activating page buttons #
    #########################################

    def show_pgm_page(self):
        self.reset_selection()
        # ui_pages is a member of self.ui,
        # which is a ContentGUIMainWindow instance
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.pgm_page)
        self.ui.pgm_button.set_active(True)
        self.ui.top_label_left.setText("MODEL SPECIFICATION")


    def show_compare_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.compare_page)
        self.ui.compare_button.set_active(True)
        self.ui.top_label_left.setText("IMPLEMENTATION COMPARISON")


    def show_cmd_log_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.cmd_log_page)
        self.ui.cmd_log_button.set_active(True)
        self.ui.top_label_left.setText("COMMAND LOG")


    def show_warnings_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.warnings_page)
        self.ui.warning_button.set_active(True)
        self.ui.top_label_left.setText("WARNINGS")
        self.ui.bottom_label_left.setText("")


    ###################################
    # Various event-related functions #
    ###################################

    def parse_cmd_update_gui(self):
        cmd_line = self.ui.ui_pages.cmd_prompt.text()  # grab input
        
        # in case the user enters many lines
        # (also ignores empty lines)
        cmd_line_list = [cl for cl in re.split("\n", cmd_line) if cl]
        
        # read commands
        # (side-effect: gui_modeling stores cmd hist)
        self.gui_modeling.parse_cmd_update_pgm(cmd_line_list, self)
        
        # update GUI cmd history
        self.refresh_cmd_history()

        self.ui.ui_pages.cmd_prompt.clear()  # clear

        # now update ListWidget #
        pgm_node_list = [node_name for \
                node_name in self.gui_modeling.pgm_obj.node_name_val_dict]
        self.ui.ui_pages.node_list.clear()  # clear first
        self.ui.ui_pages.node_list.addItems(pgm_node_list)


    def selected_node_display(self, node_pgm, do_all_samples, sample_idx=None, repl_idx=0, repl_size=1):
        display_node_pgm_value_str = str()
        display_node_pgm_stat_str = str()

        # first we do values
        # we care about a specific sample and maybe a specific replicate
        if not sample_idx == None and not do_all_samples:
            start = sample_idx * repl_size
            end = start + repl_size

            # values
            display_node_pgm_value_str = \
                node_pgm.get_start2end_str(start, end)
            
            # summary stats
            display_node_pgm_stat_str = \
                node_pgm.get_node_stats_str(start, end, repl_idx)
        
        # we get all samples
        else:
            # just calling __str__
            display_node_pgm_value_str = \
                self.gui_modeling.pgm_obj.\
                get_display_str_by_name(node_pgm.node_name)
            
            # getting all values
            display_node_pgm_stat_str = \
                node_pgm.get_node_stats_str(
                    0, len(node_pgm.value), repl_idx) # summary stats
        
        # print("Set values_content QLineEdit widget with text: " + display_node_pgm_value_str)
        # print("Set summary_content QLineEdit widget with text: " + display_node_pgm_stat_str)
        self.ui.ui_pages.values_content.setText(display_node_pgm_value_str)
        self.ui.ui_pages.summary_content.setText(display_node_pgm_stat_str)


    def selected_node_plot(self, fig_obj, fig_axes, node_pgm, do_all_samples,
                           sample_idx=None, repl_idx=0, repl_size=1):
        """
        Plot pgm node on 'node_display_fig_axes' (Axes object) scoped to 'call_gui()',
        then update canvas with new plot
        """
        try:
            # if a tree
            if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                self.draw_node_pgm(fig_axes, node_pgm, sample_idx=sample_idx, repl_idx=repl_idx)
            
            # when not a tree
            else:
                if do_all_samples: 
                    self.draw_node_pgm(fig_axes, node_pgm, repl_size=repl_size)
                else:
                    self.draw_node_pgm(fig_axes, node_pgm, sample_idx=sample_idx, repl_size=repl_size)
        
        # when it's deterministic
        except:
            self.draw_node_pgm(fig_axes, node_pgm)

        fig_obj.canvas.draw()
    

    def selected_node_read(self, node_name):
        node_pgm = self.gui_modeling.pgm_obj.get_node_pgm_by_name(node_name)
        # this is n_sim inside sampling distribution classes
        sample_size = len(node_pgm)
        repl_size = node_pgm.repl_size

        return node_pgm, sample_size, repl_size


    def do_selected_node_pgm_page(self):
        """
        Display selected node's string representation and
        plot it on canvas if possible, for pgm page
        (model specification menu)
        """
        
        # first get selected node's name #
        node_name = \
            self.ui.ui_pages.node_list.currentItem().text()

        # reading node information #
        node_pgm, sample_size, repl_size = \
            self.selected_node_read(node_name)

        # spin boxes must be up-to-date #
        self.ui.ui_pages.repl_idx_spin.setMaximum(repl_size)
        if sample_size == 0:
            self.ui.ui_pages.sample_idx_spin.setMaximum(0)    
        
        else:
            self.ui.ui_pages.sample_idx_spin.setMaximum(sample_size)

        # grab pgm_page's figure and axes #
        fig_obj = self.ui.ui_pages.pgm_page_matplotlib_widget.fig
        fig_axes = self.ui.ui_pages.pgm_page_matplotlib_widget.axes
        
        # activate radio buttons and spin elements #
        self.init_and_refresh_radio_spin(node_pgm, sample_size)
        
        do_all_samples = \
            self.ui.ui_pages.all_samples_radio.isChecked()
        
        
        #################################
        # Collect information from spin #
        #################################
        
        sample_idx = self.ui.ui_pages.sample_idx_spin.value() - 1
        repl_idx = self.ui.ui_pages.repl_idx_spin.value() - 1

        print("sample_idx = " + str(sample_idx))
        print("repl_idx = " + str(repl_idx))

        ###############
        # Now do node #
        ###############
        self.selected_node_display(node_pgm,
                                   do_all_samples,
                                   sample_idx=sample_idx,
                                   repl_idx=repl_idx,
                                   repl_size=repl_size)

        self.selected_node_plot(fig_obj,
                                fig_axes,
                                node_pgm,
                                do_all_samples,
                                sample_idx=sample_idx,
                                repl_idx=repl_idx,
                                repl_size=repl_size)      


    def do_selected_node_compare_page(self):
        """
        Chose node to get summary statistics
        and compare against user implementation
        (implementation comparison menu)
        """

        # first get selected node's name #
        selected_node_name = \
            self.ui.ui_pages.compare_node_list.currentItem().text()

        # if there is at least one
        # sampled node with a name
        if not selected_node_name in ("", None):

            # reading node information #
            selected_node_pgm, selected_node_sample_size, selected_node_repl_size = \
                self.selected_node_read(selected_node_name)

            # could be more efficient, but this
            # makes sure that stashes are always up-to-date
            scalar_output_stash, tree_output_stash = \
                pjwrite.prep_data_df(self.gui_modeling.pgm_obj)
            
            # collecting scalars #
            _, scalar_value_df_dict, scalar_repl_summary_df = \
                scalar_output_stash

            # collecting trees #
            tree_value_df_dict, tree_ann_value_df_dict, tree_rec_value_df_dict, \
                tree_rec_ann_value_df_dict, tree_summary_df_dict, tree_repl_summary_df_dict, \
                tree_living_nd_states_str_dict, tree_living_nd_states_str_nexus_dict, \
                tree_internal_nd_states_str_dict = tree_output_stash
                    
            # adding "program" column to
            # PJ's pandas DataFrame's stored inside dicts
            for repl_size in scalar_value_df_dict:
                scalar_value_df_dict[selected_node_repl_size].loc[:, "program"] = "PJ"

            for tr_nd_name in tree_summary_df_dict:
                    tree_summary_df_dict[tr_nd_name].loc[:, "program"] = "PJ"
                    
            for tr_w_repl_nd_name in tree_repl_summary_df_dict:
                tree_repl_summary_df_dict[tr_w_repl_nd_name].loc[:, "program"] = "PJ"

            # scalar was selected #
            if isinstance(selected_node_pgm.value[0], (int, float, np.float64)):
                self.ui.ui_pages.summary_stats_list.clear()
                
                if not self.is_avg_repl_check:
                    self.pj_comparison_df = scalar_value_df_dict[selected_node_repl_size]

                else:
                    self.ui.ui_pages.summary_stats_list.addItems(
                        ["average", "std. dev."]
                    )
                    self.pj_comparison_df = scalar_repl_summary_df
                

            # tree was selected #
            elif isinstance(selected_node_pgm.value[0], pjdt.AnnotatedTree):
                self.pj_comparison_df = tree_summary_df_dict[selected_node_name]
                self.ui.ui_pages.summary_stats_list.clear()
                self.ui.ui_pages.summary_stats_list.addItems(
                    [tree_stat for tree_stat in self.pj_comparison_df.keys() \
                        if tree_stat not in \
                            ("program", "sample", "replicate")]
                )

        # no sampled nodes
        else:
            print("If there are no sampled nodes in PJ model, we cannot compare against PJ. Later write an Exception for this")
            pass


    ##################
    # Functions for  #
    # reading events #
    ##################
    
    def read_execute_script(self):

        # reset everything #
        self.clean_disable_everything()

        # read file path #
        script_in_fp, filter = QFileDialog.getOpenFileName(parent=self, caption="Read script", dir=".", filter="*.pj")

        if script_in_fp:
            cmd_line_list = pjread.read_text_file(script_in_fp)

            # (side-effect: gui_modeling stores cmd hist) #
            self.gui_modeling.parse_cmd_update_pgm(cmd_line_list, self)

            # update GUI cmd history #
            self.refresh_cmd_history()

            # add node names to...
            # pgm page node list
            self.ui.ui_pages.node_list.addItems(
                [node_name for node_name in \
                    self.gui_modeling.pgm_obj.node_name_val_dict]
            )

            # compare page node list
            self.ui.ui_pages.compare_node_list.addItems(
                [nd.node_name for nd in \
                    self.gui_modeling.pgm_obj.get_sorted_node_pgm_list() if \
                        nd.is_sampled]
            )

            # coverage page node list
            # TODO


    def read_compare_csv(self):

        # read file path #
        compare_csv_fp, filter = QFileDialog.getOpenFileName(parent=self, caption="Read compare .csv", dir=".", filter="*.csv")

        if compare_csv_fp:
            if not os.path.isfile(compare_csv_fp):
                print("Could not find " + compare_csv_fp) # TODO: add exception here later

            try:
                self.other_comparison_df = pjread.read_csv_into_dataframe(compare_csv_fp)
                other_comparison_df_str = tabulate(self.other_comparison_df, self.other_comparison_df.head(), tablefmt="plain", showindex=False).lstrip()
                self.ui.ui_pages.compare_csv_textbox.setText(other_comparison_df_str)
                
            except:
                print("Could not load .csv file into pandas DataFrame. Later write an Exception for this")
            
        else:
            pass # event canceled by user

    
    #################
    # Functions for #
    # drawing       #
    #################

    def draw_node_pgm(self, axes, node_pgm, sample_idx=None, repl_idx=0, repl_size=1):
        return node_pgm.plot_node(axes, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)


    def draw_violin(self):
        node_name: str = ""
        thing_to_compare: str = ""
        compare_node_list_empty = self.ui.ui_pages.compare_node_list.count() == 0
        summary_stat_list_empty = self.ui.ui_pages.summary_stats_list.count() == 0

        # before draw, we reset violins
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()
        
        # first get selected node's name #
        node_name = \
            self.ui.ui_pages.compare_node_list.currentItem().text()

        # reading node information #
        node_pgm, sample_size, repl_size = \
            self.selected_node_read(node_name)

        if not self.pj_comparison_df.empty and \
            not self.other_comparison_df.empty:

            # scalar
            if isinstance(node_pgm.value[0], (int, float, np.float64)):
                if not self.is_avg_repl_check:
                    thing_to_compare = node_name
                
                # averaging over replicates
                else:
                    if not summary_stat_list_empty:
                        if self.ui.ui_pages.summary_stats_list.currentItem():
                            thing_to_compare = self.ui.ui_pages.summary_stats_list.currentItem().text()

                        # summary stats were not picked, doing nothing
                        else:
                            pass
            
            # custom class object (e.g., AnnotatedTree)
            else:
                if not summary_stat_list_empty:
                    if self.ui.ui_pages.summary_stats_list.currentItem():
                        thing_to_compare = self.ui.ui_pages.summary_stats_list.currentItem().text()
                
                    # summary stats were not picked, doing nothing
                    else:
                        pass

            # debugging
            # print(tabulate(self.pj_comparison_df, self.pj_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())
            # print(tabulate(self.other_comparison_df, self.other_comparison_df.head(), tablefmt="pretty", showindex=False).lstrip())

            joint_dataframe = pjorg.join_dataframes(
                self.pj_comparison_df, self.other_comparison_df,
                value_to_compare=thing_to_compare,
                summaries_avg_over_repl=False)

            # something went wrong,
            # we clear comparison figure
            if joint_dataframe.empty:                    
                print("joint_dataframe was empty")
                pass

            else:
                # initializing plot
                self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes(disabled_yticks=False)
                self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()
                fig_obj = self.ui.ui_pages.compare_page_matplotlib_widget.fig
                fig_axes = self.ui.ui_pages.compare_page_matplotlib_widget.axes
                
                pjdraw.plot_violins(fig_obj, fig_axes,
                    joint_dataframe, "program", thing_to_compare,
                    xlab="Program", ylab=thing_to_compare)


    #####################
    # Events related to #
    # initializing,     #
    # refreshing,       #
    # cleaning,         #
    # checking          #
    #####################

    def init_and_refresh_radio_spin(self, node_pgm, sample_size):

        ###################
        # Stochastic node #
        ###################

        if node_pgm.is_sampled:
            # tree #
            if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                # radio #
                # we always look one tree at a time
                self.ui.ui_pages.all_samples_radio.setCheckable(False)
                self.ui.ui_pages.all_samples_radio.setDisabled(True)
                
                self.ui.ui_pages.one_sample_radio.setEnabled(True)
                self.ui.ui_pages.one_sample_radio.setCheckable(True)
                self.ui.ui_pages.one_sample_radio.setChecked(True)

                # spin #
                self.ui.ui_pages.sample_idx_spin.setEnabled(True)
                self.ui.ui_pages.sample_idx_spin.setMinimum(1)
                self.ui.ui_pages.repl_idx_spin.setEnabled(True)
                self.ui.ui_pages.repl_idx_spin.setMinimum(1)

            # non-tree
            else:
                # TODO: coalesce the if and the else below
                
                # if it's the first click selecting node
                if not self.ui.ui_pages.all_samples_radio.isEnabled() and \
                    not self.ui.ui_pages.one_sample_radio.isEnabled():

                    # radio #
                    # one sample is the default
                    self.ui.ui_pages.one_sample_radio.setCheckable(True)
                    self.ui.ui_pages.one_sample_radio.setChecked(True)

                    # spin #
                    self.ui.ui_pages.sample_idx_spin.setMinimum(0)
                    self.ui.ui_pages.sample_idx_spin.setValue(0)
                    self.ui.ui_pages.sample_idx_spin.setDisabled(True)
                    
                    self.ui.ui_pages.repl_idx_spin.setMinimum(0)
                    self.ui.ui_pages.repl_idx_spin.setValue(0)
                    self.ui.ui_pages.repl_idx_spin.setDisabled(True)

                else:
                    self.ui.ui_pages.all_samples_radio.setEnabled(True)
                    self.ui.ui_pages.one_sample_radio.setEnabled(True)
                    
                    # if looking at all samples, we
                    # don't circle through repls
                    if self.ui.ui_pages.all_samples_radio.isChecked():
                        self.ui.ui_pages.repl_idx_spin.setMinimum(0)
                        self.ui.ui_pages.repl_idx_spin.setValue(0)
                        self.ui.ui_pages.repl_idx_spin.setDisabled(True)
                        
                        self.ui.ui_pages.sample_idx_spin.setMinimum(0)
                        self.ui.ui_pages.sample_idx_spin.setValue(0)
                        self.ui.ui_pages.sample_idx_spin.setDisabled(True)

                    # if looking at one sample at a time
                    else:
                        self.ui.ui_pages.repl_idx_spin.setEnabled(True)
                        self.ui.ui_pages.repl_idx_spin.setMinimum(1)
                        self.ui.ui_pages.sample_idx_spin.setMaximum(100)
                        self.ui.ui_pages.repl_idx_spin.setValue(1)

                        self.ui.ui_pages.sample_idx_spin.setEnabled(True)
                        self.ui.ui_pages.sample_idx_spin.setMinimum(1)
                        self.ui.ui_pages.sample_idx_spin.setMaximum(sample_size)
                        self.ui.ui_pages.sample_idx_spin.setValue(1)


        #################
        # Constant node #
        #################

        # cannot circle through replicates
        # because no 2D-nesting when 
        # assigning constants (and no
        # repls in deterministic nodes)
        else:
            # NOTE: we need to both disable
            # and uncheck so that the radio
            # button is totally reset and
            # turned off

            # radio #
            self.ui.ui_pages.one_sample_radio.setCheckable(False)
            self.ui.ui_pages.one_sample_radio.setDisabled(True)
            
            # if deterministic
            if sample_size == 0:
                self.ui.ui_pages.all_samples_radio.setCheckable(False)
                self.ui.ui_pages.all_samples_radio.setDisabled(True)
            
            else: 
                self.ui.ui_pages.all_samples_radio.setEnabled(True)
                self.ui.ui_pages.all_samples_radio.setCheckable(True)
                self.ui.ui_pages.all_samples_radio.setChecked(True)

            # spin #
            self.ui.ui_pages.sample_idx_spin.setMinimum(0)
            self.ui.ui_pages.sample_idx_spin.setValue(0)
            self.ui.ui_pages.sample_idx_spin.setDisabled(True)
            self.ui.ui_pages.repl_idx_spin.setMinimum(0)
            self.ui.ui_pages.repl_idx_spin.setValue(0)
            self.ui.ui_pages.repl_idx_spin.setDisabled(True)


    def refresh_selected_node_display_plot(self):
        # if nodes have been created and selected #
        if self.ui.ui_pages.node_list.currentItem() != None:
            self.do_selected_node_pgm_page()


    def refresh_cmd_history(self, user_reset=False):
        if user_reset:
            self.gui_modeling.cmd_log_list.append("\n## RESET: Model is being cleared at this point\n")

        cmd_hist_str = self.gui_modeling.cmd_log()
        self.ui.ui_pages.cmd_log_textbox.setText(cmd_hist_str)


    def clean_disable_everything(self):
        # pgm-related objects are reset
        self.gui_modeling.clear()

        # update GUI cmd history
        self.refresh_cmd_history(user_reset=True)

        # remove all node names from list
        self.ui.ui_pages.node_list.clear()
        self.ui.ui_pages.compare_node_list.clear()
        self.ui.ui_pages.summary_stats_list.clear()

        # reset text on display and summary
        # panels
        self.ui.ui_pages.values_content.clear()
        self.ui.ui_pages.summary_content.clear()

        # resetting and disabling buttons
        # radio #
        self.ui.ui_pages.one_sample_radio.setDisabled(True)
        self.ui.ui_pages.one_sample_radio.setCheckable(False)
        self.ui.ui_pages.one_sample_radio.setChecked(False)

        self.ui.ui_pages.all_samples_radio.setDisabled(True)
        self.ui.ui_pages.all_samples_radio.setCheckable(False)
        self.ui.ui_pages.all_samples_radio.setChecked(False)

        # spin #
        self.ui.ui_pages.sample_idx_spin.setMinimum(0)
        self.ui.ui_pages.sample_idx_spin.setValue(0)
        self.ui.ui_pages.sample_idx_spin.setDisabled(True)
        self.ui.ui_pages.repl_idx_spin.setMinimum(0)
        self.ui.ui_pages.repl_idx_spin.setValue(0)
        self.ui.ui_pages.repl_idx_spin.setDisabled(True)

        # reset plot(s)
        self.ui.ui_pages.pgm_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.pgm_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.pgm_page_matplotlib_widget.fig.canvas.draw()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.clf()
        self.ui.ui_pages.compare_page_matplotlib_widget.initialize_axes()
        self.ui.ui_pages.compare_page_matplotlib_widget.fig.canvas.draw()


    # side-effect: sets is_avg_repl_check
    # member as True or False
    def avg_repl_check(self):
        chk_button = self.sender()
        self.is_avg_repl_check = chk_button.isChecked()
        
        # clearing summary stats
        # (should only be there
        # if averaging over repls)
        if not self.is_avg_repl_check:
            self.ui.ui_pages.summary_stats_list.clear()

        # clicking check box should de-select node
        # being compared to force user to select
        self.ui.ui_pages.compare_node_list.clearSelection()


    def quit_app_button_clicked(self):
        sys.exit()


    def print_about(self):
        print("oi")


def call_gui():
    gui_app = QApplication(sys.argv)
    window = GUIMainWindow()

    timer = QTimer()
    timer.timeout.connect(lambda: None)
    timer.start(100)

    sys.exit(gui_app.exec())


# call GUI # 
if __name__ == "__main__":
    call_gui()
