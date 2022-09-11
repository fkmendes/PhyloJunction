import sys
import os
import re
import typing as ty

from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton
from PySide6.QtCore import QPropertyAnimation, QEasingCurve, QTimer

# pj imports #
from phylojunction.interface.pysidegui.content_main_window \
    import ContentGUIMainWindow
from phylojunction.pgm.pgm import ProbabilisticGraphicalModel
import phylojunction.interface.cmdbox.cmd_parse as cmdp
import phylojunction.readwrite.pj_read as pjread
import phylojunction.data.tree as pjdt


class GUIMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # initialize all model-related members #
        self.gui_modeling = GUIModeling()

        self.setWindowTitle("PhyloJunction")

        # setup #
        self.ui = ContentGUIMainWindow()
        # passing myself as argument to auxiliary setup class
        self.ui.setup_ui(self)

        # toggle button #
        self.ui.menu_button.clicked.connect(self.toggle_button_animation)

        # menu button #
        self.ui.pgm_button.clicked.connect(self.show_pgm_page)

        # pgm page button #
        self.ui.cmd_log_button.clicked.connect(self.show_cmd_log_page)

        # cmd line enter #
        self.ui.ui_pages.cmd_prompt.returnPressed.connect(self.parse_cmd_update_gui)

        # node list #
        self.ui.ui_pages.node_list.itemClicked.connect(self.do_selected_node)

        # radio button update #
        self.ui.ui_pages.one_sample_radio.toggled.connect(lambda: self.refresh_selected_node_display_plot(True))
        self.ui.ui_pages.all_samples_radio.toggled.connect(lambda: self.refresh_selected_node_display_plot(False))

        # spin boxes update #
        self.ui.ui_pages.sample_idx_spin.valueChanged.connect(self.refresh_selected_node_display_plot)
        self.ui.ui_pages.repl_idx_spin.valueChanged.connect(self.refresh_selected_node_display_plot)

        # show #
        self.show()


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
        # ui_pages is a member of self.ui, which is an instance
        # of ContentGUIMainWindow
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.pgm_page)
        self.ui.pgm_button.set_active(True)


    def show_cmd_log_page(self):
        self.reset_selection()
        self.ui.pages.setCurrentWidget(self.ui.ui_pages.cmd_log_page)
        self.ui.cmd_log_button.set_active(True)


    ###################################
    # Various event-related functions #
    ###################################

    def parse_cmd_update_gui(self):
        cmd_line = self.ui.ui_pages.cmd_prompt.text()  # grab input
        # in case the user enters multiple lines, also ignores empty lines
        cmd_line_list = [cl for cl in re.split("\n", cmd_line) if cl]
        self.gui_modeling.parse_cmd_update_pgm(cmd_line_list)
        self.ui.ui_pages.cmd_prompt.clear()  # clear

        # now update ListWidget #
        pgm_node_list = [node_name for \
                node_name in self.gui_modeling.pgm_obj.node_name_val_dict]
        self.ui.ui_pages.node_list.clear()  # clear first
        self.ui.ui_pages.node_list.addItems(pgm_node_list)


    def selected_node_read(self, node_name):
        node_pgm = self.gui_modeling.pgm_obj.get_node_pgm_by_name(node_name)
        # this is n_sim inside sampling distribution classes
        sample_size = len(node_pgm)
        repl_size = node_pgm.repl_size

        return node_pgm, sample_size, repl_size


    def draw_node_pgm(self, axes, node_pgm, sample_idx=None, repl_idx=0, repl_size=1):
        return node_pgm.plot_node(axes, sample_idx=sample_idx, repl_idx=repl_idx, repl_size=repl_size)

        
    def selected_node_display(self, node_pgm, do_all_samples, sample_idx=None, repl_idx=0, repl_size=1):
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
            display_node_pgm_value_str = self.gui_modeling.pgm_obj.get_display_str_by_name(node_pgm.node_name)
            # getting all values
            display_node_pgm_stat_str = node_pgm.get_node_stats_str(0, len(node_pgm.value), repl_idx) # summary stats
        
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


    def do_selected_node(self):
        """
        Display selected node's string representation and
        plot it on canvas if possible
        """
        
        # first get selected node's name
        node_name = self.ui.ui_pages.node_list.currentItem().text()
        
        # grab pgm_page's figure and axes #
        fig_obj = self.ui.ui_pages.pgm_page_matplotlib_widget.fig
        fig_axes = self.ui.ui_pages.pgm_page_matplotlib_widget.axes

        node_pgm, sample_size, repl_size = \
            self.selected_node_read(node_name)

        # updating spin boxes #
        self.ui.ui_pages.sample_idx_spin.setMaximum(sample_size - 1)
        self.ui.ui_pages.repl_idx_spin.setMaximum(repl_size - 1)

        # doing all samples?
        do_all_samples = \
            self.ui.ui_pages.all_samples_radio.isChecked()
        
        # if node is a tree, we always look one tree at
        # a time
        if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
            do_all_samples = False
        
        # print("do_all_samples = " + str(do_all_samples))

        # TODO:
        # sample_idx = 0  # = int(wdw["-ITH-SAMPLE-"].get()) - 1 # (offset)
        # repl_idx = 0  # = int(wdw["-ITH-REPL-"].get()) - 1 # (offset)
        sample_idx = self.ui.ui_pages.sample_idx_spin.value()
        repl_idx = self.ui.ui_pages.repl_idx_spin.value()

        self.selected_node_display(node_pgm,
                                   do_all_samples,
                                   sample_idx=sample_idx,
                                   repl_idx=0,
                                   repl_size=1)

        self.selected_node_plot(fig_obj,
                                fig_axes,
                                node_pgm,
                                do_all_samples,
                                sample_idx=sample_idx,
                                repl_idx=repl_idx,
                                repl_size=repl_size)      


    def refresh_selected_node_display_plot(self, one_sample: bool):
        if self.ui.ui_pages.node_list.currentItem() != None:
            node_pgm, sample_size, repl_size = \
                self.selected_node_read(
                    self.ui.ui_pages.node_list.currentItem().text())

            # if one sample at a time, spin element
            # works (only for sampled random variables)
            if one_sample:
                self.ui.ui_pages.repl_idx_spin.setMaximum(sample_size)
                self.ui.ui_pages.sample_idx_spin.setEnabled(True)

                # if we are looking at trees, we can cycle through replicates
                if isinstance(node_pgm.value[0], pjdt.AnnotatedTree):
                    self.ui.ui_pages.repl_idx_spin.setMaximum(repl_size)
                    self.ui.ui_pages.repl_idx_spin.setEnabled(True)

                # otherwise, all replicates will be visualized
                # as histogram (no cycling allowed)
                else:
                    self.ui.ui_pages.repl_idx_spin.setDisabled(True)

            # if we're looking at all samples, all samples
            # and replicates shown in histogram
            #
            # and if tree, only one sample at a time is allowed
            else:
                self.ui.ui_pages.repl_idx_spin.setValue(0)
                self.ui.ui_pages.repl_idx_spin.setMaximum(0)
                self.ui.ui_pages.repl_idx_spin.setDisabled(True)
                self.ui.ui_pages.sample_idx_spin.setValue(0)
                self.ui.ui_pages.sample_idx_spin.setMaximum(0)
                self.ui.ui_pages.sample_idx_spin.setDisabled(True)

            # now refresh selected node display and plot #
            self.do_selected_node()


class GUIModeling():
    def __init__(self):
        self.pgm_obj: ProbabilisticGraphicalModel = \
            ProbabilisticGraphicalModel()
        self.cmd_log_list: ty.List[str] = []
        self.cmd_log: str = ""


    def parse_cmd_update_pgm(self, cmd_line_list):
        valid_cmd_line = None
        
        for line in cmd_line_list:
            # removing whitespaces from left and right
            line = line.strip()
            
            try:
                valid_cmd_line = cmdp.cmdline2pgm(self.pgm_obj, line)
        
            except Exception as e:
                # if not event == "Simulate":
                if e.__context__:
                    # __context__ catches innermost exception 
                    # update warnings page later, with e.__context__
                    print(e.__context__)
                else:
                    print(e)  # update warnings page with e

            if valid_cmd_line:
                self.cmd_log_list.append(valid_cmd_line)
                self.cmd_log = "\n".join(self.cmd_log_list) 


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
