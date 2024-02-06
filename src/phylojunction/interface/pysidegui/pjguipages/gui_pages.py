# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'src/phylojunction/interface/pysidegui/pjguipages/pjgui_pages.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PySide6 import QtCore, QtGui, QtWidgets


class Ui_PJGUIPages(object):
    def setupUi(self, PJGUIPages):
        PJGUIPages.setObjectName("PJGUIPages")
        PJGUIPages.resize(980, 700)
        PJGUIPages.setMinimumSize(QtCore.QSize(960, 650))
        PJGUIPages.setMaximumSize(QtCore.QSize(16777215, 700))
        self.settings_page = QtWidgets.QWidget()
        self.settings_page.setObjectName("settings_page")
        self.gridLayoutWidget_4 = QtWidgets.QWidget(self.settings_page)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(10, 10, 961, 681))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.settings_page_grid_layout = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.settings_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.settings_page_grid_layout.setObjectName("settings_page_grid_layout")
        self.filename_prefix_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.filename_prefix_label.setMinimumSize(QtCore.QSize(110, 24))
        self.filename_prefix_label.setMaximumSize(QtCore.QSize(110, 24))
        self.filename_prefix_label.setObjectName("filename_prefix_label")
        self.settings_page_grid_layout.addWidget(self.filename_prefix_label, 1, 0, 1, 1)
        self.random_seed_prefix_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.random_seed_prefix_label.setMinimumSize(QtCore.QSize(110, 24))
        self.random_seed_prefix_label.setMaximumSize(QtCore.QSize(110, 24))
        self.random_seed_prefix_label.setObjectName("random_seed_prefix_label")
        self.settings_page_grid_layout.addWidget(self.random_seed_prefix_label, 2, 0, 1, 1)
        self.settings_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.settings_label.setMinimumSize(QtCore.QSize(0, 24))
        self.settings_label.setMaximumSize(QtCore.QSize(16777215, 24))
        font = QtGui.QFont()
        font.setBold(True)
        self.settings_label.setFont(font)
        self.settings_label.setObjectName("settings_label")
        self.settings_page_grid_layout.addWidget(self.settings_label, 0, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.settings_page_grid_layout.addItem(spacerItem, 3, 0, 1, 1)
        self.filename_prefix_textbox = QtWidgets.QTextEdit(self.gridLayoutWidget_4)
        self.filename_prefix_textbox.setMinimumSize(QtCore.QSize(0, 24))
        self.filename_prefix_textbox.setMaximumSize(QtCore.QSize(16777215, 24))
        self.filename_prefix_textbox.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.filename_prefix_textbox.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.filename_prefix_textbox.setObjectName("filename_prefix_textbox")
        self.settings_page_grid_layout.addWidget(self.filename_prefix_textbox, 1, 1, 1, 1)
        self.random_seed_prefix_textbox = QtWidgets.QTextEdit(self.gridLayoutWidget_4)
        self.random_seed_prefix_textbox.setMinimumSize(QtCore.QSize(0, 24))
        self.random_seed_prefix_textbox.setMaximumSize(QtCore.QSize(16777215, 24))
        self.random_seed_prefix_textbox.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.random_seed_prefix_textbox.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.random_seed_prefix_textbox.setObjectName("random_seed_prefix_textbox")
        self.settings_page_grid_layout.addWidget(self.random_seed_prefix_textbox, 2, 1, 1, 1)
        self.line = QtWidgets.QFrame(self.gridLayoutWidget_4)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.settings_page_grid_layout.addWidget(self.line, 0, 1, 1, 1)
        PJGUIPages.addWidget(self.settings_page)
        self.pgm_page = QtWidgets.QWidget()
        self.pgm_page.setMinimumSize(QtCore.QSize(980, 700))
        self.pgm_page.setMaximumSize(QtCore.QSize(16777215, 700))
        self.pgm_page.setObjectName("pgm_page")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.pgm_page)
        self.verticalLayout.setObjectName("verticalLayout")
        self.pgm_page_frame = QtWidgets.QFrame(self.pgm_page)
        self.pgm_page_frame.setMinimumSize(QtCore.QSize(980, 700))
        self.pgm_page_frame.setMaximumSize(QtCore.QSize(16777215, 700))
        self.pgm_page_frame.setStyleSheet("background-color: white;\n"
"border: 0;")
        self.pgm_page_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.pgm_page_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.pgm_page_frame.setObjectName("pgm_page_frame")
        self.gridLayoutWidget = QtWidgets.QWidget(self.pgm_page_frame)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(0, 170, 954, 451))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.pgm_fig_node_list_grid_layout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.pgm_fig_node_list_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.pgm_fig_node_list_grid_layout.setObjectName("pgm_fig_node_list_grid_layout")
        self.model_label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.model_label.setMinimumSize(QtCore.QSize(200, 25))
        self.model_label.setMaximumSize(QtCore.QSize(200, 20))
        self.model_label.setStyleSheet("font: 14pt \"Ubuntu\";\n"
"color: black;\n"
"")
        self.model_label.setAlignment(QtCore.Qt.AlignCenter)
        self.model_label.setObjectName("model_label")
        self.pgm_fig_node_list_grid_layout.addWidget(self.model_label, 0, 1, 1, 1)
        self.pgm_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget)
        self.pgm_page_matplotlib_widget.setMinimumSize(QtCore.QSize(740, 400))
        self.pgm_page_matplotlib_widget.setMaximumSize(QtCore.QSize(16777215, 400))
        self.pgm_page_matplotlib_widget.setObjectName("pgm_page_matplotlib_widget")
        self.pgm_fig_node_list_grid_layout.addWidget(self.pgm_page_matplotlib_widget, 1, 0, 1, 1)
        self.node_list_vert_layout = QtWidgets.QVBoxLayout()
        self.node_list_vert_layout.setObjectName("node_list_vert_layout")
        self.node_list = QtWidgets.QListWidget(self.gridLayoutWidget)
        self.node_list.setMinimumSize(QtCore.QSize(195, 250))
        self.node_list.setMaximumSize(QtCore.QSize(195, 380))
        self.node_list.setStyleSheet("QListWidget {\n"
"    background-color: #f1f3f5;\n"
"    color: #495057;\n"
"    border-radius: 5px;\n"
"    border: 0;\n"
"}\n"
"QListWidget::item:selected{\n"
"    background-color: #495057;\n"
"    color: white;\n"
"    border: 0;\n"
"}\n"
"")
        self.node_list.setObjectName("node_list")
        self.node_list_vert_layout.addWidget(self.node_list)
        self.clear_model = PJClearPGMQPushButton(self.gridLayoutWidget)
        self.clear_model.setEnabled(True)
        self.clear_model.setMinimumSize(QtCore.QSize(195, 25))
        self.clear_model.setMaximumSize(QtCore.QSize(195, 16777215))
        self.clear_model.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.clear_model.setStyleSheet("QPushButton:hover {\n"
"    color: #ec4a8a;\n"
"}")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("src/phylojunction/interface/pysidegui/pjguipages/../images/icons/icon_clear.svg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.clear_model.setIcon(icon)
        self.clear_model.setIconSize(QtCore.QSize(20, 20))
        self.clear_model.setObjectName("clear_model")
        self.node_list_vert_layout.addWidget(self.clear_model)
        self.pgm_fig_node_list_grid_layout.addLayout(self.node_list_vert_layout, 1, 1, 1, 1)
        self.radio_spin_hor_layout = QtWidgets.QHBoxLayout()
        self.radio_spin_hor_layout.setObjectName("radio_spin_hor_layout")
        self.one_sample_radio = QtWidgets.QRadioButton(self.gridLayoutWidget)
        self.one_sample_radio.setEnabled(False)
        self.one_sample_radio.setStyleSheet("color: black;\n"
"")
        self.one_sample_radio.setCheckable(False)
        self.one_sample_radio.setChecked(False)
        self.one_sample_radio.setObjectName("one_sample_radio")
        self.radio_spin_hor_layout.addWidget(self.one_sample_radio)
        self.all_samples_radio = QtWidgets.QRadioButton(self.gridLayoutWidget)
        self.all_samples_radio.setEnabled(False)
        self.all_samples_radio.setStyleSheet("color: black;")
        self.all_samples_radio.setCheckable(False)
        self.all_samples_radio.setAutoExclusive(True)
        self.all_samples_radio.setObjectName("all_samples_radio")
        self.radio_spin_hor_layout.addWidget(self.all_samples_radio)
        spacerItem1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Minimum)
        self.radio_spin_hor_layout.addItem(spacerItem1)
        self.sample_idx_spin = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.sample_idx_spin.setEnabled(False)
        self.sample_idx_spin.setStyleSheet("color: black;")
        self.sample_idx_spin.setAccelerated(True)
        self.sample_idx_spin.setObjectName("sample_idx_spin")
        self.radio_spin_hor_layout.addWidget(self.sample_idx_spin)
        spacerItem2 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Minimum)
        self.radio_spin_hor_layout.addItem(spacerItem2)
        self.repl_idx_spin = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.repl_idx_spin.setEnabled(False)
        self.repl_idx_spin.setStyleSheet("color: black;\n"
"")
        self.repl_idx_spin.setReadOnly(False)
        self.repl_idx_spin.setAccelerated(True)
        self.repl_idx_spin.setObjectName("repl_idx_spin")
        self.radio_spin_hor_layout.addWidget(self.repl_idx_spin)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.radio_spin_hor_layout.addItem(spacerItem3)
        self.save_pgm_node_plot = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.save_pgm_node_plot.setEnabled(True)
        self.save_pgm_node_plot.setMinimumSize(QtCore.QSize(130, 24))
        self.save_pgm_node_plot.setMaximumSize(QtCore.QSize(130, 24))
        self.save_pgm_node_plot.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.save_pgm_node_plot.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.save_pgm_node_plot.setObjectName("save_pgm_node_plot")
        self.radio_spin_hor_layout.addWidget(self.save_pgm_node_plot)
        self.pgm_fig_node_list_grid_layout.addLayout(self.radio_spin_hor_layout, 0, 0, 1, 1)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.pgm_page_frame)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 620, 951, 65))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.cmd_prompt_vert_layout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.cmd_prompt_vert_layout.setContentsMargins(0, 0, 0, 0)
        self.cmd_prompt_vert_layout.setObjectName("cmd_prompt_vert_layout")
        self.cmd_prompt_label = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.cmd_prompt_label.setMinimumSize(QtCore.QSize(940, 20))
        self.cmd_prompt_label.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.cmd_prompt_label.setStyleSheet("font: 14pt \"Ubuntu\";\n"
"color: black;")
        self.cmd_prompt_label.setAlignment(QtCore.Qt.AlignBottom|QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft)
        self.cmd_prompt_label.setObjectName("cmd_prompt_label")
        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt_label)
        self.cmd_prompt = QtWidgets.QLineEdit(self.verticalLayoutWidget)
        self.cmd_prompt.setMinimumSize(QtCore.QSize(940, 35))
        self.cmd_prompt.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setFamily("Courier")
        font.setPointSize(14)
        font.setBold(True)
        font.setItalic(False)
        self.cmd_prompt.setFont(font)
        self.cmd_prompt.setToolTip("")
        self.cmd_prompt.setStyleSheet("background-color: black;\n"
"font: 14pt \"Courier\";\n"
"font-weight: bold;\n"
"color: white;\n"
"border-radius: 5px;")
        self.cmd_prompt.setCursorPosition(0)
        self.cmd_prompt.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.cmd_prompt.setDragEnabled(True)
        self.cmd_prompt.setPlaceholderText("")
        self.cmd_prompt.setObjectName("cmd_prompt")
        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt)
        self.node_content_tabs = QtWidgets.QTabWidget(self.pgm_page_frame)
        self.node_content_tabs.setGeometry(QtCore.QRect(0, 0, 951, 171))
        self.node_content_tabs.setMinimumSize(QtCore.QSize(940, 0))
        self.node_content_tabs.setStyleSheet("color: black;\n"
"font: 14pt ;\n"
"")
        self.node_content_tabs.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.node_content_tabs.setIconSize(QtCore.QSize(16, 16))
        self.node_content_tabs.setObjectName("node_content_tabs")
        self.values_tab = QtWidgets.QWidget()
        self.values_tab.setObjectName("values_tab")
        self.values_content = QtWidgets.QTextEdit(self.values_tab)
        self.values_content.setGeometry(QtCore.QRect(0, 0, 951, 141))
        self.values_content.setMinimumSize(QtCore.QSize(940, 0))
        self.values_content.setAcceptDrops(False)
        self.values_content.setStyleSheet("border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.values_content.setReadOnly(True)
        self.values_content.setObjectName("values_content")
        self.node_content_tabs.addTab(self.values_tab, "")
        self.summary_tab = QtWidgets.QWidget()
        self.summary_tab.setObjectName("summary_tab")
        self.summary_content = QtWidgets.QTextEdit(self.summary_tab)
        self.summary_content.setGeometry(QtCore.QRect(0, 0, 951, 141))
        self.summary_content.setMinimumSize(QtCore.QSize(940, 0))
        self.summary_content.setStyleSheet("border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.summary_content.setObjectName("summary_content")
        self.node_content_tabs.addTab(self.summary_tab, "")
        self.verticalLayout.addWidget(self.pgm_page_frame)
        PJGUIPages.addWidget(self.pgm_page)
        self.compare_page = QtWidgets.QWidget()
        self.compare_page.setEnabled(True)
        self.compare_page.setObjectName("compare_page")
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.compare_page)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(10, 10, 961, 681))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.compare_page_grid_layout = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.compare_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.compare_page_grid_layout.setObjectName("compare_page_grid_layout")
        self.node_stat_vert_layout = QtWidgets.QVBoxLayout()
        self.node_stat_vert_layout.setObjectName("node_stat_vert_layout")
        self.compare_node_frame = QtWidgets.QFrame(self.gridLayoutWidget_2)
        self.compare_node_frame.setMinimumSize(QtCore.QSize(195, 350))
        self.compare_node_frame.setMaximumSize(QtCore.QSize(195, 350))
        self.compare_node_frame.setStyleSheet("border-radius: 5px;\n"
"background-color: #e6f4f4;\n"
"border: 0;")
        self.compare_node_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.compare_node_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.compare_node_frame.setObjectName("compare_node_frame")
        self.compare_node_list = QtWidgets.QListWidget(self.compare_node_frame)
        self.compare_node_list.setGeometry(QtCore.QRect(10, 30, 172, 280))
        self.compare_node_list.setMinimumSize(QtCore.QSize(172, 280))
        self.compare_node_list.setMaximumSize(QtCore.QSize(172, 280))
        self.compare_node_list.setStyleSheet("QListWidget {\n"
"    background-color: #f2f9f9;\n"
"    color: #495057;\n"
"    border: 0;\n"
"    border-radius: 0px;\n"
"}\n"
"QListWidget::item:selected{\n"
"    background-color: #495057;\n"
"    color: white;\n"
"    border: 0;\n"
"}")
        self.compare_node_list.setObjectName("compare_node_list")
        self.compare_node_label = QtWidgets.QLabel(self.compare_node_frame)
        self.compare_node_label.setGeometry(QtCore.QRect(20, 10, 155, 16))
        self.compare_node_label.setMinimumSize(QtCore.QSize(155, 16))
        self.compare_node_label.setMaximumSize(QtCore.QSize(155, 16))
        self.compare_node_label.setAlignment(QtCore.Qt.AlignCenter)
        self.compare_node_label.setObjectName("compare_node_label")
        self.avg_replicate_check_button = QtWidgets.QCheckBox(self.compare_node_frame)
        self.avg_replicate_check_button.setGeometry(QtCore.QRect(40, 320, 155, 20))
        self.avg_replicate_check_button.setMinimumSize(QtCore.QSize(155, 20))
        self.avg_replicate_check_button.setMaximumSize(QtCore.QSize(155, 20))
        font = QtGui.QFont()
        font.setPointSize(9)
        font.setItalic(False)
        self.avg_replicate_check_button.setFont(font)
        self.avg_replicate_check_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.avg_replicate_check_button.setObjectName("avg_replicate_check_button")
        self.node_stat_vert_layout.addWidget(self.compare_node_frame)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.node_stat_vert_layout.addItem(spacerItem4)
        self.compare_stats_label = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.compare_stats_label.setMinimumSize(QtCore.QSize(0, 16))
        self.compare_stats_label.setMaximumSize(QtCore.QSize(16777215, 16))
        self.compare_stats_label.setObjectName("compare_stats_label")
        self.node_stat_vert_layout.addWidget(self.compare_stats_label, 0, QtCore.Qt.AlignHCenter)
        self.summary_stats_list = QtWidgets.QListWidget(self.gridLayoutWidget_2)
        self.summary_stats_list.setMinimumSize(QtCore.QSize(195, 250))
        self.summary_stats_list.setMaximumSize(QtCore.QSize(195, 250))
        self.summary_stats_list.setStyleSheet("QListWidget{\n"
"    background-color: #f7f7f7;\n"
"    border: 2px solid lightgray;\n"
"}\n"
"QListWidget::item:selected{\n"
"    background-color: #495057;\n"
"    color: white;\n"
"    border: 0;\n"
"}")
        self.summary_stats_list.setObjectName("summary_stats_list")
        self.node_stat_vert_layout.addWidget(self.summary_stats_list, 0, QtCore.Qt.AlignHCenter)
        self.compare_page_grid_layout.addLayout(self.node_stat_vert_layout, 0, 0, 1, 1)
        self.csv_violinplot_vert_layout = QtWidgets.QVBoxLayout()
        self.csv_violinplot_vert_layout.setObjectName("csv_violinplot_vert_layout")
        self.compare_csv_button = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.compare_csv_button.setMinimumSize(QtCore.QSize(170, 24))
        self.compare_csv_button.setMaximumSize(QtCore.QSize(170, 24))
        self.compare_csv_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.compare_csv_button.setMouseTracking(True)
        self.compare_csv_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.compare_csv_button.setObjectName("compare_csv_button")
        self.csv_violinplot_vert_layout.addWidget(self.compare_csv_button, 0, QtCore.Qt.AlignHCenter)
        self.compare_csv_textbox = QtWidgets.QTextEdit(self.gridLayoutWidget_2)
        self.compare_csv_textbox.setMinimumSize(QtCore.QSize(0, 188))
        self.compare_csv_textbox.setMaximumSize(QtCore.QSize(16777215, 188))
        self.compare_csv_textbox.setStyleSheet("background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.compare_csv_textbox.setReadOnly(True)
        self.compare_csv_textbox.setObjectName("compare_csv_textbox")
        self.csv_violinplot_vert_layout.addWidget(self.compare_csv_textbox)
        self.compare_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget_2)
        self.compare_page_matplotlib_widget.setMinimumSize(QtCore.QSize(740, 400))
        self.compare_page_matplotlib_widget.setMaximumSize(QtCore.QSize(740, 400))
        self.compare_page_matplotlib_widget.setObjectName("compare_page_matplotlib_widget")
        self.csv_violinplot_vert_layout.addWidget(self.compare_page_matplotlib_widget)
        self.draw_save_hor_layout = QtWidgets.QHBoxLayout()
        self.draw_save_hor_layout.setObjectName("draw_save_hor_layout")
        self.draw_violins_button = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.draw_violins_button.setMinimumSize(QtCore.QSize(95, 24))
        self.draw_violins_button.setMaximumSize(QtCore.QSize(95, 24))
        self.draw_violins_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.draw_violins_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.draw_violins_button.setObjectName("draw_violins_button")
        self.draw_save_hor_layout.addWidget(self.draw_violins_button)
        self.save_violins_button = QtWidgets.QPushButton(self.gridLayoutWidget_2)
        self.save_violins_button.setEnabled(True)
        self.save_violins_button.setMinimumSize(QtCore.QSize(130, 24))
        self.save_violins_button.setMaximumSize(QtCore.QSize(130, 24))
        self.save_violins_button.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.save_violins_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.save_violins_button.setObjectName("save_violins_button")
        self.draw_save_hor_layout.addWidget(self.save_violins_button)
        self.csv_violinplot_vert_layout.addLayout(self.draw_save_hor_layout)
        self.compare_page_grid_layout.addLayout(self.csv_violinplot_vert_layout, 0, 1, 1, 1)
        PJGUIPages.addWidget(self.compare_page)
        self.coverage_page = QtWidgets.QWidget()
        self.coverage_page.setObjectName("coverage_page")
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.coverage_page)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(10, 10, 963, 681))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.coverage_page_grid_layout = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.coverage_page_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.coverage_page_grid_layout.setObjectName("coverage_page_grid_layout")
        self.csv_covplot_vert_layout = QtWidgets.QVBoxLayout()
        self.csv_covplot_vert_layout.setObjectName("csv_covplot_vert_layout")
        self.cov_top_hor_layout = QtWidgets.QHBoxLayout()
        self.cov_top_hor_layout.setObjectName("cov_top_hor_layout")
        self.read_hpd_csv_button = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.read_hpd_csv_button.setMinimumSize(QtCore.QSize(170, 24))
        self.read_hpd_csv_button.setMaximumSize(QtCore.QSize(170, 24))
        self.read_hpd_csv_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.read_hpd_csv_button.setObjectName("read_hpd_csv_button")
        self.cov_top_hor_layout.addWidget(self.read_hpd_csv_button, 0, QtCore.Qt.AlignHCenter)
        self.read_logfile_button = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.read_logfile_button.setMinimumSize(QtCore.QSize(150, 24))
        self.read_logfile_button.setMaximumSize(QtCore.QSize(150, 24))
        self.read_logfile_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.read_logfile_button.setObjectName("read_logfile_button")
        self.cov_top_hor_layout.addWidget(self.read_logfile_button, 0, QtCore.Qt.AlignHCenter)
        self.csv_covplot_vert_layout.addLayout(self.cov_top_hor_layout)
        self.cov_csv_coverage_grid_layout = QtWidgets.QGridLayout()
        self.cov_csv_coverage_grid_layout.setObjectName("cov_csv_coverage_grid_layout")
        self.coverage_csv_textbox = QtWidgets.QTextEdit(self.gridLayoutWidget_3)
        self.coverage_csv_textbox.setMinimumSize(QtCore.QSize(570, 170))
        self.coverage_csv_textbox.setMaximumSize(QtCore.QSize(570, 170))
        self.coverage_csv_textbox.setStyleSheet("background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"font: 11pt \"Courier\";")
        self.coverage_csv_textbox.setReadOnly(True)
        self.coverage_csv_textbox.setObjectName("coverage_csv_textbox")
        self.cov_csv_coverage_grid_layout.addWidget(self.coverage_csv_textbox, 0, 0, 1, 1, QtCore.Qt.AlignHCenter)
        self.coverage_textbox = QtWidgets.QListWidget(self.gridLayoutWidget_3)
        self.coverage_textbox.setMinimumSize(QtCore.QSize(170, 170))
        self.coverage_textbox.setMaximumSize(QtCore.QSize(170, 170))
        self.coverage_textbox.setStyleSheet("background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;")
        self.coverage_textbox.setObjectName("coverage_textbox")
        self.cov_csv_coverage_grid_layout.addWidget(self.coverage_textbox, 0, 1, 1, 1, QtCore.Qt.AlignHCenter)
        self.csv_covplot_vert_layout.addLayout(self.cov_csv_coverage_grid_layout)
        self.coverage_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget_3)
        self.coverage_page_matplotlib_widget.setMinimumSize(QtCore.QSize(740, 400))
        self.coverage_page_matplotlib_widget.setMaximumSize(QtCore.QSize(740, 400))
        self.coverage_page_matplotlib_widget.setObjectName("coverage_page_matplotlib_widget")
        self.csv_covplot_vert_layout.addWidget(self.coverage_page_matplotlib_widget, 0, QtCore.Qt.AlignHCenter)
        self.cov_draw_save_hor_layout = QtWidgets.QHBoxLayout()
        self.cov_draw_save_hor_layout.setObjectName("cov_draw_save_hor_layout")
        self.draw_cov_button = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.draw_cov_button.setMinimumSize(QtCore.QSize(95, 24))
        self.draw_cov_button.setMaximumSize(QtCore.QSize(95, 24))
        self.draw_cov_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.draw_cov_button.setObjectName("draw_cov_button")
        self.cov_draw_save_hor_layout.addWidget(self.draw_cov_button, 0, QtCore.Qt.AlignHCenter)
        self.save_cov_button = QtWidgets.QPushButton(self.gridLayoutWidget_3)
        self.save_cov_button.setMinimumSize(QtCore.QSize(130, 24))
        self.save_cov_button.setMaximumSize(QtCore.QSize(130, 24))
        self.save_cov_button.setStyleSheet("QPushButton{\n"
"    background-color: lightgray;\n"
"    border-radius: 2px;\n"
"}\n"
"QPushButton:hover{\n"
"    background-color: #f7f7f7;\n"
"}\n"
"QPushButton:pressed{\n"
"    background-color: #ffffff;\n"
"}")
        self.save_cov_button.setObjectName("save_cov_button")
        self.cov_draw_save_hor_layout.addWidget(self.save_cov_button, 0, QtCore.Qt.AlignHCenter)
        self.csv_covplot_vert_layout.addLayout(self.cov_draw_save_hor_layout)
        self.coverage_page_grid_layout.addLayout(self.csv_covplot_vert_layout, 0, 1, 1, 1)
        self.cov_node_stat_vert_layout = QtWidgets.QVBoxLayout()
        self.cov_node_stat_vert_layout.setObjectName("cov_node_stat_vert_layout")
        self.coverage_frame = QtWidgets.QFrame(self.gridLayoutWidget_3)
        self.coverage_frame.setMinimumSize(QtCore.QSize(195, 320))
        self.coverage_frame.setMaximumSize(QtCore.QSize(195, 320))
        self.coverage_frame.setStyleSheet("border-radius: 5px;\n"
"background-color: #fcf5e3;\n"
"border: 0;")
        self.coverage_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.coverage_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.coverage_frame.setObjectName("coverage_frame")
        self.coverage_node_list = QtWidgets.QListWidget(self.coverage_frame)
        self.coverage_node_list.setGeometry(QtCore.QRect(10, 30, 172, 280))
        self.coverage_node_list.setMinimumSize(QtCore.QSize(172, 280))
        self.coverage_node_list.setMaximumSize(QtCore.QSize(172, 280))
        self.coverage_node_list.setStyleSheet("QListWidget {\n"
"    background-color: #f7f5f0;\n"
"    color: #495057;\n"
"    border: 0;\n"
"    border-radius: 0px;\n"
"}\n"
"QListWidget::item:selected{\n"
"    background-color: #495057;\n"
"    color: white;\n"
"    border: 0;\n"
"}")
        self.coverage_node_list.setObjectName("coverage_node_list")
        self.coverage_node_label = QtWidgets.QLabel(self.coverage_frame)
        self.coverage_node_label.setGeometry(QtCore.QRect(20, 10, 140, 16))
        self.coverage_node_label.setMinimumSize(QtCore.QSize(140, 16))
        self.coverage_node_label.setMaximumSize(QtCore.QSize(140, 16))
        self.coverage_node_label.setAlignment(QtCore.Qt.AlignCenter)
        self.coverage_node_label.setObjectName("coverage_node_label")
        self.cov_node_stat_vert_layout.addWidget(self.coverage_frame, 0, QtCore.Qt.AlignHCenter)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.cov_node_stat_vert_layout.addItem(spacerItem5)
        self.compare_stats_label_2 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.compare_stats_label_2.setMinimumSize(QtCore.QSize(0, 16))
        self.compare_stats_label_2.setMaximumSize(QtCore.QSize(16777215, 16))
        self.compare_stats_label_2.setObjectName("compare_stats_label_2")
        self.cov_node_stat_vert_layout.addWidget(self.compare_stats_label_2, 0, QtCore.Qt.AlignHCenter)
        self.cov_summary_stats_list = QtWidgets.QListWidget(self.gridLayoutWidget_3)
        self.cov_summary_stats_list.setMinimumSize(QtCore.QSize(195, 250))
        self.cov_summary_stats_list.setMaximumSize(QtCore.QSize(195, 250))
        self.cov_summary_stats_list.setStyleSheet("QListWidget{\n"
"    background-color: #f7f7f7;\n"
"    border: 2px solid lightgray;\n"
"}\n"
"QListWidget::item:selected{\n"
"    background-color: #495057;\n"
"    color: white;\n"
"    border: 0;\n"
"}")
        self.cov_summary_stats_list.setObjectName("cov_summary_stats_list")
        self.cov_node_stat_vert_layout.addWidget(self.cov_summary_stats_list, 0, QtCore.Qt.AlignHCenter)
        self.coverage_page_grid_layout.addLayout(self.cov_node_stat_vert_layout, 0, 0, 1, 1)
        PJGUIPages.addWidget(self.coverage_page)
        self.cmd_log_page = QtWidgets.QWidget()
        self.cmd_log_page.setObjectName("cmd_log_page")
        self.verticalLayoutWidget_3 = QtWidgets.QWidget(self.cmd_log_page)
        self.verticalLayoutWidget_3.setGeometry(QtCore.QRect(10, 10, 962, 681))
        self.verticalLayoutWidget_3.setObjectName("verticalLayoutWidget_3")
        self.cmd_log_page_frame = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_3)
        self.cmd_log_page_frame.setContentsMargins(0, 0, 0, 0)
        self.cmd_log_page_frame.setObjectName("cmd_log_page_frame")
        self.cmd_log_textbox = QtWidgets.QTextEdit(self.verticalLayoutWidget_3)
        self.cmd_log_textbox.setEnabled(True)
        self.cmd_log_textbox.setMinimumSize(QtCore.QSize(960, 650))
        self.cmd_log_textbox.setMaximumSize(QtCore.QSize(16777215, 650))
        self.cmd_log_textbox.setAcceptDrops(False)
        self.cmd_log_textbox.setStyleSheet("background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;")
        self.cmd_log_textbox.setReadOnly(True)
        self.cmd_log_textbox.setObjectName("cmd_log_textbox")
        self.cmd_log_page_frame.addWidget(self.cmd_log_textbox)
        PJGUIPages.addWidget(self.cmd_log_page)
        self.warnings_page = QtWidgets.QWidget()
        self.warnings_page.setObjectName("warnings_page")
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(self.warnings_page)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(10, 10, 962, 681))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.warnings_page_frame = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.warnings_page_frame.setContentsMargins(0, 0, 0, 0)
        self.warnings_page_frame.setObjectName("warnings_page_frame")
        self.warnings_textbox = QtWidgets.QTextEdit(self.verticalLayoutWidget_2)
        self.warnings_textbox.setMinimumSize(QtCore.QSize(960, 650))
        self.warnings_textbox.setMaximumSize(QtCore.QSize(16777215, 650))
        self.warnings_textbox.setStyleSheet("background-color: #ffffff;\n"
"border: 2px solid lightgray;\n"
"border-radius: 2px;\n"
"color: red;")
        self.warnings_textbox.setReadOnly(True)
        self.warnings_textbox.setObjectName("warnings_textbox")
        self.warnings_page_frame.addWidget(self.warnings_textbox)
        PJGUIPages.addWidget(self.warnings_page)

        self.retranslateUi(PJGUIPages)
        PJGUIPages.setCurrentIndex(1)
        self.node_content_tabs.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(PJGUIPages)

    def retranslateUi(self, PJGUIPages):
        _translate = QtCore.QCoreApplication.translate
        PJGUIPages.setWindowTitle(_translate("PJGUIPages", "StackedWidget"))
        self.filename_prefix_label.setText(_translate("PJGUIPages", "File name prefix:"))
        self.random_seed_prefix_label.setText(_translate("PJGUIPages", "Random seed:"))
        self.settings_label.setText(_translate("PJGUIPages", "Output configuration"))
        self.model_label.setText(_translate("PJGUIPages", "Model nodes"))
        self.clear_model.setText(_translate("PJGUIPages", " Clear model"))
        self.one_sample_radio.setText(_translate("PJGUIPages", "One sample"))
        self.all_samples_radio.setText(_translate("PJGUIPages", "All samples"))
        self.sample_idx_spin.setPrefix(_translate("PJGUIPages", "Sample #"))
        self.repl_idx_spin.setPrefix(_translate("PJGUIPages", "Replicate #"))
        self.save_pgm_node_plot.setText(_translate("PJGUIPages", "Save plot as"))
        self.cmd_prompt_label.setText(_translate("PJGUIPages", "Command prompt"))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.values_tab), _translate("PJGUIPages", "Value(s)"))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.summary_tab), _translate("PJGUIPages", "Summary stats."))
        self.compare_node_label.setText(_translate("PJGUIPages", "Node to compare"))
        self.avg_replicate_check_button.setText(_translate("PJGUIPages", "Replicates"))
        self.compare_stats_label.setText(_translate("PJGUIPages", "Summary statistics"))
        self.compare_csv_button.setText(_translate("PJGUIPages", "Compare to .csv (...)"))
        self.draw_violins_button.setText(_translate("PJGUIPages", "Draw"))
        self.save_violins_button.setText(_translate("PJGUIPages", "Save plot as"))
        self.read_hpd_csv_button.setText(_translate("PJGUIPages", "Read HPDs .csv (...)"))
        self.read_logfile_button.setText(_translate("PJGUIPages", "Read .log file (...)"))
        self.draw_cov_button.setText(_translate("PJGUIPages", "Draw"))
        self.save_cov_button.setText(_translate("PJGUIPages", "Save plot as"))
        self.coverage_node_label.setText(_translate("PJGUIPages", "Non-det. nodes"))
        self.compare_stats_label_2.setText(_translate("PJGUIPages", "Summary statistics"))
from phylojunction.interface.pysidegui.pjguiwidgets.matplotlibwidget import MatplotlibWidget
from phylojunction.interface.pysidegui.pjguiwidgets.pj_buttons import PJClearPGMQPushButton


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    PJGUIPages = QtWidgets.QStackedWidget()
    ui = Ui_PJGUIPages()
    ui.setupUi(PJGUIPages)
    PJGUIPages.show()
    sys.exit(app.exec_())
