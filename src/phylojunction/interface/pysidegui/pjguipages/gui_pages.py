# Form implementation generated from reading ui file 'pjgui_pages.ui'
#
# Created by: PyQt6 UI code generator 6.3.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PySide6 import QtCore, QtGui, QtWidgets
from phylojunction.interface.pysidegui.pjguiwidgets.matplotlibwidget import MatplotlibWidget


class Ui_PJGUIPages(object):
    def setupUi(self, PJGUIPages):
        PJGUIPages.setObjectName("PJGUIPages")
        PJGUIPages.resize(972, 726)
        self.pgm_page = QtWidgets.QWidget()
        self.pgm_page.setObjectName("pgm_page")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.pgm_page)
        self.verticalLayout.setObjectName("verticalLayout")
        self.pgm_page_frame = QtWidgets.QFrame(self.pgm_page)
        self.pgm_page_frame.setMinimumSize(QtCore.QSize(950, 700))
        self.pgm_page_frame.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.pgm_page_frame.setStyleSheet("background-color: white;\n"
"border: 0;")
        self.pgm_page_frame.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.pgm_page_frame.setFrameShadow(QtWidgets.QFrame.Shadow.Raised)
        self.pgm_page_frame.setObjectName("pgm_page_frame")
        self.gridLayoutWidget = QtWidgets.QWidget(self.pgm_page_frame)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(0, 180, 951, 451))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.pgm_fig_node_list_grid_layout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.pgm_fig_node_list_grid_layout.setContentsMargins(0, 0, 0, 0)
        self.pgm_fig_node_list_grid_layout.setObjectName("pgm_fig_node_list_grid_layout")
        self.node_list = QtWidgets.QListWidget(self.gridLayoutWidget)
        self.node_list.setMinimumSize(QtCore.QSize(200, 400))
        self.node_list.setMaximumSize(QtCore.QSize(200, 16777215))
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
        self.pgm_fig_node_list_grid_layout.addWidget(self.node_list, 1, 1, 1, 1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.one_sample_radio = QtWidgets.QRadioButton(self.gridLayoutWidget)
        self.one_sample_radio.setStyleSheet("color: black;\n"
"")
        self.one_sample_radio.setCheckable(False)
        self.one_sample_radio.setChecked(False)
        self.one_sample_radio.setObjectName("one_sample_radio")
        self.horizontalLayout.addWidget(self.one_sample_radio)
        self.all_samples_radio = QtWidgets.QRadioButton(self.gridLayoutWidget)
        self.all_samples_radio.setStyleSheet("color: black;")
        self.all_samples_radio.setCheckable(False)
        self.all_samples_radio.setAutoExclusive(True)
        self.all_samples_radio.setObjectName("all_samples_radio")
        self.horizontalLayout.addWidget(self.all_samples_radio)
        spacerItem = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.sample_idx_spin = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.sample_idx_spin.setStyleSheet("color: black;")
        self.sample_idx_spin.setAccelerated(True)
        self.sample_idx_spin.setObjectName("sample_idx_spin")
        self.horizontalLayout.addWidget(self.sample_idx_spin)
        spacerItem1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Policy.Maximum, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.repl_idx_spin = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.repl_idx_spin.setStyleSheet("color: black;\n"
"")
        self.repl_idx_spin.setReadOnly(False)
        self.repl_idx_spin.setAccelerated(True)
        self.repl_idx_spin.setObjectName("repl_idx_spin")
        self.horizontalLayout.addWidget(self.repl_idx_spin)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.pgm_fig_node_list_grid_layout.addLayout(self.horizontalLayout, 0, 0, 1, 1)
        self.model_label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.model_label.setMinimumSize(QtCore.QSize(200, 25))
        self.model_label.setMaximumSize(QtCore.QSize(200, 20))
        self.model_label.setStyleSheet("font: 14pt \"Ubuntu\";\n"
"color: black;\n"
"")
        self.model_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.model_label.setObjectName("model_label")
        self.pgm_fig_node_list_grid_layout.addWidget(self.model_label, 0, 1, 1, 1)
        self.pgm_page_matplotlib_widget = MatplotlibWidget(self.gridLayoutWidget)
        self.pgm_page_matplotlib_widget.setMinimumSize(QtCore.QSize(740, 400))
        self.pgm_page_matplotlib_widget.setMaximumSize(QtCore.QSize(16777215, 400))
        self.pgm_page_matplotlib_widget.setObjectName("pgm_page_matplotlib_widget")
        self.pgm_fig_node_list_grid_layout.addWidget(self.pgm_page_matplotlib_widget, 1, 0, 1, 1)
        self.verticalLayoutWidget = QtWidgets.QWidget(self.pgm_page_frame)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 630, 951, 63))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.cmd_prompt_vert_layout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.cmd_prompt_vert_layout.setContentsMargins(0, 0, 0, 0)
        self.cmd_prompt_vert_layout.setObjectName("cmd_prompt_vert_layout")
        self.cmd_prompt_label = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.cmd_prompt_label.setMinimumSize(QtCore.QSize(940, 20))
        self.cmd_prompt_label.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.cmd_prompt_label.setStyleSheet("font: 14pt \"Ubuntu\";\n"
"color: black;")
        self.cmd_prompt_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignBottom|QtCore.Qt.AlignmentFlag.AlignLeading|QtCore.Qt.AlignmentFlag.AlignLeft)
        self.cmd_prompt_label.setObjectName("cmd_prompt_label")
        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt_label)
        self.cmd_prompt = QtWidgets.QLineEdit(self.verticalLayoutWidget)
        self.cmd_prompt.setMinimumSize(QtCore.QSize(940, 35))
        self.cmd_prompt.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.cmd_prompt.setToolTip("")
        self.cmd_prompt.setStyleSheet("background-color: black;\n"
"font: 14pt \"Courier\";\n"
"font-weight: bold;\n"
"color: white;\n"
"border-radius: 5px;")
        self.cmd_prompt.setCursorPosition(0)
        self.cmd_prompt.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeading|QtCore.Qt.AlignmentFlag.AlignLeft|QtCore.Qt.AlignmentFlag.AlignVCenter)
        self.cmd_prompt.setDragEnabled(True)
        self.cmd_prompt.setPlaceholderText("")
        self.cmd_prompt.setObjectName("cmd_prompt")
        self.cmd_prompt_vert_layout.addWidget(self.cmd_prompt)
        self.node_content_tabs = QtWidgets.QTabWidget(self.pgm_page_frame)
        self.node_content_tabs.setGeometry(QtCore.QRect(0, 0, 951, 181))
        self.node_content_tabs.setMinimumSize(QtCore.QSize(940, 0))
        self.node_content_tabs.setStyleSheet("color: black;\n"
"font: 14pt ;\n"
"")
        self.node_content_tabs.setTabShape(QtWidgets.QTabWidget.TabShape.Rounded)
        self.node_content_tabs.setIconSize(QtCore.QSize(16, 16))
        self.node_content_tabs.setObjectName("node_content_tabs")
        self.values_tab = QtWidgets.QWidget()
        self.values_tab.setObjectName("values_tab")
        self.values_content = QtWidgets.QTextEdit(self.values_tab)
        self.values_content.setGeometry(QtCore.QRect(0, 0, 951, 141))
        self.values_content.setMinimumSize(QtCore.QSize(940, 0))
        self.values_content.setStyleSheet("border: 4px solid lightgray;\n"
"border-radius: 4px;\n"
"font: 11pt \"Courier\";")
        self.values_content.setReadOnly(True)
        self.values_content.setObjectName("values_content")
        self.values_content.raise_()
        self.gridLayoutWidget.raise_()
        self.node_content_tabs.addTab(self.values_tab, "")
        self.summary_tab = QtWidgets.QWidget()
        self.summary_tab.setObjectName("summary_tab")
        self.summary_content = QtWidgets.QTextEdit(self.summary_tab)
        self.summary_content.setGeometry(QtCore.QRect(0, 0, 951, 141))
        self.summary_content.setMinimumSize(QtCore.QSize(940, 0))
        self.summary_content.setStyleSheet("border: 4px solid lightgray;\n"
"border-radius: 4px;\n"
"font: 11pt \"Courier\";")
        self.summary_content.setObjectName("summary_content")
        self.node_content_tabs.addTab(self.summary_tab, "")
        self.verticalLayout.addWidget(self.pgm_page_frame)
        PJGUIPages.addWidget(self.pgm_page)
        self.cmd_log_page = QtWidgets.QWidget()
        self.cmd_log_page.setObjectName("cmd_log_page")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.cmd_log_page)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(self.cmd_log_page)
        self.label_2.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        PJGUIPages.addWidget(self.cmd_log_page)

        self.retranslateUi(PJGUIPages)
        PJGUIPages.setCurrentIndex(0)
        self.node_content_tabs.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(PJGUIPages)

    def retranslateUi(self, PJGUIPages):
        _translate = QtCore.QCoreApplication.translate
        PJGUIPages.setWindowTitle(_translate("PJGUIPages", "StackedWidget"))
        self.one_sample_radio.setText(_translate("PJGUIPages", "One sample"))
        self.all_samples_radio.setText(_translate("PJGUIPages", "All samples"))
        self.sample_idx_spin.setPrefix(_translate("PJGUIPages", "Sample #"))
        self.repl_idx_spin.setPrefix(_translate("PJGUIPages", "Replicate #"))
        self.model_label.setText(_translate("PJGUIPages", "Model nodes"))
        self.cmd_prompt_label.setText(_translate("PJGUIPages", "Command prompt"))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.values_tab), _translate("PJGUIPages", "Value(s)"))
        self.node_content_tabs.setTabText(self.node_content_tabs.indexOf(self.summary_tab), _translate("PJGUIPages", "Summary stats."))
        self.label_2.setText(_translate("PJGUIPages", "Page 2"))