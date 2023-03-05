import matplotlib.pyplot as plt # type: ignore

from PySide6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure # type: ignore

class MatplotlibWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.fig = Figure(figsize=(11,4.5))
        # self.fig = Figure(figsize=(15,6))
        # self.axes = self.initialize_axes(self.fig)
        self.initialize_axes(self.fig)
        self.canvas = FigureCanvasQTAgg(self.fig)  # widget
        self.canvas.setParent(self)
        self.canvas.setMinimumSize(self.parent().size())

        fig_layout = QVBoxLayout()
        fig_layout.addWidget(self.canvas)

        self.setLayout(fig_layout)
        

    # def initialize_axes(self, fig: Figure, disabled_yticks: bool=True) -> plt.Axes:
    def initialize_axes(self, disabled_yticks: bool=True):
        # horiz coord of lower-left corner
        # vertical coord of lower-left corner
        # subplott width
        # subplot height
        # ax = fig.add_axes([0.25, 0.2, 0.5, 0.6]) 
        
        # ax = fig.add_axes([0.075, 0.25, 0.6, 0.7])
        ax = self.fig.add_axes([0.075, 0.25, 0.6, 0.7])
    
        ax.patch.set_alpha(0.0)
        ax.xaxis.set_ticks([])
        if disabled_yticks:
            ax.yaxis.set_ticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        self.axes = ax