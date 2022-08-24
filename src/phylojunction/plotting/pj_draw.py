import typing as ty
import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore
import pandas as pd # type: ignore
from matplotlib.figure import Figure # type: ignore
from tabulate import tabulate # type: ignore


def plot_violins(fig: Figure, ax: plt.Axes, df: pd.DataFrame, x: str, y: str, xlab: ty.Optional[str]=None, ylab: ty.Optional[str]=None) -> None:
    """Draw violin plots on provided Axes object, for one variable (column 1) and two or more factors (column 2)

    Args:
        ax (plt.Axes): pyplot Axes object to draw on
        df (pd.DataFrame): pandas DataFrame with values for a variable (column 1) in two or more scenarios (i.e., two or more factors, such as different simulators; column 2)
        xlab (str): x-axis label. Defaults to None.
        ylab (str): y-axis label. Defaults to None.
        color1 (str): Color for first violin plot. Defaults to None.
        color2 (str): Color for second violing plot. Defaults to None.

    Returns:
        None
    """
    
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

    custom_palette = sns.color_palette("Spectral")
    sns.set_palette(custom_palette)
    sns.violinplot(x=x, y=y, data=df, ax=ax, linewidth=0.0)
    
    if xlab != None:
        ax.set_xlabel(xlab)

    if ylab != None:
        ax.set_ylabel(ylab)

    fig.canvas.draw()
    

if __name__ == "__main__":
    # testing below
    
    df1 = pd.DataFrame({
        "Total taxon count": [1, 3, 4, 5, 5, 5, 5, 6, 7, 10],
        "Program": ["PJ" for i in range(10)]
        })

    df2 = pd.DataFrame({
        "Total taxon count": [1, 3, 4, 5, 5, 5, 5, 6, 7, 10],
        "Program": ["Other" for i in range(10)]
        })
    
    df = pd.concat([df1, df2], axis=0)

    print(tabulate(df, headers=["", "Total taxon count", "Program"], tablefmt="pretty"))

    fig, axes = plt.subplots()
    # plot violin. 'Scenario' is according to x axis, 
    # 'LMP' is y axis, data is your dataframe. ax - is axes instance
    sns.violinplot("Program", "Total taxon count", data=df, ax=axes)
    # axes.set_title("")
    # axes.yaxis.grid(True)
    # axes.set_xlabel("Simulator")
    # axes.set_ylabel("Total taxon count")

    plt.show()
