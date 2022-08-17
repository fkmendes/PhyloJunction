import matplotlib.pyplot as plt # type: ignore
import seaborn as sns # type: ignore
import pandas as pd # type: ignore
from tabulate import tabulate
import typing as ty

def plot_violins(ax: plt.Axes, df: pd.DataFrame) -> None:
    """Draw violin plots on provided Axes object, for one variable (column 1) and two or more factors (column 2)

    Args:
        ax (plt.Axes): pyplot Axes object to draw on
        df (pd.DataFrame): pandas DataFrame with values for a variable (column 1) in two or more scenarios (i.e., two or more factors, such as different simulators; column 2)

    Returns:
        None
    """
    pass

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
