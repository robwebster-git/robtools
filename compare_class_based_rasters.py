import numpy as np
from matplotlib import pyplot as plt
import rasterio as rio
import click
import pandas as pd
import seaborn as sns

# Confusion Matrix
from sklearn.metrics import confusion_matrix


@click.option('--cal', type=click.Path(),)
@click.option('--ard', type=click.Path())
@click.command()
def main(cal, ard):
    with rio.open(ard) as src:
        ard_data = src.read(1, masked=True)
        profile = src.profile.copy()
    
    with rio.open(cal) as src:
        cal_data = src.read(1, masked=True)


    ardbeg = pd.Series(ard_data[ard_data.mask == False])
    calveg = pd.Series(cal_data[ard_data.mask == False])

    result = pd.crosstab(ardbeg, calveg, rownames=['Ardbeg Class'], colnames=['Calveg Class'])

    print(result)

    sns.heatmap(result)
    plt.show()

    #c = confusion_matrix(ard_flat, cal_flat)

    #print(c)

    #plt.imshow(c, cmap='binary', interpolation='None')
    #plt.show()

    '''
    with rio.open('masked_cal.tif', 'w', **profile) as dst:
        dst.write(masked_cal.filled(-9999).astype(np.int16), 1)
    '''



if __name__ == "__main__":
    main()
